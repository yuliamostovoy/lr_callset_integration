# PAV Whole-Callset bcftools merge — Run Order

A standalone pipeline that cohort-merges **raw PAV calls including SNVs** (SNV +
indel + SV) for a very large cohort (~13k samples) using **only `bcftools
merge`** — no kanpig, scoring, or truvari collapse. It borrows the scaling
tricks from the `SV_Integration_Workpackage*` pipeline (genome chunking,
two-level hierarchical merge, `SumFileSizes` disk pre-check, preemptible VMs)
but is a separate set of files.

## Files

1. `PAV_BcftoolsMerge_1_Split.wdl` — per-sample normalize + chunk
2. `PAV_BcftoolsMerge_2_Merge.wdl` — per-chunk cohort `bcftools merge`
3. `PAV_BcftoolsMerge_3_Concat.wdl` — chunks → one BCF per primary contig
4. `PAV_BcftoolsMerge_Orchestrator.wdl` — runs 1→2→3 as one submission
5. `../hg38_split_5mb_for_bcftools_merge.csv` — default ~5 Mbp chunk partition

Run stages 1→2→3 in order (or run the orchestrator). Each stage is idempotent
(`.done` markers), so re-submitting only redoes failed units.

## No hand-built TSVs between stages

Two sources of truth drive everything:

- **The Terra data table** = samples (`sample_id` + PAV VCF URI). Consumed only
  by stage 1, via `this.samples.sample_id` / `this.samples.pav_vcf` as
  **`Array[String]`** (URIs are streamed in-task, not localized up front — this
  is what preserves the scaling method).
- **The split CSV** = chunks. Chunk id == 0-based CSV line number. Stage 2
  derives the chunk count from it (`scatter` over `range(n_chunks)`); stage 3
  derives the chunk→chromosome mapping from it.

Stage 1 auto-writes `sample_ids.txt` (canonical merge order) to its output dir;
stages 2/3 pick it up from GCS. So you configure stage 1 once against the data
table, and stages 2/3 only need *(the same CSV + the previous stage's GCS dir)*.

## What each stage does

- **Stage 1 (Split).** Root entity `sample_set`. Batches the cohort in-WDL
  (`batch_size` samples per VM). Per sample: localize → reheader to `sample_id`
  → **`bcftools norm -f ref -m -any`** (left-align + split multiallelics) → split
  into chunks with `bcftools view --targets-file` (POS-only) → upload
  `chunk_<i>/<sample>.bcf`.

  Normalization runs **before** chunking so every record is at its final
  canonical position before chunk assignment; otherwise an indel could
  left-align across a chunk boundary and fragment the same variant across
  chunks. Every sample emits a (possibly empty) BCF for every chunk, because the
  merge and concat require an identical sample set/order in every chunk.

- **Stage 2 (Merge).** One `MergeChunk` VM per chunk: verify #files == #samples,
  `SumFileSizes` disk pre-check, localize, **two-level `bcftools merge --merge
  none`** (batches of `n_files_per_merge`, then merge-of-merges), then `bcftools
  norm --do-not-normalize -m -any` to guarantee biallelic-only output. Uploads
  `chunk_<i>.bcf`.

- **Stage 3 (Concat).** `bcftools concat --naive` the chunks of each primary
  contig (chr1-22, X, Y; in CSV/position order) into one `<contig>/<contig>.bcf`.
  There is deliberately no genome-wide file: at 13k+SNV scale it would be ~1 TB.
  Per-contig BCFs are the final deliverable; tools wanting finer granularity can
  consume stage 2's per-chunk BCFs directly.

## Phasing (critical)

PAV per-haplotype phasing is preserved end to end. `bcftools norm -m -any` only
ever **splits** multiallelics (never joins), preserving the `|` separator,
haplotype assignment, and PS tags through both multiallelic splitting and indel
left-alignment. `bcftools merge --merge none` copies each sample's GT verbatim.
The stage-2 post-merge split uses `--do-not-normalize` (no left-alignment), so
positions never move out of a chunk. **Never** switch any `norm` to `-m +any`
(join) or the merge to `-m both`.

## Scaling notes

- The default ~5 Mbp CSV (`hg38_split_5mb_for_bcftools_merge.csv`, 632 chunks)
  is finer than the SV pipeline's 30 Mbp CSV because PAV-with-SNVs is ~2 orders
  of magnitude heavier per bp. Regenerate for a different size with
  `java SplitForBcftoolsMerge <ref.fai> <chunk_bp>`.
- **Pilot before the full fan-out:** run stage 2 on one dense chunk (e.g. an
  HLA/chr6 window) and one sparse chunk, read the emitted `df -h` /
  `time --verbose` output, then set `merge_disk_size_gb` / `merge_ram_size_gb`
  for the worst (densest) chunk. The `SumFileSizes` pre-check fails a chunk
  cheaply if inputs won't fit, rather than hanging a VM.
