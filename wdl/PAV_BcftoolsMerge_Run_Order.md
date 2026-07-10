# PAV Whole-Callset bcftools merge — Run Order

A standalone pipeline that cohort-merges **raw PAV calls including SNVs** (SNV +
indel + SV) for a very large cohort (~13k samples) using **only `bcftools
merge`** — no kanpig, scoring, or truvari collapse. It borrows the scaling
tricks from the `SV_Integration_Workpackage*` pipeline (genome chunking,
two-level hierarchical merge, `SumFileSizes` disk pre-check, preemptible VMs)
but is a separate set of files.

## Files

1. `PAV_BcftoolsMerge_1_Split.wdl` — per-sample normalize + chunk. Run once
   **per center** (see "Multiple centers" below).
2. `PAV_BcftoolsMerge_2_Merge.wdl` — per-chunk cohort `bcftools merge` across
   all enabled centers, with overlap resolution.
3. `PAV_BcftoolsMerge_3_Concat.wdl` — chunks → one BCF per primary contig
4. `PAV_BcftoolsMerge_Orchestrator.wdl` — runs stage 2→3 as one submission,
   picking up from already-completed per-center stage 1 runs.
5. `../hg38_split_5mb_for_bcftools_merge.csv` — default ~5 Mbp chunk partition

Run stage 1 once per enabled center, then stage 2→3 in order (or run the
orchestrator for 2→3). Each stage is idempotent (`.done` markers), so
re-submitting only redoes failed units.

## Multiple centers, and why stage 1 is not in the orchestrator

The cohort is assembled from up to **six sources**: `bi`, `ha`, `bcm`, `uw`,
and two control sets, `controls_15x` / `controls_30x`. Each source has its own
Terra data table. A single Terra submission can only iterate one data table at
a time (`this.samples.sample_id` resolves against whichever table is the
submission's root entity), so **stage 1 is run separately, once per enabled
source**, submitted directly against that source's own table, each pointed at
its own `remote_outdir` (e.g. `.../01_split_bi`, `.../01_split_ha`, ...). A
source with no samples is simply not run.

This is why the orchestrator only covers stages 2→3: there is no single moment
where all sources' stage-1 work can be kicked off together, since each is
submitted against a different table on its own schedule. Once every enabled
source's stage 1 has completed, point the orchestrator (or stage 2 directly)
at all six `remote_outdir`s (unused sources get `n_expected_samples_<source>
= 0`, and their `remote_indir_<source>` value is then ignored) to run the
merge and concat.

**Overlap resolution.** A sample present in more than one source's chunk dir
is resolved by *localization order* in stage 2, not by any content-based
rule: `bi`, then `ha`, then `uw`, then `bcm`, then the controls — each later
source silently overwrites an earlier source's file for the same `sample_id`.
So by default **`ha` wins over `bi`**, and **each later source wins over every
earlier one**. The one exception is `bi_samples_to_prefer_over_ha`, an
`Array[String]` of `sample_id`s for which `bi`'s copy is re-localized *after*
`ha`, flipping the winner back to `bi` for exactly those samples. There is no
equivalent override for `uw`/`bcm`/the controls — for those, whichever source
is localized last always wins on overlap. This exactly mirrors
`SV_Integration_Workpackage5`'s `InterCenterMerge` pattern.

## No hand-built TSVs between stages

Two sources of truth drive everything:

- **The Terra data table** = samples (`sample_id` + PAV VCF URI), one per
  center/control set. Consumed only by stage 1, via `this.samples.sample_id` /
  `this.samples.pav_vcf` as **`Array[String]`** (URIs are streamed in-task,
  not localized up front — this is what preserves the scaling method).
- **The split CSV** = chunks. Chunk id == 0-based CSV line number. Stage 2
  derives the chunk count from it (`scatter` over `range(n_chunks)`); stage 3
  derives the chunk→chromosome mapping from it.

Stage 1 auto-writes `sample_ids.txt` (its source's own sample list) to its
output dir. Stage 2 unions every enabled source's `sample_ids.txt` into the
canonical, deduplicated merge-column order automatically if
`sample_ids_file` is omitted. So you configure stage 1 once per source
against that source's data table, and stage 2/3 only need *(the same CSV +
every enabled source's GCS dir + expected sample counts)*.

## What each stage does

- **Stage 1 (Split).** Root entity `sample_set`. Run once per enabled center/
  control set, submitted against that source's own Terra table, each writing
  to its own `remote_outdir`. Batches the cohort in-WDL
  (`batch_size` samples per VM), then runs two SEPARATE tasks per batch:
  - **`NormalizeBatch`** (expensive, CSV-independent): per sample localize →
    reheader to `sample_id` → **`bcftools norm -f ref -m -any`** (left-align +
    split multiallelics) → upload `<norm_remote_dir>/<sample>.bcf`. Skips any
    sample already present in `norm_remote_dir`.
  - **`ChunkBatch`** (cheap, CSV-dependent): localize the normalized BCF → slice
    with `bcftools view --targets-file` (POS-only) → upload
    `chunk_<i>/<sample>.bcf`.

  Normalization runs **before** chunking so every record is at its final
  canonical position before chunk assignment; otherwise an indel could
  left-align across a chunk boundary and fragment the same variant across
  chunks. Every sample emits a (possibly empty) BCF for every chunk, because the
  merge and concat require an identical sample set/order in every chunk.

  **Reuse:** normalization is the CPU-heavy step and does not depend on the
  chunk partition, so it is a separate task writing to a **stable
  `norm_remote_dir`**. Point every run (pilot and full, and any chunk-size
  retune) at the SAME `norm_remote_dir`; `NormalizeBatch` existence-skips samples
  already there, so only the cheap `ChunkBatch` re-runs. This is more reliable
  than Cromwell call caching here: because these tasks upload via `gcloud`
  inside the command, the CSV content and output path are in the cache key, so a
  changed CSV or output dir misses — the existence-check does not.

- **Stage 2 (Merge).** One `MergeChunk` VM per chunk: for each enabled source,
  verify #files == `n_expected_samples_<source>` (controls only require `>=`,
  warning if over), `SumFileSizes` disk pre-check across all sources, then
  localize every source's chunk dir into the same local directory **in the
  fixed order bi → ha → (bi_samples_to_prefer_over_ha override) → uw → bcm →
  controls_15x → controls_30x** — this is the overlap-resolution step
  described above. Then **two-level `bcftools merge --merge none`** (batches
  of `n_files_per_merge`, then merge-of-merges), then `bcftools norm
  --do-not-normalize -m -any` to guarantee biallelic-only output. Uploads
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
