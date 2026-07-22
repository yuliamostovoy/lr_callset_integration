# SV Integration — Consolidated Run Order (4 workflows)

This is the table-driven consolidation of the old 12-submission SV integration
run (see `SV_Integration_Run_Order.md` for the original per-workpackage chain).
It merges the workpackages into **4 discrete Terra submissions** and removes every
hand-built inter-step file. It covers the main short-SV chain plus the non-scored
ultralong/BND cohort integration; it deliberately excludes the ultralong/BND
annotation + XGBoost scoring branches.

## The 4 workflows

| Step | Workflow | Merges | Entity / driver |
|------|----------|--------|-----------------|
| A | `SV_Integration_WorkflowA_Intrasample_Scoring.wdl` | WP1 + WP2 | `sample` table (whole cohort) |
| B | `SV_Integration_WorkflowB_Merge_Collapse.wdl` | WP3 + WP4 + WP5 + WP6 | GCS dirs + interval CSV |
| C | `SV_Integration_WorkflowC_Regenotype.wdl` | WP7 + (WP3 merge) + WP8 | one `sample_set` per instance |
| D | `SV_Integration_WorkflowD_Ultralong_Bnd.wdl` | WP12 + WP13 + WP14 + WP15 | GCS dirs + interval CSV |

Dependency order: **A → { B ∥ D } → C**. B and D both branch off A's outputs and
can run at the same time; C runs after B.

`SV_Integration_WorkflowD_PerSuffix.wdl` is a sub-workflow of D (one call per
suffix); it is not launched directly.

## No hand-built TSVs between steps

Two sources of truth drive everything, exactly like the PAV pipeline:

- **The Terra `sample` table** — one row per sample, with the per-sample input
  URIs. Consumed by Step A (as parallel `Array[String]` columns) and, via
  PED-derived **sample sets**, by Step C.
- **The interval CSV** (`hg38_split_for_bcftools_merge.csv`) — chunk id == 0-based
  line number. Steps B, C, D derive chunk counts / chunk-id lists from it in-graph.

Everything the old run built by hand is now produced inside the WDLs:

| Old manual artifact | Now |
|---------------------|-----|
| subset the table into `sv_integration_chunk_tsv` (x2, for WP1 and WP2) | `MakeManifests` in A, from `this.samples.*` |
| cohort `sample_ids` list (WP3, WP12) | `WriteSampleList` / `WriteSampleListSuffix` (list the upstream dir's `.done` / `*_<suffix>.bcf`) |
| `bcftools_chunks` comma string (WP4) | `ChunksForChromosome` (awk over the CSV) |
| `chunks_ids` files via `make_workpackage7_chunk_id_files.sh` (WP5, WP14) | `DeriveChunkIds` (reads WP4/WP13 `regions.txt` col 2) |
| family arrays + `make_terra_family_set_membership.sh` + ped | Step C binds `this.samples.*` per sample set; `MakeFamilyInputs` synthesizes the ped + merge order |
| `chunk_ids` string (WP8) | `MakeFamilyInputs` derives the genome-ordered list from the CSV |
| copy each step's `remote_outdir` into the next `remote_indir` | one submission per step; stage dirs derived from `remote_outdir`, ordered by `done`/`upstream_signal` handshakes |
| `n_expected_samples` (WP12) | auto-derived by `WriteSampleListSuffix` |

The reused WP tasks gained only additive, backward-compatible edits: a
`String done` output on producer tasks, an optional `Array[String]? upstream_signal`
(or `String upstream_signal` on WP2) on consumer tasks, a `File regions_txt`
output on WP4/WP13, and a `.done` idempotency skip on WP3. Standalone runs of the
original workpackages are unaffected.

## Step A — WorkflowA_Intrasample_Scoring (per-sample)

Root entity: `sample`. Launch on the whole `sample` table (or a `sample_set`
covering the cohort). Bind these columns (names are whatever your table uses):

| Input | Bind to |
|-------|---------|
| `sample_ids` | `this.sample_id` |
| `sample_sexes` | `this.sex` (must exist; 'M' ⇒ male ploidy) |
| `aligned_bais` / `aligned_bams` | `this.02_aligned_bai` / `this.02_aligned_bam` |
| `pbsv_tbis` / `pbsv_vcfs` | `this.03_pbsv_tbi` / `this.03_pbsv_vcf` |
| `sniffles_tbis` / `sniffles_vcfs` | `this.03_sniffles_tbi` / `this.03_sniffles_vcf` |
| `pav_beds` / `pav_tbis` / `pav_vcfs` | `this.03_pav_bed` / `this.03_pav_tbi` / `this.03_pav_vcf` (only if `has_pav=true`) |

Set one `remote_outdir`; A writes WP1 per-sample outputs to `<remote_outdir>/01_intrasample`
(consumed by D) and WP2 scored chunks to `<remote_outdir>/02_scoring` (consumed by
B). `MakeManifests` rebuilds the exact per-sample TSV the containers expect — its
paste order is fixed in the WDL and matches the container `cut` positions; you only
pick which table column feeds each input.

## Step B — WorkflowB_Merge_Collapse (cohort, no table)

`remote_indir` = A's `<remote_outdir>/02_scoring`; `remote_outdir` = a fresh dir
(stages land in `/03_merge`, `/04_shard`, `/05_collapse`, `/06_concat`). The genome-wide
callset is `<remote_outdir>/06_concat/truvari_collapsed.bcf`. `merge_mode=1` for the
main chain. `sample_ids_file` is optional (auto-derived from A's `.done` markers).

## Step D — WorkflowD_Ultralong_Bnd (cohort, no table)

`remote_indir` = A's `<remote_outdir>/01_intrasample` (holds `<sample>_ultralong.bcf`
and `<sample>_bnd.bcf`). Runs both suffixes in one submission; per-suffix callset is
`<remote_outdir>/<suffix>/15_concat/truvari_collapsed.bcf`. No annotation, no scoring.

## Step C — WorkflowC_Regenotype (per family, sample sets)

Build one Terra `sample_set` per family/group to joint-genotype (from your PED),
then launch Step C on the selected sets — Terra runs one workflow copy per set.
Bind `family_id`=`this.sample_set_id`, and `this.samples.sample_id / sample_sex /
aligned_bai / aligned_bam`. `remote_indir` = B's remote_outdir (C reads its
`/06_concat` subdir automatically). The final regenotyped BCF is
`<remote_outdir>/08_concat/merged.bcf`. The WP3 per-chunk merge
(merge_mode=2) between WP7 and WP8 is wired in automatically.

## Idempotency

Every heavy step keeps its GCS `.done` markers (and WP3 now has one too), so
resubmitting a step only reruns the units that failed.
