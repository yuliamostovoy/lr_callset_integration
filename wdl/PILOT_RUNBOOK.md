# SV Integration consolidated pipeline — Pilot Runbook

Cheapest-first pilots for the consolidated workflows (A → { B ∥ A2 } → D → C),
before scaling to a full cohort. See `SV_Integration_Consolidated_Run_Order.md`
for the design. Each step is chained to the previous one only by a GCS dir.

Concrete values below are from the PacBio pilot in workspace
`talk-LR-gnomADLR_supplement/Talk_LISA`; swap paths for your own run.

## Prerequisites

- Workflows imported into Terra (Dockstore sync from repo `main`, or direct WDL).
  Workflow D also needs `SV_Integration_WorkflowD_PerSuffix.wdl` present.
- Pilot interval CSV uploaded (5-line test set already in repo root:
  `hg38_test_chunks_for_bcftools_merge.csv`). It covers **chr1, chr2, chr6, chr19,
  chrM**.
- One pilot output root, e.g. `gs://<bucket>/…/pilot`.
- Reuse the same file inputs across steps: `reference_fa`/`reference_fai`,
  `reference_agp`, `standard_chromosomes_bed`, `autosomes_bed`, `ploidy_bed_*`,
  `training_resource_*`, `training_python_script`, `scoring_python_script`,
  `hyperparameters_json`.

**Chromosome-list gotcha:** in Steps B and D set the chromosome list to exactly
the chromosomes present in your CSV (start with `chr1, chr2, chr6, chr19`; drop
`chrM`). WP6/WP15 hard-error on a chromosome that produced zero collapse chunks.

**ONT disk gotcha (Step A):** WP1 localizes each sample's full BAM, and its disk
default (`Intrasample.disk_size_gb = 256`) is sized for HiFi. ONT whole-genome
BAMs are larger (this cohort ranged up to ~391 GiB), so samples with big BAMs
fail during BAM download with "job stopped before the command finished" (VM
killed on a full disk — no `rc`/`monitoring.log` written). For an ONT cohort,
override `SV_Integration_WorkflowA_Intrasample_Scoring.Intrasample.disk_size_gb`
to comfortably clear the largest BAM (~500 for ~400 GiB BAMs; check with
`gsutil ls -l <bam-glob>`). WP1 is idempotent per sample (`<sample>.done`
markers), so resubmitting to the same `remote_outdir` with the larger disk only
recomputes the samples that hadn't finished.

**ONT memory gotcha (Step A):** WP1's Kanpig step also runs hotter on ONT. The
memory default (`Intrasample.ram_size_gb = 8`, HiFi-tuned) OOM-kills Kanpig on
dense ONT samples — signature: shard `rc = 137` (SIGKILL), `monitoring.log` shows
memory pinned at the VM total, and stderr ends mid-Kanpig with `Command
terminated by signal 9`. For ONT, override
`SV_Integration_WorkflowA_Intrasample_Scoring.Intrasample.ram_size_gb` to ~24.
Kanpig memory scales with variant/chunk count × threads, so if the largest
samples still OOM, bump further or lower `Intrasample.n_cpu` (fewer Kanpig
threads = fewer concurrent variant graphs = less peak memory). Same per-sample
idempotency applies on resubmit.

---

## Pilot A — `SV_Integration_WorkflowA_Intrasample_Scoring` (per-sample)

Make a small `sample_set` (2–3 samples) and run **with inputs defined by data
table**, root entity = that set (array inputs require a set root).

Table-column bindings (right-hand side = your table's column names):

| Input | Expression |
|---|---|
| `sample_ids` | `this.samples.sample_id` |
| `sample_sexes` | `this.samples.sex` |
| `aligned_bais` / `aligned_bams` | `this.samples.02_aligned_bai` / `this.samples.02_aligned_bam` |
| `pbsv_tbis` / `pbsv_vcfs` | `this.samples.03_pbsv_tbi` / `this.samples.03_pbsv_vcf` |
| `sniffles_tbis` / `sniffles_vcfs` | `this.samples.03_sniffles_tbi` / `this.samples.03_sniffles_vcf` |
| `pav_beds` / `pav_tbis` / `pav_vcfs` | `this.samples.03_pav_*` (only if `has_pav=true`) |

Literals:
- `has_pav` = `true`/`false` (PacBio pilot used `false` → leave the 3 pav_* unbound)
- `batch_size` = `1` for the pilot
- `remote_outdir` = `gs://<bucket>/…/pilot/A`  → A writes `/A/01_intrasample`
  (consumed by D) and `/A/02_scoring` (consumed by B)
- `split_for_bcftools_merge_csv` = the pilot CSV
- `requester_pays_project` if BAMs are in a requester-pays bucket (pilot used
  `longreadassembly`)
- plus the shared file inputs

**Check:** `gsutil ls gs://<bucket>/…/pilot/A/02_scoring/` shows `<sample>.done`
per sample and `chunk_0/`…`chunk_N/` each with one `<sample>.bcf` per sample.

---

## Pilot B — `SV_Integration_WorkflowB_Merge_Collapse` (cohort, no table)

Run **with inputs defined by file paths**. All literals:

- `remote_indir` = `gs://<bucket>/…/pilot/A`  (A's remote_outdir; B reads its `/02_scoring` subdir automatically)
- `remote_outdir` = `gs://<bucket>/…/pilot/B`
- `split_for_bcftools_merge_csv` = the pilot CSV
- `chromosomes` = `["chr1","chr2","chr6","chr19"]`
- `merge_mode` = `1`
- `sample_ids_file` = **unset** (auto-derived from A's `.done` markers)
- rest default

**Check:** `gsutil ls gs://<bucket>/…/pilot/B/06_concat/truvari_collapsed.bcf`.
Confirm the auto-derived sample list matches your cohort and the record count is sane.

---

## Pilot A2 — `SV_Integration_WorkflowA2_UltralongBnd_AnnotateScore` (per-sample; between A and D)

Per-sample ultralong/BND **annotation + XGBoost scoring** (upstream pretrained,
coverage-selected models), **score-only — no variant filtering** (left to the
user, like the short-SV branch). Run **with inputs defined by data table**, root
entity = a `sample_set` (array inputs require a set root), same as Pilot A.

Table-column bindings (right-hand side = your table's column names):

| Input | Expression |
|---|---|
| `sample_ids` | `this.samples.sample_id` |
| `mean_coverages` | `this.samples.mean_coverage` |
| `aligned_bais` / `aligned_bams` | `this.samples.02_aligned_bai` / `this.samples.02_aligned_bam` |

Literals:
- `remote_indir` = `gs://<bucket>/…/pilot/A`  (A's remote_outdir; A2 reads its `/01_intrasample` subdir automatically)
- `remote_outdir` = `gs://<bucket>/…/pilot/A2`  → A2 writes `/01a_annotated` (intermediate) and `/01b_scored` (consumed by D)
- `batch_size` = `1` for the pilot
- `requester_pays_project` = your billing project **if the BAMs are in a requester-pays bucket** (the PacBio pilot used `longreadassembly`); else leave blank
- Annotate resources: `reference_fa` / `reference_fai`, `feature_extraction_py`, `tr_bed`, `segdup_bed`, `gc_content_bed`, `mei_bed_gz` / `mei_bed_tbi`
- Score resources: the **12 pretrained model bundles** `{del,ins,dup,insdup,inv,bnd}_indel_{scorer_15x_pkl,scorer_30x_pkl,calibrationScores_15x_hdf5,calibrationScores_30x_hdf5}`, `scoring_python_script` (= `scripts/xgb.py`), `UltralongInsdups2Ins_java` / `AddSvlenToSymbolicAlt_java` (from `scripts/`)
- `mean_coverages` feeds BOTH the annotate coverage-bin normalization and the 15x-vs-30x model pick (threshold 22.5×)
- rest default (annotation feature lists are baked into the score step, matched to the models)

**Check:** `gsutil ls gs://<bucket>/…/pilot/A2/01b_scored/` shows `<sample>_ultralong.bcf`,
`<sample>_bnd.bcf` (+`.csi`) and `<sample>.done` per sample; records carry
`FORMAT/SCORE` and `FORMAT/CALIBRATION_SENSITIVITY`.

---

## Pilot D — `SV_Integration_WorkflowD_Ultralong_Bnd` (cohort, no table; after A2)

Run **with inputs defined by file paths**:

- `remote_indir` = `gs://<bucket>/…/pilot/A2`  (A2's remote_outdir)
- `intrasample_subdir` = `"01b_scored"`  (so D reads the **scored** per-sample BCFs; default `01_intrasample` reads A's unscored calls)
- `remote_outdir` = `gs://<bucket>/…/pilot/D`
- `suffixes` = `["ultralong"]` first, then re-run with `["ultralong","bnd"]`
- `chromosomes` = `["chr1","chr2","chr6","chr19"]`
- `reference_fa` / `reference_fai`
- `n_expected_samples` = **unset** (auto-derived per suffix)
- rest default

**Check:** `gsutil ls gs://<bucket>/…/pilot/D/ultralong/15_concat/truvari_collapsed.bcf`
(and `/bnd/…` once added). The collapsed records carry `FORMAT/SCORE` +
`FORMAT/CALIBRATION_SENSITIVITY` (per sample); apply your own cutoff downstream.

---

## Pilot C — `SV_Integration_WorkflowC_Regenotype` (per family) + Workflow E (gather)

**Regenotyping is a two-submission flow:** C force-calls each family (one Terra
instance per family-set, all writing to one SHARED `remote_outdir`), then
Workflow E merges every family's per-sample chunks into the single cohort VCF.
C by itself does NOT emit a cohort VCF — E does.

**Step C** — run **with inputs defined by data table**, root entity = a
`sample_set` (one set == one family). Select all the family-sets you want; Terra
scatters one instance each, all sharing `remote_outdir`.

| Input | Expression |
|---|---|
| `family_id` | `this.<set_type>_id` (or the set name literally) |
| `sample_ids` | `this.samples.sample_id` |
| `sample_sexes` | `this.samples.sex` |
| `aligned_bais` / `aligned_bams` | `this.samples.02_aligned_bai` / `this.samples.02_aligned_bam` |

Literals:
- `remote_indir` = `gs://<bucket>/…/pilot/B`  (B's remote_outdir; C reads its `/06_concat` subdir automatically)
- `remote_outdir` = `gs://<bucket>/…/pilot/C_regeno`  (SHARED across all family-sets)
- `ped` = the cohort PED. **Mandatory**: each `family_id` (set id) must list
  exactly that set's members in PED cols 1–2, or the run errors (guards against a
  mis-built set). Passed through to WP7.
- `split_for_bcftools_merge_csv` = the pilot CSV
- `reference_fa` / `reference_fai`, `ploidy_bed_female` / `ploidy_bed_male`, `autosomes_bed`
- `requester_pays_project` if needed

**Check (C):** `gsutil ls gs://<bucket>/…/pilot/C_regeno/` shows `<sample>.done`
+ `chunk_*/` for every family member across all sets.

**Step E** — `SV_Integration_WorkflowE_Regenotype_Merge`, **inputs by file paths**,
after all C instances finish:
- `remote_indir` = `gs://<bucket>/…/pilot/C_regeno`  (C's shared `remote_outdir`)
- `remote_outdir` = `gs://<bucket>/…/pilot/C_cohort`
- `split_for_bcftools_merge_csv` = same CSV; `merge_mode` = 2 (default)
- leave `sample_ids_file` unset (auto-derived from the `.done` markers)

**Check (E):** `gsutil ls gs://<bucket>/…/pilot/C_cohort/concat/merged.bcf`;
`bcftools query -l` == all samples across all families.

## Ultralong cuteFC regenotyping test (`SV_Integration_WorkflowD_Cutefc_Test`) + E

Same two-step shape, for the fork-under-test on the ultralong branch. The cohort
ultralong callset must already exist from a **Workflow D (ultralong)** run.

**Force-call** — data table, root = family `sample_set` (one per family):
`family_id`/`sample_ids`/`aligned_bais`/`aligned_bams` bound as in C (no
`sample_sexes` — cuteFC ignores ploidy); `ped` mandatory (same validation);
`remote_indir` = Workflow D's `remote_outdir` (reads `/ultralong/15_concat`);
`remote_outdir` = a shared dir; `cutefc_docker_image` = `quay.io/ymostovoy/lr-ultralong:latest`.

**Gather** — run Workflow E on that shared `remote_outdir` (merge_mode 2) →
`concat/merged.bcf`. Same E definition as the main branch (merge-by-ID is
type-agnostic); ultralong records are bulkier, so bump `merge_disk_size_gb` /
`merge_ram_size_gb` if a chunk is tight.

---

## Scaling up

Once all pilots produce sane BCFs, swap: the pilot CSV → full
`hg38_split_for_bcftools_merge.csv`; the chromosome lists → all 24; the pilot
`sample_set` → the real cohort; and raise `batch_size` in A and A2. Steps are idempotent
via GCS `.done` markers, so resubmitting only reruns failed units.

## Troubleshooting a failed submission

1. `terra_get_submission` → find the failed workflow id + which task.
2. `terra_get_workflow_metadata` → `failures` message + `callsSummary`.
3. The MCP `stderr_tail` is only the **head** of stderr; with `set -x` the real
   error is at the END. Read it directly:
   `gsutil cat <workflowRoot>/call-<Task>/shard-<N>/stderr | tail -60`
   and check `…/shard-<N>/monitoring.log` for OOM/disk.
