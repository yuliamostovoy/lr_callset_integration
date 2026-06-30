# SV Integration Run Order

This is a high-level runbook for the three SV integration branches. The short-SV
main branch produces the regular integrated cohort. The BND and ultralong
branches add feature annotation and XGBoost scoring for those special record
classes. Scores are retained; they are not used to remove variants.

## Main Short-SV Branch

Run these workflows in order:

1. `SV_Integration_Workpackage1_intrasample_merge.wdl`
2. `SV_Integration_Workpackage2_Main_scoring.wdl`
3. `SV_Integration_Workpackage3_Main_bcftools_merge.wdl`
4. `SV_Integration_Workpackage4_Main_shard.wdl`
5. `SV_Integration_Workpackage5_Main_truvari_collapse.wdl`
6. `SV_Integration_Workpackage6_Main_concat_shards.wdl`
7. `SV_Integration_Workpackage7_Main_joint_genotype_families.wdl`
8. `SV_Integration_Workpackage11.wdl`

What they do:

- `SV_Integration_Workpackage1_intrasample_merge.wdl`: per-sample integration. Canonicalizes PAV,
  PBSV, and Sniffles calls; removes unsupported regions and short records;
  collapses/merges intra-sample calls; runs Kanpig; writes per-sample `kanpig`,
  `training`, `bnd`, and `ultralong` outputs.
- `SV_Integration_Workpackage2_Main_scoring.wdl`: per-sample XGBoost scoring for the regular
  Kanpig SV VCF. It copies `SCORE` and `CALIBRATION_SENSITIVITY` into FORMAT,
  keeps all records, and scatters each sample into merge chunks.
- `SV_Integration_Workpackage3_Main_bcftools_merge.wdl`: cohort-level `bcftools merge` for each
  regular SV chunk.
- `SV_Integration_Workpackage4_Main_shard.wdl`: divides merged regular SV chromosomes into
  Truvari-collapse chunks.
- `SV_Integration_Workpackage5_Main_truvari_collapse.wdl`: runs Truvari collapse on each regular 
  chunk.
- `SV_Integration_Workpackage6_Main_concat_shards.wdl`: concatenates collapsed regular SV chunks
  across chromosomes, assigns final IDs, adds `N_DISCOVERY_SAMPLES`, and writes
  the genome-wide `truvari_collapsed.bcf` consumed by the family regenotyping
  workflow.
- `SV_Integration_Workpackage7_Main_joint_genotype_families.wdl`: builds family-specific candidate
  VCFs from the WP8 collapsed cohort VCF, keeps sites present in at least one
  family member, regenotypes each family member with Kanpig, and scatters each
  sample's regenotyped calls into merge chunks. This is the production rare
  disease path and replaces the older personalized VCF workflow.
- `SV_Integration_Workpackage11.wdl`: concatenates regenotyped merge chunks into
  final merged regenotyped BCFs.

Optional/test variant:

- `SV_Integration_Workpackage7_trios.wdl`: builds trio-specific candidate VCFs
  from complete mother/father/proband trios and runs `kanpig trio`. By default
  it uploads proband chunks only, with optional parent upload. This is a
  test-mode sibling of the family workflow, not the default production path.

## BND Scoring Branch

Run these workflows after `SV_Integration_Workpackage1.wdl` has produced
per-sample BND outputs:

1. `SV_Integration_BndAnnotate.wdl`
2. `SV_Integration_BndGetTrainingIntervals.wdl`
3. `UltralongMerge.wdl` for the annotated BND VCFs
4. `UltralongMerge.wdl` for the BND training VCFs
5. `UltralongScore.wdl` with `svtype = "bnd"`

What they do:

- `SV_Integration_BndAnnotate.wdl`: converts each per-sample BND BCF to VCF,
  cleans BND REF/ALT/QUAL/FILTER fields, adds BAM-derived breakpoint features,
  optional cuteFC features, and feature-extraction annotations, then writes
  `${sample_id}_bnd.vcf.gz`.
- `SV_Integration_BndGetTrainingIntervals.wdl`: matches annotated BND calls to
  assembly-derived svim-asm BND truth using Truvari and writes
  `${sample_id}_bnd_training.vcf.gz`.
- `UltralongMerge.wdl` on annotated BNDs: concatenates/reheaders all
  `${sample_id}_bnd.vcf.gz` files into a single BND query VCF for scoring.
- `UltralongMerge.wdl` on BND training records: concatenates/reheaders all
  `${sample_id}_bnd_training.vcf.gz` files into the BND training resource VCF.
- `UltralongScore.wdl`: trains and applies the BND scoring model. It writes
  scored BND VCFs with score fields retained in FORMAT and does not filter by
  score.

If a non-scored BND cohort merge is needed, `SV_Integration_Workpackage5_Bnd.wdl`
still performs the older simple all-sample BND merge from per-sample BND VCFs.

## Ultralong Scoring Branch

Run these workflows after the ultralong inputs and BAMs are available:

1. `SV_Integration_UltralongAnnotate.wdl`
2. `SV_Integration_UltralongGetTrainingIntervals.wdl`
3. `UltralongMerge.wdl` for each annotated ultralong SVTYPE
4. `UltralongMerge.wdl` for each corresponding training SVTYPE
5. `UltralongScore.wdl` for each SVTYPE

What they do:

- `SV_Integration_UltralongAnnotate.wdl`: canonicalizes per-sample ultralong
  records, optionally converts INS-derived DUP evidence, adds BAM-derived and
  genotyper-derived annotations, and splits outputs by SVTYPE. Expected output
  suffixes are `del`, `dup`, `inv`, `ins`, and `insdup`.
- `SV_Integration_UltralongGetTrainingIntervals.wdl`: compares annotated
  ultralong records to assembly-derived truth and writes per-sample
  `${sample_id}_${svtype}_training.vcf.gz` files.
- `UltralongMerge.wdl` on annotated calls: for each `svtype`, merges
  `${sample_id}_${svtype}.vcf.gz` into one query VCF for scoring.
- `UltralongMerge.wdl` on training calls: for each `svtype`, merges
  `${sample_id}_${svtype}_training.vcf.gz` into the matching training resource
  VCF.
- `UltralongScore.wdl`: trains and applies the scoring model for the requested
  `svtype`. Run separately for each ultralong class to score. It retains all
  records and writes `SCORE` and `CALIBRATION_SENSITIVITY` into FORMAT.

For the older ultralong cohort-discovery path, use:

1. `SV_Integration_Workpackage12.wdl`
2. `SV_Integration_Workpackage13.wdl`
3. `SV_Integration_Workpackage14.wdl`
4. `SV_Integration_Workpackage15.wdl`

Those workflows merge per-sample ultralong or BND BCFs, split chromosomes into
collapse chunks, run Truvari collapse, and concatenate final collapsed BCFs.
They are separate from the annotation/scoring workflows above.
