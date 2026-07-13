# SV Integration WDL Rename Log

The main short-SV branch WDLs used to be numbered non-consecutively (WP1, 3, 5,
6, 7, 8, 9, 11) because intermediate workpackages had been dropped over time.
They were renamed to a consecutive WP1–WP8 sequence that matches the run order
in `SV_Integration_Run_Order.md`.

The file renames were done first; this log also records the follow-up pass that
made the internal `workflow` names (and stale in-file / doc references) match the
new file names. Task names were all generic (`Impl`, `SingleChromosome`,
`AllChromosomes`) and did not change.

## File name and workflow name mapping

| Run order | Old file | New file | Old `workflow` name | New `workflow` name |
|-----------|----------|----------|---------------------|---------------------|
| 1 | `SV_Integration_Workpackage1_KANPIG_DIRECT_REMOTE_ACCESS.wdl` | `SV_Integration_Workpackage1_intrasample_merge.wdl` | `SV_Integration_Workpackage1` | `SV_Integration_Workpackage1` *(unchanged)* |
| 2 | `SV_Integration_Workpackage3.wdl` | `SV_Integration_Workpackage2_Main_scoring.wdl` | `SV_Integration_Workpackage3` | `SV_Integration_Workpackage2` |
| 3 | `SV_Integration_Workpackage5.wdl` | `SV_Integration_Workpackage3_Main_bcftools_merge.wdl` | `SV_Integration_Workpackage5` | `SV_Integration_Workpackage3` |
| 4 | `SV_Integration_Workpackage6.wdl` | `SV_Integration_Workpackage4_Main_shard.wdl` | `SV_Integration_Workpackage6` | `SV_Integration_Workpackage4` |
| 5 | `SV_Integration_Workpackage7.wdl` | `SV_Integration_Workpackage5_Main_truvari_collapse.wdl` | `SV_Integration_Workpackage7` | `SV_Integration_Workpackage5` |
| 6 | `SV_Integration_Workpackage8.wdl` | `SV_Integration_Workpackage6_Main_concat_shards.wdl` | `SV_Integration_Workpackage8` | `SV_Integration_Workpackage6` |
| 7 | `SV_Integration_Workpackage9_families.wdl` | `SV_Integration_Workpackage7_Main_joint_genotype_families.wdl` | `SV_Integration_Workpackage9_families` | `SV_Integration_Workpackage7_families` |
| 8 | `SV_Integration_Workpackage11.wdl` | `SV_Integration_Workpackage8_Main_concat_regenotyped_shards.wdl` | `SV_Integration_Workpackage11` | `SV_Integration_Workpackage8` |

Optional / test-mode variants of the family regenotyping step (WP7):

| Old file | New file | Old `workflow` name | New `workflow` name |
|----------|----------|---------------------|---------------------|
| `SV_Integration_Workpackage9_trios.wdl` | `SV_Integration_Workpackage7_trios.wdl` | `SV_Integration_Workpackage9_trios` | `SV_Integration_Workpackage7_trios` |
| *(new, this session)* | `SV_Integration_Workpackage7_Main_joint_genotype_families_cutefc.wdl` | — | `SV_Integration_Workpackage7_families_cutefc` |

## Other reference fixes made for consistency

- The WP7 family/trios/cutefc workflows described their input cohort VCF as the
  "WP8" cohort. Under the old numbering the genome-wide `truvari_collapsed.bcf`
  was produced by `SV_Integration_Workpackage8` (concat shards), which is now
  **WP6** (`SV_Integration_Workpackage6_Main_concat_shards.wdl`). Those comments
  and `parameter_meta` entries were updated from `WP8` to `WP6`.
- `SV_Integration_Run_Order.md` step 8 referenced the old
  `SV_Integration_Workpackage11.wdl`; it now points to
  `SV_Integration_Workpackage8_Main_concat_regenotyped_shards.wdl`.
- `SV_Integration_Workpackage7_trios.wdl` called itself a "test-mode sibling of
  Workpackage9_families"; updated to "Workpackage7_families".
