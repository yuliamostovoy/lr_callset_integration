# SV Integration Workpackage Summary

This summarizes the major processing steps run by the `SV_Integration_Workpackage*.wdl` workflows in this directory.

## Workpackage 1: single-sample integration

Source: `SV_Integration_Workpackage1.wdl`

1. Build reference non-gap BED from the AGP file.
2. For each sample, localize caller VCFs and BAM/index as needed.
3. Canonicalize PAV, PBSV, and Sniffles calls: normalize symbolic records, clean REF/ALT/QUAL, filter chromosomes/gaps/lengths.
4. Intra-sample merge standard SVs with `bcftools merge` and `truvari collapse`.
5. Separately merge ultralong calls and BND calls.
6. Copy caller support fields into INFO.
7. Run `kanpig gt` on the intra-sample merged SV VCF.
8. Copy Kanpig fields into INFO, and generate Kanpig BED/CSV QC.
9. Extract training records against the training resource.
10. Copy support fields into ultralong/BND outputs and upload `kanpig`, `training`, `ultralong`, and `bnd` outputs.

Variants:

- `SV_Integration_Workpackage1_COMPRESSED.wdl` is the older/compressed-output version of the same flow.
- `SV_Integration_Workpackage1_KANPIG_DIRECT_REMOTE_ACCESS.wdl` is the same flow but runs Kanpig against a remote BAM path instead of localizing the BAM first.

## Workpackage 3: per-sample XGBoost scoring and chunking

Source: `SV_Integration_Workpackage3.wdl`

1. For each sample, localize the Workpackage 1 Kanpig VCF and training VCF.
2. Extract annotations with GATK `ExtractVariantAnnotations`.
3. Train a per-sample annotation model with `TrainVariantAnnotationsModel`.
4. Score all sample variants with `ScoreVariantAnnotations`.
5. Copy `SUPP_*`, `SCORE`, and `CALIBRATION_SENSITIVITY` from INFO to FORMAT.
6. Emit debug counts for score thresholds.
7. Split the scored BCF into `bcftools merge` chunks and upload each chunk.

## Workpackage 5: cohort merge of regular SV chunks

Source: `SV_Integration_Workpackage5.wdl`

1. For one chunk, verify expected sample counts across BI/HA/BCM/UW/control inputs.
2. Check disk space against remote file sizes.
3. Localize all per-sample BCFs for that chunk, resolving BI-over-HA preferences.
4. Do a two-level `bcftools merge`.
5. If standard merge mode, split multiallelics with `bcftools norm`.
6. If ID-merge mode, remove records that are REF in all samples and write ALT-count QC.
7. Upload merged chunk BCF/index.

## Workpackage 5_Bnd: older/simple BND merge

Source: `SV_Integration_Workpackage5_Bnd.wdl`

1. Download all per-sample BND VCFs.
2. Merge with `bcftools merge`.
3. Split multiallelics.
4. Remove exact duplicates.
5. Upload final BND VCF.

## Workpackage 6: make Truvari-collapse chunks for regular SVs

Source: `SV_Integration_Workpackage6.wdl`

1. Localize all `bcftools merge` chunks for one chromosome.
2. Concatenate them into one chromosome BCF.
3. Query POS/REF/ALT.
4. Use `TruvariDivide2` to compute collapse regions.
5. Slice chromosome BCF into collapse chunks.
6. Run consistency checks, then upload chunks and `regions.txt`.

Variant:

- `SV_Integration_Workpackage6_TRUVARI_DIVIDE.wdl` does the older version using `truvari divide` directly on a whole-chromosome VCF.

## Workpackage 7: run Truvari collapse on regular SV chunks

Source: `SV_Integration_Workpackage7.wdl`

1. Optionally localize an included BED.
2. For each collapse chunk, download BCF/index.
3. Set QUAL to number of samples with ALT genotype.
4. Run `truvari collapse` with `--keep maxqual --gt off`.
5. Sort/index collapsed output.
6. Remove long Truvari-generated IDs.
7. Upload collapsed chunk and done marker.

## Workpackage 8: concatenate collapsed regular SVs

Source: `SV_Integration_Workpackage8.wdl`

1. Scatter by chromosome.
2. Per chromosome: download collapsed chunks and concatenate.
3. Assign globally unique chromosome-prefixed IDs.
4. Add `N_DISCOVERY_SAMPLES`.
5. Save per-chromosome `truvari_collapsed.bcf`.
6. Across chromosomes: concatenate per-chromosome BCFs into genome-wide `truvari_collapsed.bcf`.

Variant:

- `SV_Integration_Workpackage8_Prime.wdl` performs the older frequent/infrequent split from an existing cohort Truvari VCF instead of per-chromosome chunk outputs.

## Workpackage 9: personalized cohort VCF regenotyping

Source: `SV_Integration_Workpackage9.wdl`

1. Download genome-wide `frequent` and `infrequent` BCFs.
2. Split `infrequent.bcf` by sample.
3. For each sample, localize BAM/index.
4. Keep that sample's ALT infrequent records and concatenate with frequent records to make a personalized VCF.
5. Run `kanpig gt` on the personalized VCF.
6. Sort/index and emit QC.
7. Split regenotyped VCF into merge chunks and upload.

## Workpackage 11: concatenate regenotyped merge chunks

Source: `SV_Integration_Workpackage11.wdl`

1. Localize specified chunk BCFs.
2. Concatenate them in order with `bcftools concat --naive`.
3. Index and upload `merged.bcf`.

## Workpackage 12: cohort merge ultralong or BND per-sample BCFs

Source: `SV_Integration_Workpackage12.wdl`

1. Verify expected sample counts across input datasets.
2. Check disk space and localize all per-sample `*_ultralong.bcf` or `*_bnd.bcf`.
3. First-level whole-genome merge in batches of samples.
4. Normalize/split multiallelics.
5. Split each batch merge by chromosome.
6. Second-level per-chromosome merge.
7. Normalize per chromosome.
8. Upload one merged BCF/index per chromosome.

## Workpackage 13: split ultralong/BND chromosome BCFs for collapse

Source: `SV_Integration_Workpackage13.wdl`

1. Download one chromosome BCF.
2. For ultralong, optionally remove symbolic `<INS>` records.
3. Query POS/REF/ALT.
4. Use `TruvariDivide2Ultralong` to compute regions.
5. Slice chromosome BCF into collapse chunks.
6. Optionally verify record counts and IDs.
7. Upload chunk BCFs.

## Workpackage 14: run Truvari collapse on ultralong chunks

Source: `SV_Integration_Workpackage14.wdl`

1. Optionally localize an included BED.
2. For each chunk, download BCF/index.
3. Reset symbolic ALT strings to avoid SVLEN interference.
4. Set QUAL to number of ALT samples.
5. Run `truvari collapse`, defaulting to sequence similarity off.
6. Sort/index.
7. Rewrite IDs to compact chunk-local IDs.
8. Upload collapsed chunk and done marker.

## Workpackage 15: concatenate collapsed ultralong/BND chunks

Source: `SV_Integration_Workpackage15.wdl`

1. Scatter by chromosome.
2. Per chromosome: download collapsed chunks.
3. Concatenate chunks, or rename the single chunk if only one exists.
4. Assign unique chromosome-prefixed IDs.
5. Add `N_DISCOVERY_SAMPLES`.
6. Upload per-chromosome `truvari_collapsed.bcf`.
7. Across chromosomes: concatenate per-chromosome BCFs into final genome-wide `truvari_collapsed.bcf`.

## Missing numbered SV integration workpackages

I did not find `SV_Integration_Workpackage2.wdl`, `SV_Integration_Workpackage4.wdl`, or `SV_Integration_Workpackage10.wdl` in this directory.
