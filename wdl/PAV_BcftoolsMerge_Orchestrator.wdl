version 1.0

import "PAV_BcftoolsMerge_2_Merge.wdl" as merge_wf
import "PAV_BcftoolsMerge_3_Concat.wdl" as concat_wf


# End-to-end PAV whole-callset (SNV+indel+SV) `bcftools merge` pipeline,
# STAGES 2+3 ONLY: runs stage 2 (merge) -> stage 3 (concat) as one submission.
#
# Stage 1 (split) is deliberately NOT part of this orchestrator. Each center
# has its own Terra data table, and a single Terra submission can only iterate
# one data table at a time (`this.samples.sample_id` resolves against
# whichever table is the submission's root entity) -- so stage 1 is run
# separately, once per center, submitted directly against that center's own
# table (see `PAV_BcftoolsMerge_1_Split.wdl`), each pointed at its own
# `remote_outdir`. Once every enabled center's stage-1 run has completed, THIS
# workflow picks up from their output directories and runs stage 2 -- which is
# center-aware; see `PAV_BcftoolsMerge_2_Merge.wdl` for the exact
# overlap-resolution rule (localization order + the bi_samples_to_prefer_over_ha
# override) -- followed by stage 3.
#
# Stage ordering is enforced by a real data dependency: Concat consumes
# Merge's `done` signal.
#
# Each stage is also runnable on its own (see the individual WDLs) for reruns,
# monitoring, or resource retuning at 13k-sample scale.
#
workflow PAV_BcftoolsMerge_Orchestrator {
    input {
        File split_for_bcftools_merge_csv

        String remote_indir_bi
        String remote_indir_ha
        String remote_indir_bcm
        String remote_indir_uw
        String remote_indir_controls_15x
        String remote_indir_controls_30x

        Array[String] bi_samples_to_prefer_over_ha

        Int n_expected_samples_bi
        Int n_expected_samples_ha
        Int n_expected_samples_bcm
        Int n_expected_samples_uw
        Int n_expected_samples_controls_15x
        Int n_expected_samples_controls_30x

        String remote_outdir

        File? sample_ids_file

        Int n_files_per_merge = 100
        Int merge_disk_size_gb = 256
        Int merge_ram_size_gb = 32
        Int merge_n_cpu = 8

        Array[String] contigs = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        Int concat_disk_size_gb = 512

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        split_for_bcftools_merge_csv: "The SAME CSV every enabled center's stage-1 run used. Its line count is the number of chunks; chunk id == 0-based line number."
        remote_indir_bi: "Without final slash. bi's stage-1 output dir (from a separate PAV_BcftoolsMerge_1_Split run against bi's own Terra table). Set n_expected_samples_bi=0 and this to anything if bi is unused."
        remote_indir_ha: "Without final slash. ha's stage-1 output dir. Set n_expected_samples_ha=0 if unused."
        remote_indir_bcm: "Without final slash. bcm's stage-1 output dir. Set n_expected_samples_bcm=0 if unused."
        remote_indir_uw: "Without final slash. uw's stage-1 output dir. Set n_expected_samples_uw=0 if unused."
        remote_indir_controls_15x: "Without final slash. 15x-controls stage-1 output dir. Set n_expected_samples_controls_15x=0 if unused."
        remote_indir_controls_30x: "Without final slash. 30x-controls stage-1 output dir. Set n_expected_samples_controls_30x=0 if unused."
        bi_samples_to_prefer_over_ha: "sample_ids present in BOTH bi and ha for which bi's copy should win instead of the default (ha wins on overlap). No effect on samples not present in both."
        n_expected_samples_bi: "Exact expected file count in bi's chunk dirs. 0 disables bi entirely."
        n_expected_samples_ha: "Exact expected file count in ha's chunk dirs. 0 disables ha entirely."
        n_expected_samples_bcm: "Exact expected file count in bcm's chunk dirs. 0 disables bcm entirely."
        n_expected_samples_uw: "Exact expected file count in uw's chunk dirs. 0 disables uw entirely."
        n_expected_samples_controls_15x: "Minimum expected file count in the 15x-controls chunk dirs (more than expected only warns). 0 disables this source entirely."
        n_expected_samples_controls_30x: "Minimum expected file count in the 30x-controls chunk dirs (more than expected only warns). 0 disables this source entirely."
        remote_outdir: "Without final slash. Stages write to `<remote_outdir>/02_merge`, `/03_concat`."
        sample_ids_file: "Optional. Merge column order (the UNION of sample_ids across all enabled sources). If omitted, stage 2 derives it by unioning `sample_ids.txt` from every enabled source."
    }

    String merge_outdir = remote_outdir + "/02_merge"
    String concat_outdir = remote_outdir + "/03_concat"

    call merge_wf.PAV_BcftoolsMerge_2_Merge as Merge {
        input:
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir_bi = remote_indir_bi,
            remote_indir_ha = remote_indir_ha,
            remote_indir_bcm = remote_indir_bcm,
            remote_indir_uw = remote_indir_uw,
            remote_indir_controls_15x = remote_indir_controls_15x,
            remote_indir_controls_30x = remote_indir_controls_30x,
            bi_samples_to_prefer_over_ha = bi_samples_to_prefer_over_ha,
            n_expected_samples_bi = n_expected_samples_bi,
            n_expected_samples_ha = n_expected_samples_ha,
            n_expected_samples_bcm = n_expected_samples_bcm,
            n_expected_samples_uw = n_expected_samples_uw,
            n_expected_samples_controls_15x = n_expected_samples_controls_15x,
            n_expected_samples_controls_30x = n_expected_samples_controls_30x,
            remote_outdir = merge_outdir,
            sample_ids_file = sample_ids_file,
            n_files_per_merge = n_files_per_merge,
            merge_disk_size_gb = merge_disk_size_gb,
            merge_ram_size_gb = merge_ram_size_gb,
            merge_n_cpu = merge_n_cpu,
            docker_image = docker_image
    }

    call concat_wf.PAV_BcftoolsMerge_3_Concat as Concat {
        input:
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = merge_outdir,
            remote_outdir = concat_outdir,
            contigs = contigs,
            concat_disk_size_gb = concat_disk_size_gb,
            upstream_signal = Merge.done,
            docker_image = docker_image
    }

    output {
        File merged_sample_ids_file = Merge.merged_sample_ids_file
        Array[String] per_contig_bcfs = Concat.per_contig_bcfs
    }
}
