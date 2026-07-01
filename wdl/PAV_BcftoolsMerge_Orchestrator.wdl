version 1.0

import "PAV_BcftoolsMerge_1_Split.wdl" as split_wf
import "PAV_BcftoolsMerge_2_Merge.wdl" as merge_wf
import "PAV_BcftoolsMerge_3_Concat.wdl" as concat_wf


# End-to-end PAV whole-callset (SNV+indel+SV) `bcftools merge` pipeline: runs
# stage 1 (split) -> stage 2 (merge) -> stage 3 (concat) as one submission.
#
# Input is just the Terra data table (sample_id + PAV VCF URI) plus the split
# CSV; there are no hand-built TSVs between stages. Stage ordering is enforced by
# real data dependencies:
#   - Merge consumes Split's `sample_ids_file` (so it waits for Split and reuses
#     the sample list without an extra download).
#   - Concat consumes Merge's `done` signal.
#
# Each stage is also runnable on its own (see the individual WDLs) for reruns,
# monitoring, or resource retuning at 13k-sample scale.
#
workflow PAV_BcftoolsMerge_Orchestrator {
    input {
        Array[String] sample_ids
        Array[String] pav_vcfs

        File reference_fa
        File reference_fai
        File split_for_bcftools_merge_csv

        String norm_remote_dir
        String remote_outdir

        Int batch_size = 100
        Boolean left_align_and_split = true
        String check_ref = "w"

        Int n_files_per_merge = 100
        Int merge_disk_size_gb = 256
        Int merge_ram_size_gb = 32
        Int merge_n_cpu = 8

        Array[String] contigs = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        Int concat_disk_size_gb = 512

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.samples.sample_id`. Parallel to `pav_vcfs`."
        pav_vcfs: "From the data table, e.g. `this.samples.pav_vcf`. Array[String] gs:// URIs (streamed in-task)."
        split_for_bcftools_merge_csv: "0-based half-open chunk partition; chunk id == 0-based line number. Ship the ~5 Mbp CSV for the SNV-inclusive callset."
        norm_remote_dir: "Without final slash. STABLE location for normalized per-sample BCFs, reused across runs. Keep it constant (independent of `remote_outdir`) so re-chunking / pilot->full never re-normalizes."
        remote_outdir: "Without final slash. Stages write to `<remote_outdir>/01_split`, `/02_merge`, `/03_concat`."
    }

    String split_outdir = remote_outdir + "/01_split"
    String merge_outdir = remote_outdir + "/02_merge"
    String concat_outdir = remote_outdir + "/03_concat"

    call split_wf.PAV_BcftoolsMerge_1_Split as Split {
        input:
            sample_ids = sample_ids,
            pav_vcfs = pav_vcfs,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            batch_size = batch_size,
            norm_remote_dir = norm_remote_dir,
            remote_outdir = split_outdir,
            left_align_and_split = left_align_and_split,
            check_ref = check_ref,
            docker_image = docker_image
    }

    call merge_wf.PAV_BcftoolsMerge_2_Merge as Merge {
        input:
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = split_outdir,
            remote_outdir = merge_outdir,
            sample_ids_file = Split.sample_ids_file,
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
        File sample_ids_file = Split.sample_ids_file
        Array[String] per_contig_bcfs = Concat.per_contig_bcfs
    }
}
