version 1.0

import "SV_Integration_WorkflowD_PerSuffix.wdl" as perSuffix
import "SV_Integration_Workpackage7_Main_joint_genotype_families_cutefc.wdl" as wp7cutefc
import "SV_Integration_Workpackage3_Main_bcftools_merge.wdl" as wp3
import "SV_Integration_Workpackage8_Main_concat_regenotyped_shards.wdl" as wp8


# TEST-ONLY workflow (not production). One submission that:
#   1. Builds the cohort ULTRALONG callset from the per-sample _ultralong.bcf
#      files (WP12->13->14->15, via the PerSuffix sub-workflow), then
#   2. Family-regenotypes it with cuteFC (WP7 cuteFC variant) -- pointed at the
#      forked cuteFC docker under test, and then
#   3. Merges the regenotyped TEST samples back into one cohort VCF
#      (bcftools merge by ID = WP3 merge_mode=2, then WP8 concat).
#
# Only the samples/families provided here are force-called and merged; the rest
# of the cohort is not involved. Stages are ordered by done/upstream_signal
# handshakes. Root the submission on your sample set so the BAMs bind from the
# data table (like Workflow C).
#
workflow SV_Integration_WorkflowD_Cutefc_Test {
    input {
        # --- Cohort ultralong build (WP12-15) ---
        String remote_indir
        String intrasample_subdir = "01_intrasample"
        File split_for_bcftools_merge_csv
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        Int? n_expected_samples
        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        Int consistency_checks = 1
        Int max_resolve = 100000
        Boolean use_bed = false
        Int chunk_ids_per_file = 100
        Int concat_all_naive = 1
        String truvari_matching_parameters = "--refdist 500 --pctseq 0.95 --pctsize 0.95 --pctovl 0.0"

        # --- cuteFC family regenotyping (WP7 cuteFC variant) ---
        Array[String] family_ids
        File ped
        Array[String] sample_ids
        Array[String] aligned_bais
        Array[String] aligned_bams
        File reference_fa
        File reference_fai
        File autosomes_bed
        String requester_pays_project = ""
        String cutefc_params_cohort = "--max_size -1 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"

        # --- output + dockers ---
        String remote_outdir
        String wp_docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
        String cutefc_docker_image = "quay.io/ymostovoy/lr-ultralong:latest"
    }
    parameter_meta {
        remote_indir: "Dir holding the per-sample <sample>_ultralong.bcf files. Workflow A's remote_outdir (reads /01_intrasample) or a legacy WP1 dir with intrasample_subdir=\"\"."
        intrasample_subdir: "Subdir of remote_indir holding <sample>_ultralong.bcf. Default '01_intrasample'; set '' to read remote_indir directly."
        family_ids: "Families to regenotype (each matches column 1 of ped). Bind from your table / provide directly."
        ped: "6-column PED (family_id, sample_id, ...). Defines the test families."
        sample_ids: "Sample IDs, parallel to aligned_bais/aligned_bams. e.g. this.samples.sample_id."
        aligned_bais: "e.g. this.samples.aligned_bai."
        aligned_bams: "e.g. this.samples.aligned_bam."
        remote_outdir: "Test root. Cohort build -> /ultralong/{12_merge,...,15_concat}; force-call -> /regenotype; merge -> /regenotype_merge; final regenotyped cohort VCF -> /regenotype_concat/merged.bcf."
        cutefc_docker_image: "The forked-cuteFC ultralong image to test."
    }

    String indir_root = sub(remote_indir, "/+$", "")
    String ul_indir = if intrasample_subdir == "" then indir_root else indir_root + "/" + intrasample_subdir
    String outdir = sub(remote_outdir, "/+$", "")
    String regeno_dir = outdir + "/regenotype"
    String merge_dir = outdir + "/regenotype_merge"
    String concat_dir = outdir + "/regenotype_concat"

    Int n_chunks = length(read_lines(split_for_bcftools_merge_csv))

    # 1. Cohort ULTRALONG callset (WP12->15).
    call perSuffix.SV_Integration_WorkflowD_PerSuffix as BuildCohort {
        input:
            suffix = "ultralong",
            remote_indir = ul_indir,
            remote_outdir_suffix = outdir + "/ultralong",
            chromosomes = chromosomes,
            n_expected_samples = n_expected_samples,
            truvari_chunk_min_records = truvari_chunk_min_records,
            truvari_collapse_refdist = truvari_collapse_refdist,
            consistency_checks = consistency_checks,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            truvari_matching_parameters = truvari_matching_parameters,
            max_resolve = max_resolve,
            use_bed = use_bed,
            chunk_ids_per_file = chunk_ids_per_file,
            concat_all_naive = concat_all_naive,
            docker_image = wp_docker_image
    }

    # 2. cuteFC family regenotyping of the cohort ultralong sites (fork docker).
    call wp7cutefc.Impl as Regenotype {
        input:
            family_ids = family_ids,
            ped = ped,
            sample_ids = sample_ids,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = BuildCohort.cohort_dir,
            remote_outdir = regeno_dir,
            requester_pays_project = requester_pays_project,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            autosomes_bed = autosomes_bed,
            cutefc_params_cohort = cutefc_params_cohort,
            docker_image = cutefc_docker_image,
            upstream_signal = [BuildCohort.done]
    }

    # 3a. Merge order = the test samples that were force-called.
    call WriteSampleList {
        input:
            regeno_dir = regeno_dir,
            upstream_signal = Regenotype.done,
            docker_image = wp_docker_image
    }
    call MakeChunkIdsCsv {
        input:
            n_chunks = n_chunks,
            docker_image = wp_docker_image
    }

    # 3b. bcftools merge by ID, one VM per chunk.
    scatter (chunk_id in range(n_chunks)) {
        call wp3.Impl as MergeChunk {
            input:
                chunk_id = chunk_id,
                sample_ids = WriteSampleList.sample_ids_file,
                remote_indir = regeno_dir,
                merge_mode = 2,
                remote_outdir = merge_dir,
                docker_image = wp_docker_image,
                upstream_signal = [Regenotype.done]
        }
    }

    # 3c. Concat merged chunks -> one regenotyped cohort VCF (test samples only).
    call wp8.Impl as Concat {
        input:
            chunk_ids = MakeChunkIdsCsv.csv,
            remote_indir = merge_dir,
            remote_outdir = concat_dir,
            docker_image = wp_docker_image,
            upstream_signal = MergeChunk.done
    }

    output {
    }
}


# Merge order / cohort sample list = the samples cuteFC force-called, taken from
# the <sample>.done markers it wrote. Gated on the regenotype stage via
# upstream_signal.
#
task WriteSampleList {
    input {
        String regeno_dir
        String upstream_signal
        String docker_image
    }

    command <<<
        set -euxo pipefail
        echo "gated on: ~{upstream_signal}" 1>&2
        gcloud storage ls ~{regeno_dir}/'*.done' | sed 's#.*/##; s#\.done$##' | sort -u > sample_ids.txt
        wc -l sample_ids.txt 1>&2
        if [ ! -s sample_ids.txt ]; then
            echo "ERROR: no <sample>.done markers found under ~{regeno_dir}."
            exit 1
        fi
    >>>

    output {
        File sample_ids_file = "sample_ids.txt"
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}


# Comma-separated 0..n_chunks-1 for the final WP8 concat.
#
task MakeChunkIdsCsv {
    input {
        Int n_chunks
        String docker_image
    }

    command <<<
        set -euxo pipefail
        seq 0 $(( ~{n_chunks} - 1 )) | paste -sd, - > csv.txt
    >>>

    output {
        String csv = read_string("csv.txt")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}
