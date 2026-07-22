version 1.0

import "SV_Integration_WorkflowD_PerSuffix.wdl" as perSuffix


# CONSOLIDATED STEP D of the SV integration pipeline: non-scored cohort
# integration of the per-sample ultralong and BND calls (merge -> shard ->
# truvari collapse -> concat; WP12->WP13->WP14->WP15), for BOTH suffixes, in a
# SINGLE Terra submission. NO annotation, NO XGBoost scoring.
#
# Branches off Workflow A's per-sample outputs and runs in parallel with
# Workflow B. Each suffix is handled by the SV_Integration_WorkflowD_PerSuffix
# sub-workflow so no single file exceeds two scatter levels.
#
workflow SV_Integration_WorkflowD_Ultralong_Bnd {
    input {
        Array[String] suffixes = ["ultralong","bnd"]

        String remote_indir
        String remote_outdir

        String chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        Array[String] chromosomes_array = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        Int? n_expected_samples

        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        Int consistency_checks = 1

        File reference_fa
        File reference_fai
        String truvari_matching_parameters = "--refdist 500 --pctseq 0.95 --pctsize 0.95 --pctovl 0.0"
        Int max_resolve = 100000
        Boolean use_bed = false
        Int chunk_ids_per_file = 100

        Int concat_all_naive = 1

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        suffixes: "Which per-sample call classes to integrate. Default: both ultralong and bnd."
        remote_indir: "Workflow A intra-sample output dir holding <sample>_ultralong.bcf and <sample>_bnd.bcf per sample."
        remote_outdir: "Per suffix, stage outputs go to /<suffix>/{12_merge,13_shard,14_collapse,15_concat}; the genome-wide callset is /<suffix>/15_concat/truvari_collapsed.bcf."
        chromosomes: "Comma-separated (WP12 signature)."
        chromosomes_array: "Same chromosomes as an array (WP15 signature)."
        n_expected_samples: "OPTIONAL. Auto-derived per suffix when omitted."
    }

    String indir = sub(remote_indir, "/+$", "")
    String outdir = sub(remote_outdir, "/+$", "")

    scatter (suffix in suffixes) {
        call perSuffix.SV_Integration_WorkflowD_PerSuffix as PerSuffix {
            input:
                suffix = suffix,
                remote_indir = indir,
                remote_outdir_suffix = outdir + "/" + suffix,
                chromosomes = chromosomes,
                chromosomes_array = chromosomes_array,
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
                docker_image = docker_image
        }
    }

    output {
    }
}
