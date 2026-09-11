version 1.0

import "SV_Integration_BndGetTrainingIntervalsPrime.wdl" as base


# DATA-TABLE FRONT-END for `SV_Integration_BndGetTrainingIntervalsPrime`.
#
# Like `SV_Integration_BndGetTrainingIntervals_FromTable`, but for the "Prime"
# labeler, which reads breakpoints straight from the assembly alignments instead
# of svim-asm's BND calls. The upstream workflow's only used `samples_tsv` column
# is the sample ID; this wrapper takes `sample_ids` as a Terra data-table column,
# builds the one-column TSV, batches, and scatters the unmodified workflow.
#
# Per sample it pulls, by naming convention:
#   remote_indir_query/<sample_id>_bnd.vcf.gz        (BndAnnotate output)
#   remote_indir_truth/<sample_id>_breakpoints.csv   (BndBuildTruth output)
#
workflow SV_Integration_BndGetTrainingIntervalsPrime_FromTable {
    input {
        # --- Per-sample column from the Terra data table ---
        Array[String] sample_ids

        Int batch_size = 20

        # --- Underlying-workflow pass-through ---
        Int breakpoint_max_distance = 500
        File reference_agp

        String remote_indir_query
        String remote_indir_truth
        String remote_outdir

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`."
        batch_size: "Number of samples processed sequentially per VM."
        remote_indir_query: "BndAnnotate output dir holding `<sample_id>_bnd.vcf.gz`. Trailing slashes stripped."
        remote_indir_truth: "BndBuildTruth output dir holding `<sample_id>_breakpoints.csv`. Trailing slashes stripped."
        remote_outdir: "Where `<sample_id>_bnd_training*.vcf.gz` is written. Trailing slashes stripped. Feeds Merge (svtype=bnd, suffix='_training')."
    }

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            batch_size = batch_size,
            docker_image = docker_image
    }

    scatter (manifest in MakeManifests.manifests) {
        call base.SV_Integration_BndGetTrainingIntervalsPrime as GetIntervals {
            input:
                samples_tsv = manifest,
                breakpoint_max_distance = breakpoint_max_distance,
                reference_agp = reference_agp,
                remote_indir_query = sub(remote_indir_query, "/+$", ""),
                remote_indir_truth = sub(remote_indir_truth, "/+$", ""),
                remote_outdir = sub(remote_outdir, "/+$", ""),
                docker_image = docker_image
        }
    }

    output {
        String bnd_training_intervals_outdir = sub(remote_outdir, "/+$", "")
    }
}


# Builds the one-column (`ID`/line) samples TSV from `sample_ids`, then splits
# into `batch_size`-row manifests (one per VM).
task MakeManifests {
    input {
        Array[String] sample_ids
        Int batch_size
        String docker_image
    }

    command <<<
        set -euxo pipefail

        cp ~{write_lines(sample_ids)} all.tsv
        split --lines=~{batch_size} --numeric-suffixes=0 --suffix-length=6 all.tsv batch_
        ls batch_* 1>&2
    >>>

    output {
        Array[File] manifests = glob("batch_*")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}
