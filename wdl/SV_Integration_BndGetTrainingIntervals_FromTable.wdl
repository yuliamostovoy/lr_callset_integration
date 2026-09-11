version 1.0

import "SV_Integration_BndGetTrainingIntervals.wdl" as base


# DATA-TABLE FRONT-END for `SV_Integration_BndGetTrainingIntervals`.
#
# The upstream workflow ingests samples through a hand-built `samples_tsv` whose
# only used column is the sample ID. This wrapper takes `sample_ids` as a Terra
# data-table column, builds the one-column samples TSV, batches into
# `batch_size`-sample VMs, and scatters the unmodified workflow over the batches.
#
# Per sample it pulls, by naming convention:
#   remote_indir_svimasm/<sample_id>_canonized.vcf.gz  (SvimAsm output)
#   remote_indir_query/<sample_id>_bnd.vcf.gz           (BndAnnotate output)
#
workflow SV_Integration_BndGetTrainingIntervals_FromTable {
    input {
        # --- Per-sample column from the Terra data table ---
        Array[String] sample_ids

        Int batch_size = 20

        # --- Underlying-workflow pass-through ---
        Int remove_orientations = 1
        Int use_interval_svtypes = 0
        Int interval_svtypes_min_sv_length = 1000
        Int truvari_bnddist = 500

        String remote_indir_query
        String remote_indir_svimasm
        String remote_outdir

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`."
        batch_size: "Number of samples processed sequentially per VM."
        remote_indir_query: "BndAnnotate output dir holding `<sample_id>_bnd.vcf.gz`. Trailing slashes stripped."
        remote_indir_svimasm: "SvimAsm output dir holding `<sample_id>_canonized.vcf.gz`. Trailing slashes stripped."
        remote_outdir: "Where `<sample_id>_bnd_training.vcf.gz` is written. Trailing slashes stripped. Feeds Merge (svtype=bnd, suffix='_training')."
    }

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            batch_size = batch_size,
            docker_image = docker_image
    }

    scatter (manifest in MakeManifests.manifests) {
        call base.SV_Integration_BndGetTrainingIntervals as GetIntervals {
            input:
                samples_tsv = manifest,
                remove_orientations = remove_orientations,
                use_interval_svtypes = use_interval_svtypes,
                interval_svtypes_min_sv_length = interval_svtypes_min_sv_length,
                truvari_bnddist = truvari_bnddist,
                remote_indir_query = sub(remote_indir_query, "/+$", ""),
                remote_indir_svimasm = sub(remote_indir_svimasm, "/+$", ""),
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
