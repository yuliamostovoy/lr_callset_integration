version 1.0

import "SV_Integration_UltralongGetTrainingIntervals.wdl" as base


# DATA-TABLE FRONT-END for `SV_Integration_UltralongGetTrainingIntervals`.
#
# The upstream workflow ingests samples through a hand-built `samples_tsv` of
# `ID, dipcall_bed` (column 2 is used only when `match_to_gaps=1`). This wrapper
# takes those as parallel Terra data-table columns, reassembles the exact samples
# TSV, batches into `batch_size`-sample VMs, and scatters the unmodified
# `SV_Integration_UltralongGetTrainingIntervals` workflow over the batches.
#
# Per sample it pulls, by naming convention:
#   remote_indir_svimasm/<sample_id>_svimasm_{del,dup,ins,ins_dup,inv}.vcf.gz  (BuildTruth output)
#   remote_indir_query/<sample_id>_{del,dup,insdup,ins,inv}.vcf.gz              (Annotate output)
# so point those two at the BuildTruth and Annotate output dirs respectively.
# `dipcall_beds` may be left empty ([]) unless match_to_gaps=1.
#
workflow SV_Integration_UltralongGetTrainingIntervals_FromTable {
    input {
        # --- Per-sample columns from the Terra data table ---
        Array[String] sample_ids
        Array[String] dipcall_beds = []

        Int batch_size = 20

        # --- Underlying-workflow pass-through ---
        String remote_indir_query
        String remote_indir_svimasm
        String remote_outdir

        Int convert_ins_to_dup = 1
        Int match_to_gaps = 0
        Int match_insdups_to_dups = 1
        Int match_dups_to_insdups = 0
        Int match_ins_to_dup = 1
        Int match_ins_to_insdup = 0
        Int match_ins_to_dup_slack_bp = 200

        File reference_fai

        Int truvari_refdist = 500
        Float truvari_pctsize_loose = 0.4
        Float truvari_pctsize_strict = 0.9
        Float truvari_pctovl_loose = 0.4

        Int max_read_length = 25000

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`."
        dipcall_beds: "Optional. gs:// dipcall confident-region BEDs, one per sample. Leave empty ([]) unless match_to_gaps=1."
        batch_size: "Number of samples processed sequentially per VM."
        remote_indir_query: "Annotate output dir holding `<sample_id>_{del,dup,insdup,ins,inv}.vcf.gz`. Trailing slashes stripped."
        remote_indir_svimasm: "BuildTruth output dir holding `<sample_id>_svimasm_*.vcf.gz`. Trailing slashes stripped."
        remote_outdir: "Where the per-SVTYPE `_training.vcf.gz` labels are written. Trailing slashes stripped. Feeds Merge (suffix='_training')."
    }

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            dipcall_beds = dipcall_beds,
            batch_size = batch_size,
            docker_image = docker_image
    }

    scatter (manifest in MakeManifests.manifests) {
        call base.SV_Integration_UltralongGetTrainingIntervals as GetIntervals {
            input:
                samples_tsv = manifest,
                remote_indir_query = sub(remote_indir_query, "/+$", ""),
                remote_indir_svimasm = sub(remote_indir_svimasm, "/+$", ""),
                remote_outdir = sub(remote_outdir, "/+$", ""),
                convert_ins_to_dup = convert_ins_to_dup,
                match_to_gaps = match_to_gaps,
                match_insdups_to_dups = match_insdups_to_dups,
                match_dups_to_insdups = match_dups_to_insdups,
                match_ins_to_dup = match_ins_to_dup,
                match_ins_to_insdup = match_ins_to_insdup,
                match_ins_to_dup_slack_bp = match_ins_to_dup_slack_bp,
                reference_fai = reference_fai,
                truvari_refdist = truvari_refdist,
                truvari_pctsize_loose = truvari_pctsize_loose,
                truvari_pctsize_strict = truvari_pctsize_strict,
                truvari_pctovl_loose = truvari_pctovl_loose,
                max_read_length = max_read_length,
                docker_image = docker_image
        }
    }

    output {
        String training_intervals_outdir = sub(remote_outdir, "/+$", "")
    }
}


# Rebuilds the `ID<TAB>dipcall_bed` samples TSV from the parallel data-table
# columns, then splits into `batch_size`-row manifests. Empty `dipcall_beds` ->
# column 2 is a `.` placeholder (ignored unless match_to_gaps=1).
task MakeManifests {
    input {
        Array[String] sample_ids
        Array[String] dipcall_beds
        Int batch_size
        String docker_image
    }

    command <<<
        set -euxo pipefail

        N_ID=$(wc -l < ~{write_lines(sample_ids)})
        N_BED=$(wc -l < ~{write_lines(dipcall_beds)})
        if [ ${N_BED} -eq 0 ]; then
            awk 'BEGIN { OFS="\t"; } { print $1, "."; }' ~{write_lines(sample_ids)} > all.tsv
        elif [ ${N_BED} -eq ${N_ID} ]; then
            paste ~{write_lines(sample_ids)} ~{write_lines(dipcall_beds)} > all.tsv
        else
            echo "ERROR: dipcall_beds has ${N_BED} rows != ${N_ID} sample_ids (pass none, or one per sample)." 1>&2
            exit 1
        fi

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
