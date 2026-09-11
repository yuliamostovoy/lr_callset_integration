version 1.0

import "SV_Integration_BndAnnotate.wdl" as base


# DATA-TABLE FRONT-END for `SV_Integration_BndAnnotate`.
#
# The upstream workflow ingests samples through a hand-built `chunk_tsv` whose
# columns are positional (`ID, mean_coverage, bai, bam`) and easy to transpose.
# This wrapper takes those four as parallel Terra data-table columns instead,
# reassembles the exact chunk TSV the underlying task expects (correct column
# order guaranteed), batches the cohort into `batch_size`-sample VMs, and scatters
# the unmodified `SV_Integration_BndAnnotate` workflow over the batches.
#
# The per-sample query BCF is NOT in the table: the underlying task pulls it from
# `remote_indir/<sample_id>_bnd.bcf` by naming convention. Point `remote_indir` at
# the folder your BND intra-sample merge wrote to (e.g. the output dir of
# `SV_Integration_WorkflowA_UltralongBndOnly`). Everything else is a pass-through of
# the underlying workflow's inputs.
#
# NOTE on `feature_extraction_py`: this must be the SAME file you pass to the A2
# combined annotate at scoring time, or train/score features diverge (see
# TRAINING_RUNBOOK.md "Feature parity").
#
workflow SV_Integration_BndAnnotate_FromTable {
    input {
        # --- Per-sample columns from the Terra data table (all parallel) ---
        Array[String] sample_ids
        Array[String] mean_coverages
        Array[String] bais
        Array[String] bams

        Int batch_size = 20

        # --- Underlying-workflow pass-through ---
        String remote_indir
        String remote_outdir

        File reference_fa
        File reference_fai

        String min_mapq = "0,60"
        Int custom_breakpoint_window_bp = 500
        Int custom_breakpoint_window_slack_bp = 150
        Int custom_min_clip_length = 200
        String custom_min_indel_length = "500,1000,5000"
        Int custom_adjacency_slack_bp = 300
        File feature_extraction_py
        Int use_cutefc = 0
        File tr_bed
        File segdup_bed
        File gc_content_bed

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
        Int preemptible_number = 3
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`. Parallel to every other per-sample array."
        mean_coverages: "e.g. `this.mean_coverage`. Mean aligned-read coverage per sample."
        bais: "e.g. `this.aligned_bai`. gs:// path to the aligned-read BAM index."
        bams: "e.g. `this.aligned_bam`. gs:// path to the aligned-read BAM."
        batch_size: "Number of samples processed sequentially per VM."
        remote_indir: "Folder holding the per-sample query BCFs (`<sample_id>_bnd.bcf`), i.e. the BND intra-sample merge output dir. Trailing slashes stripped."
        remote_outdir: "Where the annotated BND VCFs are written. Trailing slashes stripped. This becomes the next stage's `remote_indir_query`."
        feature_extraction_py: "MUST match the file passed to A2 combined annotate at scoring time (see TRAINING_RUNBOOK.md feature-parity note)."
    }

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            mean_coverages = mean_coverages,
            bais = bais,
            bams = bams,
            batch_size = batch_size,
            docker_image = docker_image
    }

    scatter (manifest in MakeManifests.manifests) {
        call base.SV_Integration_BndAnnotate as Annotate {
            input:
                chunk_tsv = manifest,
                remote_indir = sub(remote_indir, "/+$", ""),
                remote_outdir = sub(remote_outdir, "/+$", ""),
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                min_mapq = min_mapq,
                custom_breakpoint_window_bp = custom_breakpoint_window_bp,
                custom_breakpoint_window_slack_bp = custom_breakpoint_window_slack_bp,
                custom_min_clip_length = custom_min_clip_length,
                custom_min_indel_length = custom_min_indel_length,
                custom_adjacency_slack_bp = custom_adjacency_slack_bp,
                feature_extraction_py = feature_extraction_py,
                use_cutefc = use_cutefc,
                tr_bed = tr_bed,
                segdup_bed = segdup_bed,
                gc_content_bed = gc_content_bed,
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
    }

    output {
        String annotated_outdir = sub(remote_outdir, "/+$", "")
    }
}


# Rebuilds the `ID<TAB>mean_coverage<TAB>bai<TAB>bam` chunk TSV from the parallel
# data-table columns (this order is exactly what the underlying task slices by
# position), then splits it into `batch_size`-row manifests (one per VM).
task MakeManifests {
    input {
        Array[String] sample_ids
        Array[String] mean_coverages
        Array[String] bais
        Array[String] bams
        Int batch_size
        String docker_image
    }

    command <<<
        set -euxo pipefail

        N_ID=$(wc -l < ~{write_lines(sample_ids)})
        N_COV=$(wc -l < ~{write_lines(mean_coverages)})
        N_BAI=$(wc -l < ~{write_lines(bais)})
        N_BAM=$(wc -l < ~{write_lines(bams)})
        for V in ${N_COV} ${N_BAI} ${N_BAM}; do
            if [ ${V} -ne ${N_ID} ]; then
                echo "ERROR: a per-sample column has ${V} rows != ${N_ID} sample_ids." 1>&2
                exit 1
            fi
        done

        paste ~{write_lines(sample_ids)} ~{write_lines(mean_coverages)} \
              ~{write_lines(bais)} ~{write_lines(bams)} > all.tsv

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
