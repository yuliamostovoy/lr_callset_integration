version 1.0

import "SV_Integration_UltralongAnnotate.wdl" as base


# DATA-TABLE FRONT-END for `SV_Integration_UltralongAnnotate`.
#
# The upstream workflow ingests samples through a hand-built `chunk_tsv` whose
# columns are positional (`ID, mean_coverage, bai, bam`) and easy to transpose.
# This wrapper takes those four as parallel Terra data-table columns instead,
# reassembles the exact chunk TSV the underlying task expects (correct column
# order guaranteed), batches the cohort into `batch_size`-sample VMs, and scatters
# the unmodified `SV_Integration_UltralongAnnotate` workflow over the batches.
#
# The per-sample query BCF is NOT in the table: the underlying task pulls it from
# `remote_indir/<sample_id>_ultralong.bcf` by naming convention. Point `remote_indir`
# at the folder your ultralong intra-sample merge wrote to (e.g. the output dir of
# `SV_Integration_WorkflowA_UltralongBndOnly`). Everything else is a pass-through of
# the underlying workflow's inputs.
#
workflow SV_Integration_UltralongAnnotate_FromTable {
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
        Int ins2dup_bin_length = 100
        Float ins2dup_bin_coverage_ratio = 1.5
        Int convert_ins_to_dup = 1
        Int custom_breakpoint_window_bp = 500
        Int custom_min_clip_length = 200
        String custom_min_indel_length = "500,1000,5000"
        Int custom_adjacency_slack_bp = 300
        Int use_cutefc = 0
        File tr_bed
        File segdup_bed
        File gc_content_bed
        Float repeat_overlap_fraction = 0.8
        File mei_bed_gz
        File mei_bed_tbi
        Int mei_bed_tolerance = 300

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
        Int preemptible_number = 3

        Int n_cpu = 4
        Int ram_size_gb = 16
        Int disk_size_gb = 50
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`. Parallel to every other per-sample array."
        mean_coverages: "e.g. `this.mean_coverage`. Mean aligned-read coverage per sample."
        bais: "e.g. `this.aligned_bai`. gs:// path to the aligned-read BAM index."
        bams: "e.g. `this.aligned_bam`. gs:// path to the aligned-read BAM."
        batch_size: "Number of samples processed sequentially per VM."
        remote_indir: "Folder holding the per-sample query BCFs (`<sample_id>_ultralong.bcf`), i.e. the ultralong intra-sample merge output dir. Trailing slashes stripped."
        remote_outdir: "Where the per-SVTYPE annotated VCFs are written. Trailing slashes stripped. This becomes the next stage's `remote_indir_query`."
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
        call base.SV_Integration_UltralongAnnotate as Annotate {
            input:
                chunk_tsv = manifest,
                remote_indir = sub(remote_indir, "/+$", ""),
                remote_outdir = sub(remote_outdir, "/+$", ""),
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                min_mapq = min_mapq,
                ins2dup_bin_length = ins2dup_bin_length,
                ins2dup_bin_coverage_ratio = ins2dup_bin_coverage_ratio,
                convert_ins_to_dup = convert_ins_to_dup,
                custom_breakpoint_window_bp = custom_breakpoint_window_bp,
                custom_min_clip_length = custom_min_clip_length,
                custom_min_indel_length = custom_min_indel_length,
                custom_adjacency_slack_bp = custom_adjacency_slack_bp,
                use_cutefc = use_cutefc,
                tr_bed = tr_bed,
                segdup_bed = segdup_bed,
                gc_content_bed = gc_content_bed,
                repeat_overlap_fraction = repeat_overlap_fraction,
                mei_bed_gz = mei_bed_gz,
                mei_bed_tbi = mei_bed_tbi,
                mei_bed_tolerance = mei_bed_tolerance,
                docker_image = docker_image,
                preemptible_number = preemptible_number,
                n_cpu = n_cpu,
                ram_size_gb = ram_size_gb,
                disk_size_gb = disk_size_gb
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
