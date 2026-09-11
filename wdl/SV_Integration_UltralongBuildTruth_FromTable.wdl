version 1.0

import "SV_Integration_UltralongBuildTruth.wdl" as base


# DATA-TABLE FRONT-END for `SV_Integration_UltralongBuildTruth`.
#
# The upstream workflow ingests samples through a hand-built `samples_tsv` of
# `ID, dipcall_bed` (column 2 is used only when `svimasm_ins_use_gaps=1`). This
# wrapper takes those as parallel Terra data-table columns, reassembles the exact
# samples TSV, batches into `batch_size`-sample VMs, and scatters the unmodified
# `SV_Integration_UltralongBuildTruth` workflow over the batches.
#
# The per-sample svim-asm truth VCF is pulled by naming convention from
# `remote_indir_svimasm/<sample_id>_canonized.vcf.gz` — point it at the SvimAsm
# output dir. `dipcall_beds` may be left empty (`[]`) when not using gaps.
#
workflow SV_Integration_UltralongBuildTruth_FromTable {
    input {
        # --- Per-sample columns from the Terra data table ---
        Array[String] sample_ids
        Array[String] dipcall_beds = []

        Int batch_size = 20

        # --- Underlying-workflow pass-through ---
        String remote_indir_svimasm
        String remote_outdir

        Int svimasm_min_sv_length = 5000
        Int svimasm_convert_ins_to_dup = 1
        Int svimasm_ins_use_gaps = 0
        Int svimasm_ins_use_gaps_slack_bp = 200
        Float svimasm_ins_use_gaps_length_similarity = 0.9
        Int svimasm_ins_use_remap = 1
        Int svimasm_ins_remap_max_length = 2000000
        Float svimasm_ins_remap_cov_threshold = 0.8

        File reference_fa
        File reference_fai

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap:latest"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`."
        dipcall_beds: "Optional. gs:// dipcall confident-region BEDs, one per sample, e.g. `this.dipcall_bed`. Leave empty ([]) unless svimasm_ins_use_gaps=1."
        batch_size: "Number of samples processed sequentially per VM."
        remote_indir_svimasm: "SvimAsm output dir holding `<sample_id>_canonized.vcf.gz`. Trailing slashes stripped."
        remote_outdir: "Where the per-SVTYPE svim-asm truth VCFs are written. Trailing slashes stripped. Becomes the next stage's `remote_indir_svimasm`."
        svimasm_min_sv_length: "Truth-set size floor (bp). Lower from 5000 to 2000 to include 2-10 kb SVs (see TRAINING_RUNBOOK.md)."
    }

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            dipcall_beds = dipcall_beds,
            batch_size = batch_size,
            docker_image = docker_image
    }

    scatter (manifest in MakeManifests.manifests) {
        call base.SV_Integration_UltralongBuildTruth as BuildTruth {
            input:
                samples_tsv = manifest,
                remote_indir_svimasm = sub(remote_indir_svimasm, "/+$", ""),
                remote_outdir = sub(remote_outdir, "/+$", ""),
                svimasm_min_sv_length = svimasm_min_sv_length,
                svimasm_convert_ins_to_dup = svimasm_convert_ins_to_dup,
                svimasm_ins_use_gaps = svimasm_ins_use_gaps,
                svimasm_ins_use_gaps_slack_bp = svimasm_ins_use_gaps_slack_bp,
                svimasm_ins_use_gaps_length_similarity = svimasm_ins_use_gaps_length_similarity,
                svimasm_ins_use_remap = svimasm_ins_use_remap,
                svimasm_ins_remap_max_length = svimasm_ins_remap_max_length,
                svimasm_ins_remap_cov_threshold = svimasm_ins_remap_cov_threshold,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                docker_image = docker_image
        }
    }

    output {
        String truth_outdir = sub(remote_outdir, "/+$", "")
    }
}


# Rebuilds the `ID<TAB>dipcall_bed` samples TSV from the parallel data-table
# columns, then splits it into `batch_size`-row manifests. When `dipcall_beds`
# is empty, column 2 is a `.` placeholder (the underlying task ignores it unless
# svimasm_ins_use_gaps=1).
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
