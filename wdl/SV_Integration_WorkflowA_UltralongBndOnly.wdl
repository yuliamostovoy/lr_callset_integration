version 1.0

import "SV_Integration_Workpackage1_ExtractUltralongBnd2kb.wdl" as extract


# CONSOLIDATED STEP A (ULTRALONG/BND-ONLY variant) of the SV integration
# pipeline: per-sample intra-sample merge of the 2kb+ ULTRALONG and BND callsets
# ONLY — no short/main branch, no Kanpig re-genotyping, no scoring — in a SINGLE
# Terra submission driven directly by a data table.
#
# This is the lightweight sibling of `SV_Integration_WorkflowA_Intrasample_Scoring`.
# Where full Workflow A runs WP1 (merge + Kanpig) then WP2 (short-SV scoring +
# chunk split), this one runs ONLY `SV_Integration_Workpackage1_ExtractUltralongBnd2kb`,
# which stops before the shorter SVs, BAM localization, and Kanpig. Use it when you
# only need the per-sample ultralong/BND calls — e.g. to (re)generate the query
# inputs for the ultralong/BND model TRAINING pipeline or for Workflow A2 scoring.
#
# PAV: set `has_pav = false` (the default here) to merge ONLY pbsv + sniffles and
# leave the `pav_*` inputs empty — the PAV-free merge used for a cohort without PAV.
# The output schema still carries INFO/SUPP_PAV, set to zero, so downstream
# annotate/score never break on a missing field (you still drop SUPP_PAV from the
# model feature lists, since a constant column is useless — see TRAINING_RUNBOOK.md).
#
# `MakeManifests` reassembles the exact per-sample chunk TSV the WP1 container
# expects from parallel data-table columns (inserting `.` placeholders for the
# sex/bai/bam columns the extraction never reads), batches the cohort into
# `batch_size`-sample VMs, and one scatter runs the extraction per batch.
#
# Outputs (per sample, written directly to `remote_outdir`):
#   <sample>_ultralong.bcf(.csi)   — 2kb+ ultralong intra-sample merge
#   <sample>_bnd.bcf(.csi)         — BND intra-sample merge
# Point the training pipeline / Workflow A2 at this `remote_outdir` (i.e. set their
# remote_indir + intrasample_subdir so they read these files).
#
workflow SV_Integration_WorkflowA_UltralongBndOnly {
    input {
        # --- Per-sample columns from the Terra `sample` data table ---
        Array[String] sample_ids
        Array[String] pbsv_tbis
        Array[String] pbsv_vcfs
        Array[String] sniffles_tbis
        Array[String] sniffles_vcfs
        Array[String] pav_beds = []
        Array[String] pav_tbis = []
        Array[String] pav_vcfs = []

        Boolean has_pav = false
        Int batch_size = 20

        # --- GCS output dir (trailing slashes stripped automatically) ---
        String remote_outdir

        # --- Extraction (WP1) parameters ---
        String region = "all"
        String requester_pays_project = ""
        Int min_ultralong_sv_length = 2000
        Int ultralong_collapse_mode = 0
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File reference_agp
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`. Parallel to every other per-sample array."
        pbsv_tbis: "e.g. `this.pbsv_tbi`. gs:// URIs (localized in-task by the extraction)."
        pbsv_vcfs: "e.g. `this.pbsv_vcf`."
        sniffles_tbis: "e.g. `this.sniffles_tbi`."
        sniffles_vcfs: "e.g. `this.sniffles_vcf`."
        pav_beds: "Only needed when has_pav=true. e.g. `this.pav_bed`."
        pav_tbis: "Only needed when has_pav=true."
        pav_vcfs: "Only needed when has_pav=true."
        has_pav: "false (default) = merge only pbsv+sniffles and leave pav_* empty (PAV-free). true = also merge PAV (all three pav_* arrays required)."
        batch_size: "Number of samples processed sequentially per VM."
        remote_outdir: "Where the per-sample <sample>_{ultralong,bnd}.bcf are written. Trailing slashes are stripped. Point training/A2 here."
        requester_pays_project: "Google Cloud project to bill for requester-pays buckets (the caller VCFs). Leave blank for non-requester-pays buckets."
        min_ultralong_sv_length: "Lower size bound (bp) for the ultralong branch. Default 2000."
        ultralong_collapse_mode: "Passed through to the WP1 extraction."
    }

    String outdir = sub(remote_outdir, "/+$", "")

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            pbsv_tbis = pbsv_tbis,
            pbsv_vcfs = pbsv_vcfs,
            sniffles_tbis = sniffles_tbis,
            sniffles_vcfs = sniffles_vcfs,
            pav_beds = pav_beds,
            pav_tbis = pav_tbis,
            pav_vcfs = pav_vcfs,
            has_pav = has_pav,
            batch_size = batch_size,
            docker_image = docker_image
    }

    scatter (manifest in MakeManifests.manifests) {
        call extract.SV_Integration_Workpackage1_UltralongBnd2kb as Extract {
            input:
                sv_integration_chunk_tsv = manifest,
                has_pav = has_pav,
                region = region,
                remote_outdir = outdir,
                requester_pays_project = requester_pays_project,
                min_ultralong_sv_length = min_ultralong_sv_length,
                ultralong_collapse_mode = ultralong_collapse_mode,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                standard_chromosomes_bed = standard_chromosomes_bed,
                reference_agp = reference_agp,
                docker_image = docker_image
        }
    }

    output {
        String intrasample_outdir = outdir
    }
}


# Rebuilds, from the parallel data-table columns, the per-sample chunk TSV the WP1
# extraction container slices by position, then splits it into batches of
# `batch_size` rows (one manifest per VM).
#
# The extraction reads (comma-separated after `tr`): field1=sample_id, then the
# caller columns. Fields 2-4 (sample_sex, aligned_bai, aligned_bam in the full
# `sv_integration_hg38` schema) are NEVER read by this extraction, so we emit `.`
# placeholders to keep this workflow's inputs minimal while preserving the column
# positions the container expects:
#   has_pav=true  (11 cols): id . . . pav_bed pav_tbi pav_vcf pbsv_tbi pbsv_vcf sniffles_tbi sniffles_vcf
#   has_pav=false ( 8 cols): id . . . pbsv_tbi pbsv_vcf sniffles_tbi sniffles_vcf
# Column order MUST NOT change.
#
task MakeManifests {
    input {
        Array[String] sample_ids
        Array[String] pbsv_tbis
        Array[String] pbsv_vcfs
        Array[String] sniffles_tbis
        Array[String] sniffles_vcfs
        Array[String] pav_beds
        Array[String] pav_tbis
        Array[String] pav_vcfs
        Boolean has_pav
        Int batch_size
        String docker_image
    }

    command <<<
        set -euxo pipefail

        N_ID=$(wc -l < ~{write_lines(sample_ids)})
        N_PBSV_TBI=$(wc -l < ~{write_lines(pbsv_tbis)})
        N_PBSV_VCF=$(wc -l < ~{write_lines(pbsv_vcfs)})
        N_SNIF_TBI=$(wc -l < ~{write_lines(sniffles_tbis)})
        N_SNIF_VCF=$(wc -l < ~{write_lines(sniffles_vcfs)})
        for V in ${N_PBSV_TBI} ${N_PBSV_VCF} ${N_SNIF_TBI} ${N_SNIF_VCF}; do
            if [ ${V} -ne ${N_ID} ]; then
                echo "ERROR: a per-sample column has ${V} rows != ${N_ID} sample_ids."
                exit 1
            fi
        done

        if [ ~{true="1" false="0" has_pav} -eq 1 ]; then
            N_PAV_BED=$(wc -l < ~{write_lines(pav_beds)})
            N_PAV_TBI=$(wc -l < ~{write_lines(pav_tbis)})
            N_PAV_VCF=$(wc -l < ~{write_lines(pav_vcfs)})
            for V in ${N_PAV_BED} ${N_PAV_TBI} ${N_PAV_VCF}; do
                if [ ${V} -ne ${N_ID} ]; then
                    echo "ERROR: has_pav=true but a PAV column has ${V} rows != ${N_ID} sample_ids."
                    exit 1
                fi
            done
            paste ~{write_lines(sample_ids)} \
                  ~{write_lines(pav_beds)} ~{write_lines(pav_tbis)} ~{write_lines(pav_vcfs)} \
                  ~{write_lines(pbsv_tbis)} ~{write_lines(pbsv_vcfs)} \
                  ~{write_lines(sniffles_tbis)} ~{write_lines(sniffles_vcfs)} \
                  | awk 'BEGIN { FS="\t"; OFS="\t"; } { print $1, ".", ".", ".", $2, $3, $4, $5, $6, $7, $8 }' > all.tsv
        else
            paste ~{write_lines(sample_ids)} \
                  ~{write_lines(pbsv_tbis)} ~{write_lines(pbsv_vcfs)} \
                  ~{write_lines(sniffles_tbis)} ~{write_lines(sniffles_vcfs)} \
                  | awk 'BEGIN { FS="\t"; OFS="\t"; } { print $1, ".", ".", ".", $2, $3, $4, $5 }' > all.tsv
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
