version 1.0

import "SV_Integration_Workpackage1_intrasample_merge.wdl" as wp1
import "SV_Integration_Workpackage2_Main_scoring.wdl" as wp2


# CONSOLIDATED STEP A of the SV integration pipeline: per-sample intra-sample
# merge + Kanpig (WP1) followed by per-sample XGBoost scoring + chunk split
# (WP2), in a SINGLE Terra submission driven directly by a data table.
#
# It replaces the old "run WP1, hand-build a chunk TSV, run WP2" sequence. The
# sample columns are bound straight from the Terra `sample` table as parallel
# Array[String] inputs (the PAV pattern); `MakeManifests` reassembles the exact
# per-sample TSV that the WP1/WP2 containers expect, batches the cohort into
# `batch_size`-sample VMs, and one scatter runs WP1 then WP2 per batch.
#
# WP1 and WP2 use DIFFERENT docker images, so they remain two scattered task
# calls sequenced by a `done`/`upstream_signal` handshake, not one task.
#
workflow SV_Integration_WorkflowA_Intrasample_Scoring {
    input {
        # --- Per-sample columns from the Terra `sample` data table ---
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bais
        Array[String] aligned_bams
        Array[String] pbsv_tbis
        Array[String] pbsv_vcfs
        Array[String] sniffles_tbis
        Array[String] sniffles_vcfs
        Array[String] pav_beds = []
        Array[String] pav_tbis = []
        Array[String] pav_vcfs = []

        Boolean has_pav = true
        Int batch_size = 20

        # --- GCS output dirs (no final slash) ---
        String intrasample_remote_outdir
        String scoring_remote_outdir

        # --- WP1 parameters ---
        String region = "all"
        String requester_pays_project = ""
        Int min_sv_length = 20
        Int max_sv_length = 10000
        String kanpig_params_singlesample = "--neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
        Int ultralong_collapse_mode = 0
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File autosomes_bed
        File reference_agp
        File ploidy_bed_female
        File ploidy_bed_male
        String wp1_docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"

        # --- WP2 parameters ---
        File split_for_bcftools_merge_csv
        String filter_string = "none"
        Array[String] annotations = ["KS_1","KS_2","SQ","GQ","DP","AD_NON_ALT","AD_ALL","GT_COUNT","SUPP_PAV","SUPP_SNIFFLES","SUPP_PBSV","SVLEN"]
        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        String wp2_docker_image = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`. Parallel to every other per-sample array."
        sample_sexes: "e.g. `this.sex`. 'M' selects the male ploidy BED in Kanpig; anything else selects female. The table MUST contain this column."
        aligned_bais: "e.g. `this.02_aligned_bai`. gs:// URIs (streamed in-task)."
        aligned_bams: "e.g. `this.02_aligned_bam`."
        pbsv_tbis: "e.g. `this.03_pbsv_tbi`."
        pbsv_vcfs: "e.g. `this.03_pbsv_vcf`."
        sniffles_tbis: "e.g. `this.03_sniffles_tbi`."
        sniffles_vcfs: "e.g. `this.03_sniffles_vcf`."
        pav_beds: "e.g. `this.03_pav_bed`. Required only when has_pav=true (leave empty otherwise)."
        pav_tbis: "e.g. `this.03_pav_tbi`. Required only when has_pav=true."
        pav_vcfs: "e.g. `this.03_pav_vcf`. Required only when has_pav=true."
        has_pav: "true: 11-column manifest incl. PAV. false: 8-column manifest (pbsv+sniffles only)."
        batch_size: "Number of samples processed sequentially per VM."
        intrasample_remote_outdir: "Without final slash. WP1 per-sample outputs (kanpig/training/bnd/ultralong) land here; WP2 reads its input from here."
        scoring_remote_outdir: "Without final slash. WP2 writes chunk_<i>/<sample>.bcf and <sample>.done here; consumed by Workflow B."
        split_for_bcftools_merge_csv: "The interval partition CSV; chunk id == 0-based line number. Same file the merge/collapse stages use."
        training_resource_bed: "Shared by WP1 and WP2."
    }

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            sample_sexes = sample_sexes,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            pbsv_tbis = pbsv_tbis,
            pbsv_vcfs = pbsv_vcfs,
            sniffles_tbis = sniffles_tbis,
            sniffles_vcfs = sniffles_vcfs,
            pav_beds = pav_beds,
            pav_tbis = pav_tbis,
            pav_vcfs = pav_vcfs,
            has_pav = has_pav,
            batch_size = batch_size,
            docker_image = wp1_docker_image
    }

    scatter (manifest in MakeManifests.manifests) {
        call wp1.Impl as Intrasample {
            input:
                sv_integration_chunk_tsv = manifest,
                has_pav = has_pav,
                region = region,
                remote_outdir = intrasample_remote_outdir,
                requester_pays_project = requester_pays_project,
                min_sv_length = min_sv_length,
                max_sv_length = max_sv_length,
                kanpig_params_singlesample = kanpig_params_singlesample,
                ultralong_collapse_mode = ultralong_collapse_mode,
                training_resource_vcf_gz = training_resource_vcf_gz,
                training_resource_tbi = training_resource_tbi,
                training_resource_bed = training_resource_bed,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                standard_chromosomes_bed = standard_chromosomes_bed,
                autosomes_bed = autosomes_bed,
                reference_agp = reference_agp,
                ploidy_bed_female = ploidy_bed_female,
                ploidy_bed_male = ploidy_bed_male,
                docker_image = wp1_docker_image
        }
        call wp2.Impl as Scoring {
            input:
                sv_integration_chunk_tsv = manifest,
                split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
                filter_string = filter_string,
                remote_indir = intrasample_remote_outdir,
                remote_outdir = scoring_remote_outdir,
                training_resource_bed = training_resource_bed,
                annotations = annotations,
                training_python_script = training_python_script,
                scoring_python_script = scoring_python_script,
                hyperparameters_json = hyperparameters_json,
                docker_image = wp2_docker_image,
                upstream_signal = Intrasample.done
        }
    }

    output {
    }
}


# Rebuilds, from the parallel data-table columns, the tab-separated per-sample
# TSV that the WP1/WP2 containers slice by position, then splits it into batches
# of `batch_size` rows (one manifest per VM). Column order is fixed by the
# container `cut` positions and MUST NOT change:
#   has_pav=true  (11): sample_id sex aligned_bai aligned_bam pav_bed pav_tbi pav_vcf pbsv_tbi pbsv_vcf sniffles_tbi sniffles_vcf
#   has_pav=false  (8): sample_id sex aligned_bai aligned_bam pbsv_tbi pbsv_vcf sniffles_tbi sniffles_vcf
#
task MakeManifests {
    input {
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bais
        Array[String] aligned_bams
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
        N_SEX=$(wc -l < ~{write_lines(sample_sexes)})
        N_BAI=$(wc -l < ~{write_lines(aligned_bais)})
        N_BAM=$(wc -l < ~{write_lines(aligned_bams)})
        N_PBSV_TBI=$(wc -l < ~{write_lines(pbsv_tbis)})
        N_PBSV_VCF=$(wc -l < ~{write_lines(pbsv_vcfs)})
        N_SNIF_TBI=$(wc -l < ~{write_lines(sniffles_tbis)})
        N_SNIF_VCF=$(wc -l < ~{write_lines(sniffles_vcfs)})
        for V in ${N_SEX} ${N_BAI} ${N_BAM} ${N_PBSV_TBI} ${N_PBSV_VCF} ${N_SNIF_TBI} ${N_SNIF_VCF}; do
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
            paste ~{write_lines(sample_ids)} ~{write_lines(sample_sexes)} \
                  ~{write_lines(aligned_bais)} ~{write_lines(aligned_bams)} \
                  ~{write_lines(pav_beds)} ~{write_lines(pav_tbis)} ~{write_lines(pav_vcfs)} \
                  ~{write_lines(pbsv_tbis)} ~{write_lines(pbsv_vcfs)} \
                  ~{write_lines(sniffles_tbis)} ~{write_lines(sniffles_vcfs)} > all.tsv
        else
            paste ~{write_lines(sample_ids)} ~{write_lines(sample_sexes)} \
                  ~{write_lines(aligned_bais)} ~{write_lines(aligned_bams)} \
                  ~{write_lines(pbsv_tbis)} ~{write_lines(pbsv_vcfs)} \
                  ~{write_lines(sniffles_tbis)} ~{write_lines(sniffles_vcfs)} > all.tsv
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
