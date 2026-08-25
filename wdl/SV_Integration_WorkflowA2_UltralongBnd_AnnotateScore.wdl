version 1.0

import "SV_Integration_UltralongAndBndAnnotate.wdl" as annotate
import "SV_Integration_Workpackage3_UltralongAndBnd.wdl" as score


# CONSOLIDATED STEP A2 of the SV integration pipeline: per-sample ULTRALONG/BND
# feature annotation (upstream `SV_Integration_UltralongAndBndAnnotate`) followed
# by per-sample XGBoost SCORING with the upstream pretrained, coverage-selected
# models (upstream `SV_Integration_Workpackage3_UltralongAndBnd`, run SCORE-ONLY:
# no variant filtering), in a SINGLE Terra submission driven directly by a data
# table.
#
# It sits between Workflow A (which writes the per-sample <sample>_ultralong.bcf /
# <sample>_bnd.bcf to /01_intrasample) and Workflow D (merge -> collapse -> concat).
# `MakeManifests` reassembles the per-sample chunk TSV the annotate/score
# containers expect from parallel data-table columns, batches the cohort into
# `batch_size`-sample VMs, and one scatter runs Annotate then Score per batch.
#
# Annotate and Score use DIFFERENT docker images, so they remain two scattered
# calls sequenced by a `done`/`upstream_signal` handshake, not one task.
#
# Outputs (per sample): /01a_annotated/<sample>_{del,ins,dup,insdup,inv,bnd}.vcf.gz
# (intermediate) and /01b_scored/<sample>_{ultralong,bnd}.bcf carrying FORMAT
# SCORE + CALIBRATION_SENSITIVITY. Point Workflow D at this remote_outdir with
# intrasample_subdir="01b_scored". Downstream filtering is left to the user.
#
workflow SV_Integration_WorkflowA2_UltralongBnd_AnnotateScore {
    input {
        # --- Per-sample columns from the Terra `sample` data table ---
        Array[String] sample_ids
        Array[String] mean_coverages
        Array[String] aligned_bais
        Array[String] aligned_bams
        Int batch_size = 20

        # --- GCS dirs (trailing slashes stripped automatically) ---
        String remote_indir
        String intrasample_subdir = "01_intrasample"
        String remote_outdir

        # --- Annotate resources ---
        File reference_fa
        File reference_fai
        File feature_extraction_py
        File tr_bed
        File segdup_bed
        File gc_content_bed
        File mei_bed_gz
        File mei_bed_tbi
        String requester_pays_project = ""
        String annotate_docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"

        # --- Score resources: pretrained per-SVTYPE x per-coverage models ---
        File del_indel_scorer_15x_pkl
        File ins_indel_scorer_15x_pkl
        File dup_indel_scorer_15x_pkl
        File insdup_indel_scorer_15x_pkl
        File inv_indel_scorer_15x_pkl
        File bnd_indel_scorer_15x_pkl
        File del_indel_calibrationScores_15x_hdf5
        File ins_indel_calibrationScores_15x_hdf5
        File dup_indel_calibrationScores_15x_hdf5
        File insdup_indel_calibrationScores_15x_hdf5
        File inv_indel_calibrationScores_15x_hdf5
        File bnd_indel_calibrationScores_15x_hdf5
        File del_indel_scorer_30x_pkl
        File ins_indel_scorer_30x_pkl
        File dup_indel_scorer_30x_pkl
        File insdup_indel_scorer_30x_pkl
        File inv_indel_scorer_30x_pkl
        File bnd_indel_scorer_30x_pkl
        File del_indel_calibrationScores_30x_hdf5
        File ins_indel_calibrationScores_30x_hdf5
        File dup_indel_calibrationScores_30x_hdf5
        File insdup_indel_calibrationScores_30x_hdf5
        File inv_indel_calibrationScores_30x_hdf5
        File bnd_indel_calibrationScores_30x_hdf5

        File scoring_python_script
        File UltralongInsdups2Ins_java
        File AddSvlenToSymbolicAlt_java
        String score_docker_image = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`. Parallel to every other per-sample array."
        mean_coverages: "e.g. `this.mean_coverage`. Per-sample mean depth: feeds the annotate coverage-bin normalization AND the 15x-vs-30x model pick (threshold 22.5x)."
        aligned_bais: "e.g. `this.02_aligned_bai`. gs:// URIs (localized in-task by the annotate step)."
        aligned_bams: "e.g. `this.02_aligned_bam`."
        batch_size: "Number of samples processed sequentially per VM."
        requester_pays_project: "Google Cloud project to bill for requester-pays buckets (the BAMs). Leave blank for non-requester-pays buckets."
        remote_indir: "Workflow A's remote_outdir (this stage reads its /01_intrasample subdir automatically). Trailing slashes are stripped."
        intrasample_subdir: "Subdir of remote_indir holding the per-sample <sample>_{ultralong,bnd}.bcf. Default '01_intrasample'; set to '' to read remote_indir directly."
        remote_outdir: "Stage outputs go to /01a_annotated (intermediate) and /01b_scored (consumed by Workflow D via intrasample_subdir=\"01b_scored\"). Trailing slashes are stripped."
    }

    String indir_root = sub(remote_indir, "/+$", "")
    String intrasample_dir = if intrasample_subdir == "" then indir_root else indir_root + "/" + intrasample_subdir
    String outdir = sub(remote_outdir, "/+$", "")
    String annotated_dir = outdir + "/01a_annotated"
    String scored_dir = outdir + "/01b_scored"

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            mean_coverages = mean_coverages,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            batch_size = batch_size,
            docker_image = annotate_docker_image
    }

    scatter (manifest in MakeManifests.manifests) {
        call annotate.SV_Integration_UltralongAndBndAnnotate as Annotate {
            input:
                chunk_tsv = manifest,
                remote_indir = intrasample_dir,
                remote_outdir = annotated_dir,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                feature_extraction_py = feature_extraction_py,
                tr_bed = tr_bed,
                segdup_bed = segdup_bed,
                gc_content_bed = gc_content_bed,
                mei_bed_gz = mei_bed_gz,
                mei_bed_tbi = mei_bed_tbi,
                requester_pays_project = requester_pays_project,
                docker_image = annotate_docker_image
        }
        call score.SV_Integration_Workpackage3_UltralongAndBnd as Score {
            input:
                sv_integration_chunk_tsv = manifest,
                remote_indir = annotated_dir,
                remote_outdir = scored_dir,
                sample_coverages_csv = MakeManifests.sample_coverages_csv,

                del_indel_scorer_15x_pkl = del_indel_scorer_15x_pkl,
                ins_indel_scorer_15x_pkl = ins_indel_scorer_15x_pkl,
                dup_indel_scorer_15x_pkl = dup_indel_scorer_15x_pkl,
                insdup_indel_scorer_15x_pkl = insdup_indel_scorer_15x_pkl,
                inv_indel_scorer_15x_pkl = inv_indel_scorer_15x_pkl,
                bnd_indel_scorer_15x_pkl = bnd_indel_scorer_15x_pkl,
                del_indel_calibrationScores_15x_hdf5 = del_indel_calibrationScores_15x_hdf5,
                ins_indel_calibrationScores_15x_hdf5 = ins_indel_calibrationScores_15x_hdf5,
                dup_indel_calibrationScores_15x_hdf5 = dup_indel_calibrationScores_15x_hdf5,
                insdup_indel_calibrationScores_15x_hdf5 = insdup_indel_calibrationScores_15x_hdf5,
                inv_indel_calibrationScores_15x_hdf5 = inv_indel_calibrationScores_15x_hdf5,
                bnd_indel_calibrationScores_15x_hdf5 = bnd_indel_calibrationScores_15x_hdf5,
                del_indel_scorer_30x_pkl = del_indel_scorer_30x_pkl,
                ins_indel_scorer_30x_pkl = ins_indel_scorer_30x_pkl,
                dup_indel_scorer_30x_pkl = dup_indel_scorer_30x_pkl,
                insdup_indel_scorer_30x_pkl = insdup_indel_scorer_30x_pkl,
                inv_indel_scorer_30x_pkl = inv_indel_scorer_30x_pkl,
                bnd_indel_scorer_30x_pkl = bnd_indel_scorer_30x_pkl,
                del_indel_calibrationScores_30x_hdf5 = del_indel_calibrationScores_30x_hdf5,
                ins_indel_calibrationScores_30x_hdf5 = ins_indel_calibrationScores_30x_hdf5,
                dup_indel_calibrationScores_30x_hdf5 = dup_indel_calibrationScores_30x_hdf5,
                insdup_indel_calibrationScores_30x_hdf5 = insdup_indel_calibrationScores_30x_hdf5,
                inv_indel_calibrationScores_30x_hdf5 = inv_indel_calibrationScores_30x_hdf5,
                bnd_indel_calibrationScores_30x_hdf5 = bnd_indel_calibrationScores_30x_hdf5,

                scoring_python_script = scoring_python_script,
                UltralongInsdups2Ins_java = UltralongInsdups2Ins_java,
                AddSvlenToSymbolicAlt_java = AddSvlenToSymbolicAlt_java,
                docker_image = score_docker_image,
                upstream_signal = Annotate.done
        }
    }

    output {
        String scored_outdir = scored_dir
    }
}


# Rebuilds, from the parallel data-table columns, the per-sample chunk TSV the
# annotate/score containers slice by position, then splits it into batches of
# `batch_size` rows (one manifest per VM). The annotate container reads (comma-
# separated after tr): field1=ID, field2=mean_coverage, field4=bai, field5=bam
# (field3 is unused, emitted as `.`); the score container reads only field1=ID
# and takes coverage from `sample_coverages.csv`. Column order MUST NOT change.
#
task MakeManifests {
    input {
        Array[String] sample_ids
        Array[String] mean_coverages
        Array[String] aligned_bais
        Array[String] aligned_bams
        Int batch_size
        String docker_image
    }

    command <<<
        set -euxo pipefail

        N_ID=$(wc -l < ~{write_lines(sample_ids)})
        N_COV=$(wc -l < ~{write_lines(mean_coverages)})
        N_BAI=$(wc -l < ~{write_lines(aligned_bais)})
        N_BAM=$(wc -l < ~{write_lines(aligned_bams)})
        for V in ${N_COV} ${N_BAI} ${N_BAM}; do
            if [ ${V} -ne ${N_ID} ]; then
                echo "ERROR: a per-sample column has ${V} rows != ${N_ID} sample_ids."
                exit 1
            fi
        done

        # Per-sample chunk TSV: ID <tab> mean_coverage <tab> . <tab> bai <tab> bam
        paste ~{write_lines(sample_ids)} ~{write_lines(mean_coverages)} \
              ~{write_lines(aligned_bais)} ~{write_lines(aligned_bams)} \
              | awk 'BEGIN { FS="\t"; OFS="\t"; } { print $1, $2, ".", $3, $4 }' > all.tsv

        # Cohort-wide coverage CSV for the score step: SAMPLE_ID,COVERAGE
        paste ~{write_lines(sample_ids)} ~{write_lines(mean_coverages)} \
              | awk 'BEGIN { FS="\t"; OFS=","; } { print $1, $2 }' > sample_coverages.csv

        split --lines=~{batch_size} --numeric-suffixes=0 --suffix-length=6 all.tsv batch_
        ls batch_* 1>&2
    >>>

    output {
        Array[File] manifests = glob("batch_*")
        File sample_coverages_csv = "sample_coverages.csv"
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}
