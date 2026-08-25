version 1.0


# A version of `SV_Integration_Workpackage3.wdl` that takes in input the 
# ultralong calls annotated by `SV_Integration_UltralongAndBndAnnotate.wdl`,
# filters them using global models, and produces output compatible with 
# `SV_Integration_Workpackage12.wdl`.
#
# Fork modification (SCORE-ONLY): NO variant filtering is applied. One scored
# VCF is produced for each sample and each class (ultralong, BND), carrying the
# XGBoost SCORE and CALIBRATION_SENSITIVITY in FORMAT. Downstream filtering is
# left to the user, exactly like the short-SV branch (WP2 Main_scoring).
#
# Remark: INSDUPs are converted back to INS, to preserve as much information as 
# possible downstream. The breakpoints of the original INSDUPs are saved in
# INFO.
#
# Remark: the entire ultralong pipeline supports symbolic <INS>, which are never
# removed.
#
workflow SV_Integration_Workpackage3_UltralongAndBnd {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir

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

        File sample_coverages_csv

        Array[String] annotations_interval = [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
                                               "START_TR","END_TR","START_SD","END_SD","START_GC","END_GC","INTERVAL_TR","INTERVAL_SD","INTERVAL_GC",
                                               "BIN_BEFORE_COVERAGE_0","BIN_LEFT_COVERAGE_0","BIN_1_COVERAGE_0","BIN_2_COVERAGE_0","BIN_3_COVERAGE_0","BIN_4_COVERAGE_0","BIN_5_COVERAGE_0","BIN_6_COVERAGE_0","BIN_7_COVERAGE_0","BIN_8_COVERAGE_0","BIN_9_COVERAGE_0","BIN_10_COVERAGE_0","BIN_RIGHT_COVERAGE_0","BIN_AFTER_COVERAGE_0",
                                               "BIN_LEFT_MAPQ","BIN_RIGHT_MAPQ","BIN_LEFT_SECONDARY","BIN_RIGHT_SECONDARY",
                                               "LL_0","LR_0","RL_0","RR_0",
                                               "LL_RL_1_0","LL_RL_2_0","LL_RL_3_0","LL_RL_4_0",
                                               "LL_RR_1_0","LL_RR_2_0","LL_RR_3_0","LL_RR_4_0",
                                               "LR_RL_1_0","LR_RL_2_0","LR_RL_3_0","LR_RL_4_0",
                                               "LR_RR_1_0","LR_RR_2_0","LR_RR_3_0","LR_RR_4_0",
                                               "L_INS_0_500","L_DEL_START_0_500","L_DEL_END_0_500",
                                               "R_INS_0_500","R_DEL_START_0_500","R_DEL_END_0_500",
                                               "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK" 
                                             ]
        Array[String] annotations_point = [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
                                            "START_TR","START_SD","START_GC",
                                            "BIN_POS_0",
                                            "BIN_POINT_MAPQ","BIN_POINT_SECONDARY",
                                            "PL_0","PR_0",
                                            "PL_PL_1_0","PL_PL_2_0","PL_PL_3_0","PL_PL_4_0",
                                            "PL_PR_1_0","PL_PR_2_0","PL_PR_3_0","PL_PR_4_0",
                                            "PR_PR_1_0","PR_PR_2_0","PR_PR_3_0","PR_PR_4_0",
                                            "P_INS_0_500","P_DEL_START_0_500","P_DEL_END_0_500",
                                            "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"
                                          ]
        Array[String] annotations_bnd = [ "GT_COUNT","SUPP_SNIFFLES","SUPP_PBSV",
                                          "START_TR","END_TR","START_SD","END_SD","START_GC","END_GC",
                                          "BIN_POS_0_0","BIN_POS_1_0","BIN_POS_2_0","BIN_POS_3_0",
                                          "BIN_POINT_MAPQ_0","BIN_POINT_MAPQ_1","BIN_POINT_MAPQ_2","BIN_POINT_MAPQ_3",
                                          "BIN_POINT_SECONDARY_0","BIN_POINT_SECONDARY_1","BIN_POINT_SECONDARY_2","BIN_POINT_SECONDARY_3",
                                          "C0R_0","C1L_0","C2R_0","C3L_0",
                                          "C0R_C2R_1_0","C0R_C2R_2_0","C0R_C2R_3_0","C0R_C2R_4_0",
                                          "C0R_C3L_1_0","C0R_C3L_2_0","C0R_C3L_3_0","C0R_C3L_4_0",
                                          "C1L_C2R_1_0","C1L_C2R_2_0","C1L_C2R_3_0","C1L_C2R_4_0",
                                          "C1L_C3L_1_0","C1L_C3L_2_0","C1L_C3L_3_0","C1L_C3L_4_0",
                                          "C0_INS_0_500","C0_DEL_START_0_500","C0_DEL_END_0_500",
                                          "C1_INS_0_500","C1_DEL_START_0_500","C1_DEL_END_0_500",
                                          "C2_INS_0_500","C2_DEL_START_0_500","C2_DEL_END_0_500",
                                          "C3_INS_0_500","C3_DEL_START_0_500","C3_DEL_END_0_500",
                                          "FEX_DEPTH_RATIO","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"
                                        ]
        File scoring_python_script

        String docker_image = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
        File UltralongInsdups2Ins_java
        File AddSvlenToSymbolicAlt_java

        String upstream_signal = ""
    }
    parameter_meta {
        remote_indir: "Trailing slashes are stripped automatically."
        remote_outdir: "Trailing slashes are stripped automatically."
        sample_coverages_csv: "One line per sample, with columns: `SAMPLE_ID,COVERAGE`. Used to select the appropriate model for each sample."
        upstream_signal: "Ordering-only handshake for orchestrator workflows: set to the `done` output of the upstream annotate step so scoring starts only after it. Ignored by standalone runs."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = sub(remote_indir, "/+$", ""),
            remote_outdir = sub(remote_outdir, "/+$", ""),

            del_indel_scorer_15x_pkl = del_indel_scorer_15x_pkl,
            ins_indel_scorer_15x_pkl = ins_indel_scorer_15x_pkl,
            dup_indel_scorer_15x_pkl = dup_indel_scorer_15x_pkl,
            insdup_indel_scorer_15x_pkl = insdup_indel_scorer_15x_pkl,
            inv_indel_scorer_15x_pkl = inv_indel_scorer_15x_pkl,
            bnd_indel_scorer_15x_pkl = bnd_indel_scorer_15x_pkl,

            del_indel_scorer_30x_pkl = del_indel_scorer_30x_pkl,
            ins_indel_scorer_30x_pkl = ins_indel_scorer_30x_pkl,
            dup_indel_scorer_30x_pkl = dup_indel_scorer_30x_pkl,
            insdup_indel_scorer_30x_pkl = insdup_indel_scorer_30x_pkl,
            inv_indel_scorer_30x_pkl = inv_indel_scorer_30x_pkl,
            bnd_indel_scorer_30x_pkl = bnd_indel_scorer_30x_pkl,

            del_indel_calibrationScores_15x_hdf5 = del_indel_calibrationScores_15x_hdf5,
            ins_indel_calibrationScores_15x_hdf5 = ins_indel_calibrationScores_15x_hdf5,
            dup_indel_calibrationScores_15x_hdf5 = dup_indel_calibrationScores_15x_hdf5,
            insdup_indel_calibrationScores_15x_hdf5 = insdup_indel_calibrationScores_15x_hdf5,
            inv_indel_calibrationScores_15x_hdf5 = inv_indel_calibrationScores_15x_hdf5,
            bnd_indel_calibrationScores_15x_hdf5 = bnd_indel_calibrationScores_15x_hdf5,

            del_indel_calibrationScores_30x_hdf5 = del_indel_calibrationScores_30x_hdf5,
            ins_indel_calibrationScores_30x_hdf5 = ins_indel_calibrationScores_30x_hdf5,
            dup_indel_calibrationScores_30x_hdf5 = dup_indel_calibrationScores_30x_hdf5,
            insdup_indel_calibrationScores_30x_hdf5 = insdup_indel_calibrationScores_30x_hdf5,
            inv_indel_calibrationScores_30x_hdf5 = inv_indel_calibrationScores_30x_hdf5,
            bnd_indel_calibrationScores_30x_hdf5 = bnd_indel_calibrationScores_30x_hdf5,

            sample_coverages_csv = sample_coverages_csv,

            annotations_interval = annotations_interval,
            annotations_point = annotations_point,
            annotations_bnd = annotations_bnd,
            scoring_python_script = scoring_python_script,

            docker_image = docker_image,
            UltralongInsdups2Ins_java = UltralongInsdups2Ins_java,
            AddSvlenToSymbolicAlt_java = AddSvlenToSymbolicAlt_java,
            upstream_signal = upstream_signal
    }

    output {
    }
}


# Remark: we use gsutil instead of gcloud since we found the latter to have
# issues in practice (maybe the gcloud version in the docker is not up to
# date?). 
#
# Memory bottlenecks (measured on a 4GB VM):
#
# ScoreVariantAnnotations               200 MB
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir

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

        File sample_coverages_csv

        Array[String] annotations_interval
        Array[String] annotations_point
        Array[String] annotations_bnd
        File scoring_python_script

        String docker_image
        File UltralongInsdups2Ins_java
        File AddSvlenToSymbolicAlt_java
        String upstream_signal = ""
        Int n_cpu = 2
        Int ram_size_gb = 4
        Int disk_size_gb = 20
        Int preemptible_number = 4
    }
    parameter_meta {
    }
    
    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export GATK_LOCAL_JAR="/root/gatk.jar"
        RAM_PER_THREAD_MB=$(( ~{ram_size_gb} * 1024 - 500 ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local REMOTE_DIR=$2
            
            gsutil cp ${REMOTE_DIR}/${SAMPLE_ID}_del.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_ins.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_dup.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_insdup.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_inv.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_bnd.vcf.'gz*' .
        }

        
        # Deletes all files and directories related to the sample.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -rf ./${SAMPLE_ID}_*
        }


        # Copies the following fields from INFO to FORMAT, so that they are
        # preserved by the inter-sample merge downstream:
        #
        # SUPP_*, SCORE, CALIBRATION_SENSITIVITY
        #
        # Remark: since we don't re-genotype the cohort-level VCF downstream, at
        # this stage we should copy other INFO fields that may be useful later,
        # e.g. the number of supporting reads.
        #
        # @param 2 A VCF where all IDs are distinct. This is guaranteed by the
        # workpackages upstream.
        #
        function CopyInfoToFormat() {
            local SAMPLE_ID=$1
            local SVTYPE=$2
            
            echo '##FORMAT=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            echo '##FORMAT=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            echo '##FORMAT=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            echo '##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="Score according to the XGBoost model">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            echo '##FORMAT=<ID=CALIBRATION_SENSITIVITY,Number=1,Type=Float,Description="Calibration sensitivity according to the model applied by ScoreVariantAnnotations">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%SUPP_PBSV\t%SUPP_SNIFFLES\t%SUPP_PAV\t%SCORE\t%CALIBRATION_SENSITIVITY\n' ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz | bgzip -c > ${SAMPLE_ID}_${SVTYPE}_format.tsv.gz
            tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_${SVTYPE}_format.tsv.gz
            local COLUMNS='CHROM,POS,REF,ALT,~ID,FORMAT/SUPP_PBSV,FORMAT/SUPP_SNIFFLES,FORMAT/SUPP_PAV,FORMAT/SCORE,FORMAT/CALIBRATION_SENSITIVITY'
            bcftools annotate --threads ${N_THREADS} --header-lines ${SAMPLE_ID}_${SVTYPE}_header.txt --annotations ${SAMPLE_ID}_${SVTYPE}_format.tsv.gz --columns ${COLUMNS} --output-type z ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz --output ${SAMPLE_ID}_${SVTYPE}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz* ; mv ${SAMPLE_ID}_${SVTYPE}_out.vcf.gz ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz | head -n 1 || echo "0") 1>&2
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_${SVTYPE}_header.txt ${SAMPLE_ID}_${SVTYPE}_format.tsv.gz*
        }


        # Assumes that `CopyInfoToFormat()` has already been executed.
        #
        function PrintDebugInformation() {
            local SAMPLE_ID=$1
            local SVTYPE=$2
            
            rm -f ${SAMPLE_ID}_${SVTYPE}_xgboost.csv
            local N_RECORDS_BEFORE_FILTERING=$(bcftools index --nrecords ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz)
            for THRESHOLD in 0.7 0.8 0.9 0.95 ; do
                local N_RECORDS_AFTER_FILTERING=$( bcftools query --format '%ID\n' --include "FORMAT/CALIBRATION_SENSITIVITY<=${THRESHOLD}" ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz | wc -l )
                if [ ${N_RECORDS_BEFORE_FILTERING} -eq 0 ]; then
                    local PERCENT=0
                else
                    local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_FILTERING} / ${N_RECORDS_BEFORE_FILTERING}" | bc )
                fi
                echo "${N_RECORDS_AFTER_FILTERING},${N_RECORDS_BEFORE_FILTERING},${PERCENT},Number of records with CALIBRATION_SENSITIVITY<=${THRESHOLD}" >> ${SAMPLE_ID}_${SVTYPE}_xgboost.csv
            done
        }


        function Score() {
            local SAMPLE_ID=$1
            local COVERAGE=$2
            local SVTYPE=$3
            local RAM_PER_THREAD_MB=$4

            # Scoring
            if [ ${SVTYPE} = "ins" ]; then
                local ANNOTATIONS="~{sep=" -A " annotations_point}"
            elif [ ${SVTYPE} = "bnd" ]; then
                local ANNOTATIONS="~{sep=" -A " annotations_bnd}"
            else
                local ANNOTATIONS="~{sep=" -A " annotations_interval}"
            fi
            gatk --java-options "-Xmx${RAM_PER_THREAD_MB}m" ScoreVariantAnnotations -V ${SAMPLE_ID}_${SVTYPE}.vcf.gz -O ${SAMPLE_ID}_${SVTYPE}_score -A ${ANNOTATIONS} --model-prefix ${COVERAGE}.${SVTYPE} --model-backend PYTHON_SCRIPT --python-script ~{scoring_python_script} --mode INDEL --mnp-type INDEL --ignore-all-filters --verbosity DEBUG
            CopyInfoToFormat ${SAMPLE_ID} ${SVTYPE}
            PrintDebugInformation ${SAMPLE_ID} ${SVTYPE}

            # Converting INSDUP to INS
            if [ ${SVTYPE} = "insdup" ]; then
                java UltralongInsdups2Ins ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz | bcftools sort --max-mem ${RAM_PER_THREAD_MB}M --output-type z --output ${SAMPLE_ID}_${SVTYPE}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz*
                mv ${SAMPLE_ID}_${SVTYPE}_out.vcf.gz ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz
            fi

            # Score-only (fork modification): NO variant filtering. Emit the
            # scored BCF as-is; downstream filtering on FORMAT/CALIBRATION_SENSITIVITY
            # is left to the user (matching the short-SV branch). The per-SVTYPE
            # threshold-pass counts written by PrintDebugInformation are retained
            # for diagnostics only and do not remove any records.
            bcftools view --threads ${N_THREADS} --output-type b ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz --output ${SAMPLE_ID}_${SVTYPE}.bcf
            bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_${SVTYPE}.bcf
            rm -f ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz*
            cat ${SAMPLE_ID}_${SVTYPE}_xgboost.csv 1>&2
        }




        # ---------------------------- Main program ----------------------------

        # Compiling input scripts (ugly, but we don't want to change the
        # docker).
        mv ~{UltralongInsdups2Ins_java} UltralongInsdups2Ins.java
        mv ~{AddSvlenToSymbolicAlt_java} AddSvlenToSymbolicAlt.java
        javac *.java

        # Enforcing a consistent naming scheme on all model files
        mv ~{del_indel_scorer_15x_pkl} 15x.del.indel.scorer.pkl
        mv ~{ins_indel_scorer_15x_pkl} 15x.ins.indel.scorer.pkl
        mv ~{dup_indel_scorer_15x_pkl} 15x.dup.indel.scorer.pkl
        mv ~{insdup_indel_scorer_15x_pkl} 15x.insdup.indel.scorer.pkl
        mv ~{inv_indel_scorer_15x_pkl} 15x.inv.indel.scorer.pkl
        mv ~{bnd_indel_scorer_15x_pkl} 15x.bnd.indel.scorer.pkl
        mv ~{del_indel_calibrationScores_15x_hdf5} 15x.del.indel.calibrationScores.hdf5
        mv ~{ins_indel_calibrationScores_15x_hdf5} 15x.ins.indel.calibrationScores.hdf5
        mv ~{dup_indel_calibrationScores_15x_hdf5} 15x.dup.indel.calibrationScores.hdf5
        mv ~{insdup_indel_calibrationScores_15x_hdf5} 15x.insdup.indel.calibrationScores.hdf5
        mv ~{inv_indel_calibrationScores_15x_hdf5} 15x.inv.indel.calibrationScores.hdf5
        mv ~{bnd_indel_calibrationScores_15x_hdf5} 15x.bnd.indel.calibrationScores.hdf5

        mv ~{del_indel_scorer_30x_pkl} 30x.del.indel.scorer.pkl
        mv ~{ins_indel_scorer_30x_pkl} 30x.ins.indel.scorer.pkl
        mv ~{dup_indel_scorer_30x_pkl} 30x.dup.indel.scorer.pkl
        mv ~{insdup_indel_scorer_30x_pkl} 30x.insdup.indel.scorer.pkl
        mv ~{inv_indel_scorer_30x_pkl} 30x.inv.indel.scorer.pkl
        mv ~{bnd_indel_scorer_30x_pkl} 30x.bnd.indel.scorer.pkl
        mv ~{del_indel_calibrationScores_30x_hdf5} 30x.del.indel.calibrationScores.hdf5
        mv ~{ins_indel_calibrationScores_30x_hdf5} 30x.ins.indel.calibrationScores.hdf5
        mv ~{dup_indel_calibrationScores_30x_hdf5} 30x.dup.indel.calibrationScores.hdf5
        mv ~{insdup_indel_calibrationScores_30x_hdf5} 30x.insdup.indel.calibrationScores.hdf5
        mv ~{inv_indel_calibrationScores_30x_hdf5} 30x.inv.indel.calibrationScores.hdf5
        mv ~{bnd_indel_calibrationScores_30x_hdf5} 30x.bnd.indel.calibrationScores.hdf5 

        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE || [ -n "${LINE}" ]; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ "${TEST}" != "0" ]; then
                continue
            fi
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}
            
            # Filtering
            awk -F ',' -v sample="${SAMPLE_ID}" '$1 == sample { print $2 }' ~{sample_coverages_csv} > ${SAMPLE_ID}_coverage.txt
            N_ROWS=$(wc -l < ${SAMPLE_ID}_coverage.txt)
            if [ ${N_ROWS} -ne 1 ]; then
                echo "Expected exactly one row for sample ${SAMPLE_ID}, found ${N_ROWS}" 1>&2
                exit 1
            fi
            COVERAGE=$(head -n 1 ${SAMPLE_ID}_coverage.txt)
            if [ $(echo "${COVERAGE} < 22.5" | bc) -eq 1 ]; then
                PREFIX="15x"
            else
                PREFIX="30x"
            fi
            Score ${SAMPLE_ID} ${PREFIX} del ${RAM_PER_THREAD_MB}
            Score ${SAMPLE_ID} ${PREFIX} ins ${RAM_PER_THREAD_MB}
            Score ${SAMPLE_ID} ${PREFIX} dup ${RAM_PER_THREAD_MB}
            Score ${SAMPLE_ID} ${PREFIX} insdup ${RAM_PER_THREAD_MB}
            Score ${SAMPLE_ID} ${PREFIX} inv ${RAM_PER_THREAD_MB}
            Score ${SAMPLE_ID} ${PREFIX} bnd ${RAM_PER_THREAD_MB}

            # Assembling a single scored ultralong BCF (the 5 ultralong SVTYPEs).
            # Adding SVLEN to symbolic ALTs, to avoid overcollapse in `bcftools
            # merge` downstream. BND is kept as its own ${SAMPLE_ID}_bnd.bcf.
            bcftools concat --threads ${N_THREADS} --allow-overlaps --output-type v ${SAMPLE_ID}_del.bcf ${SAMPLE_ID}_ins.bcf ${SAMPLE_ID}_insdup.bcf ${SAMPLE_ID}_dup.bcf ${SAMPLE_ID}_inv.bcf --output ${SAMPLE_ID}_ultralong.vcf
            rm -f ${SAMPLE_ID}_del.bcf* ${SAMPLE_ID}_ins.bcf* ${SAMPLE_ID}_insdup.bcf* ${SAMPLE_ID}_dup.bcf* ${SAMPLE_ID}_inv.bcf*
            java AddSvlenToSymbolicAlt ${SAMPLE_ID}_ultralong.vcf > ${SAMPLE_ID}_ultralong_out.vcf
            rm -f ${SAMPLE_ID}_ultralong.vcf ; mv ${SAMPLE_ID}_ultralong_out.vcf ${SAMPLE_ID}_ultralong_in.vcf
            bcftools sort --max-mem ${RAM_PER_THREAD_MB}M --output-type b ${SAMPLE_ID}_ultralong_in.vcf --output ${SAMPLE_ID}_ultralong.bcf
            rm -f ${SAMPLE_ID}_ultralong_in.vcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_ultralong.bcf

            # Uploading. Output names (${SAMPLE_ID}_ultralong.bcf / ${SAMPLE_ID}_bnd.bcf)
            # match the per-sample inputs expected by Workflow D (PerSuffix).
            gsutil mv ./${SAMPLE_ID}_ultralong.'bcf*' ~{remote_outdir}/
            gsutil mv ./${SAMPLE_ID}_bnd.'bcf*' ~{remote_outdir}/
            gsutil mv ./${SAMPLE_ID}_'*_xgboost.csv' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gsutil mv ${SAMPLE_ID}.done ~{remote_outdir}/ && echo 0 || echo 1
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done 3< chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
    }
}
