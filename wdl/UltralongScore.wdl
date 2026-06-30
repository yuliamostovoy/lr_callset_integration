version 1.0


# For DEL, INV, DUP, INSDUP:
#
# annotations_custom = [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "BIN_BEFORE_COVERAGE","BIN_LEFT_COVERAGE","BIN_1_COVERAGE","BIN_2_COVERAGE","BIN_3_COVERAGE","BIN_4_COVERAGE","BIN_5_COVERAGE","BIN_6_COVERAGE","BIN_7_COVERAGE","BIN_8_COVERAGE","BIN_9_COVERAGE","BIN_10_COVERAGE","BIN_RIGHT_COVERAGE","BIN_AFTER_COVERAGE","BIN_LEFT_MAPQ","BIN_RIGHT_MAPQ","BIN_LEFT_SECONDARY","BIN_RIGHT_SECONDARY","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4"
#                      ]
# annotations_fex =    [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"
#                      ]
# annotations_cutefc = [ "SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND"
#                      ]
# annotations_all =    [ "SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "BIN_BEFORE_COVERAGE","BIN_LEFT_COVERAGE","BIN_1_COVERAGE","BIN_2_COVERAGE","BIN_3_COVERAGE","BIN_4_COVERAGE","BIN_5_COVERAGE","BIN_6_COVERAGE","BIN_7_COVERAGE","BIN_8_COVERAGE","BIN_9_COVERAGE","BIN_10_COVERAGE","BIN_RIGHT_COVERAGE","BIN_AFTER_COVERAGE","BIN_LEFT_MAPQ","BIN_RIGHT_MAPQ","BIN_LEFT_SECONDARY","BIN_RIGHT_SECONDARY","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4",
#                        "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK",
#                        "CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND" 
#                      ]
# annotations_all_except_genotyper = [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                                      "BIN_BEFORE_COVERAGE","BIN_LEFT_COVERAGE","BIN_1_COVERAGE","BIN_2_COVERAGE","BIN_3_COVERAGE","BIN_4_COVERAGE","BIN_5_COVERAGE","BIN_6_COVERAGE","BIN_7_COVERAGE","BIN_8_COVERAGE","BIN_9_COVERAGE","BIN_10_COVERAGE","BIN_RIGHT_COVERAGE","BIN_AFTER_COVERAGE","BIN_LEFT_MAPQ","BIN_RIGHT_MAPQ","BIN_LEFT_SECONDARY","BIN_RIGHT_SECONDARY","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4",
#                                      "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK" 
#                                    ]
# For INS:
#
# annotations_custom = [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "BIN_POS","BIN_POINT_MAPQ","BIN_POINT_SECONDARY","PL","PR","PL_PL_1","PL_PL_2","PL_PL_3","PL_PL_4","PL_PR_1","PL_PR_2","PL_PR_3","PL_PR_4","PR_PR_1","PR_PR_2","PR_PR_3","PR_PR_4"
#                      ]
# annotations_fex =    [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"
#                      ]
# annotations_cutefc = [ "SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND"
#                      ]
# annotations_all =    [ "SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "BIN_POS","BIN_POINT_MAPQ","BIN_POINT_SECONDARY","PL","PR","PL_PL_1","PL_PL_2","PL_PL_3","PL_PL_4","PL_PR_1","PL_PR_2","PL_PR_3","PL_PR_4","PR_PR_1","PR_PR_2","PR_PR_3","PR_PR_4",
#                        "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK",
#                        "CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND"
#                      ]
# annotations_all_except_genotyper = [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                                      "BIN_POS","BIN_POINT_MAPQ","BIN_POINT_SECONDARY","PL","PR","PL_PL_1","PL_PL_2","PL_PL_3","PL_PL_4","PL_PR_1","PL_PR_2","PL_PR_3","PL_PR_4","PR_PR_1","PR_PR_2","PR_PR_3","PR_PR_4",
#                                      "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"
#                                    ]
#
# For BND:
#
# annotations_custom = [ "GT_COUNT","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "BIN_POS_0","BIN_POS_1","BIN_POINT_MAPQ_0","BIN_POINT_MAPQ_1","BIN_POINT_SECONDARY_0","BIN_POINT_SECONDARY_1","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4"
#                      ]
# annotations_fex =    [ "GT_COUNT","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "FEX_DEPTH_RATIO","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"
#                      ]
# annotations_cutefc = [ "SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND"
#                      ]
# annotations_all =    [ "SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                        "BIN_POS_0","BIN_POS_1","BIN_POINT_MAPQ_0","BIN_POINT_MAPQ_1","BIN_POINT_SECONDARY_0","BIN_POINT_SECONDARY_1","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4",
#                        "FEX_DEPTH_RATIO","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK",
#                        "CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND" 
#                      ]
# annotations_all_except_genotyper = [ "GT_COUNT","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
#                                      "BIN_POS_0","BIN_POS_1","BIN_POINT_MAPQ_0","BIN_POINT_MAPQ_1","BIN_POINT_SECONDARY_0","BIN_POINT_SECONDARY_1","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4",
#                                      "FEX_DEPTH_RATIO","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK" 
#                                    ]
#
workflow SV_Integration_UltralongScore {
    input {
        String svtype

        File input_vcf_gz
        File input_vcf_gz_tbi
        File resource_vcf_gz
        File resource_vcf_gz_tbi
        String remote_outdir
                
        File? training_resource_bed
        String exclude_chromosomes_string = " "
        File reference_fa
        File reference_fai

        Array[String]? annotations_custom
        Array[String]? annotations_fex
        Array[String]? annotations_cutefc
        Array[String] annotations_all
        Array[String]? annotations_all_except_genotyper
        Int annotations_have_gt_count

        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        
        String docker_image = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
    }
    parameter_meta {
        input_vcf_gz: "Assumed to contain all the annotations used in this workflow."
        remote_outdir: "Without final slash"
        exclude_chromosomes_string: "Example: -XL chr1 -XL chr2 -XL chr3 -XL chr4 -XL chr5"
        hyperparameters_json: "Parameters for `gatk TrainVariantAnnotationsModel`."
    }
    
    if (defined(annotations_custom)) {
        call Score as score_custom {
            input:
                svtype = svtype,
                id = svtype + "_custom",
                annotations = select_first([annotations_custom]),
                annotations_have_gt_count = annotations_have_gt_count,
                input_vcf_gz = input_vcf_gz,
                input_vcf_gz_tbi = input_vcf_gz_tbi,
                resource_vcf_gz = resource_vcf_gz,
                resource_vcf_gz_tbi = resource_vcf_gz_tbi,
                remote_outdir = remote_outdir,
                training_resource_bed = training_resource_bed,
                exclude_chromosomes_string = exclude_chromosomes_string,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                training_python_script = training_python_script,
                scoring_python_script = scoring_python_script,
                hyperparameters_json = hyperparameters_json,
                docker_image = docker_image
        }
    }
    if (defined(annotations_fex)) {
        call Score as score_fex {
            input:
                svtype = svtype,
                id = svtype + "_fex",
                annotations = select_first([annotations_fex]),
                annotations_have_gt_count = annotations_have_gt_count,
                input_vcf_gz = input_vcf_gz,
                input_vcf_gz_tbi = input_vcf_gz_tbi,
                resource_vcf_gz = resource_vcf_gz,
                resource_vcf_gz_tbi = resource_vcf_gz_tbi,
                remote_outdir = remote_outdir,
                training_resource_bed = training_resource_bed,
                exclude_chromosomes_string = exclude_chromosomes_string,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                training_python_script = training_python_script,
                scoring_python_script = scoring_python_script,
                hyperparameters_json = hyperparameters_json,
                docker_image = docker_image
        }
    }
    if (defined(annotations_cutefc)) {
        call Score as score_cutefc {
            input:
                svtype = svtype,
                id = svtype + "_cutefc",
                annotations = select_first([annotations_cutefc]),
                annotations_have_gt_count = annotations_have_gt_count,
                input_vcf_gz = input_vcf_gz,
                input_vcf_gz_tbi = input_vcf_gz_tbi,
                resource_vcf_gz = resource_vcf_gz,
                resource_vcf_gz_tbi = resource_vcf_gz_tbi,
                remote_outdir = remote_outdir,
                training_resource_bed = training_resource_bed,
                exclude_chromosomes_string = exclude_chromosomes_string,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                training_python_script = training_python_script,
                scoring_python_script = scoring_python_script,
                hyperparameters_json = hyperparameters_json,
                docker_image = docker_image
        }
    }
    call Score as score_all {
        input:
            svtype = svtype,
            id = svtype + "_all",
            annotations = annotations_all,
            annotations_have_gt_count = annotations_have_gt_count,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            exclude_chromosomes_string = exclude_chromosomes_string,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    if (defined(annotations_all_except_genotyper)) {
        call Score as score_all_except_genotyper {
            input:
                svtype = svtype,
                id = svtype + "_all_except_genotyper",
                annotations = select_first([annotations_all_except_genotyper]),
                annotations_have_gt_count = annotations_have_gt_count,
                input_vcf_gz = input_vcf_gz,
                input_vcf_gz_tbi = input_vcf_gz_tbi,
                resource_vcf_gz = resource_vcf_gz,
                resource_vcf_gz_tbi = resource_vcf_gz_tbi,
                remote_outdir = remote_outdir,
                training_resource_bed = training_resource_bed,
                exclude_chromosomes_string = exclude_chromosomes_string,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                training_python_script = training_python_script,
                scoring_python_script = scoring_python_script,
                hyperparameters_json = hyperparameters_json,
                docker_image = docker_image
        }
    }
    
    output {
    }
}


# Performance on a 2-core, 16GB VM, all HPRC+HGSVC samples:
#
# TOOL                                  CPU            RAM              TIME
# ExtractVariantAnnotations             ???             8G               20s
# TrainVariantAnnotationsModel          ???           200M               10s
# ScoreVariantAnnotations               ???           800M               20s
#
task Score {
    input {
        String svtype
        String id
        
        File input_vcf_gz
        File input_vcf_gz_tbi
        File resource_vcf_gz
        File resource_vcf_gz_tbi
        String remote_outdir
        
        File? training_resource_bed
        String exclude_chromosomes_string
        File reference_fa
        File reference_fai

        Array[String] annotations
        Int annotations_have_gt_count
        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 10
        Int disk_size_gb = 10
    }
    parameter_meta {
    }
    
    String docker_dir = "/root"
    
    command <<<
        set -euxo pipefail
        
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        export GATK_LOCAL_JAR="/root/gatk.jar"
        



        # ----------------------- Steps of the pipeline ------------------------

        # Adds field `INFO/GT_COUNT` to the input VCF, which is overwritten.
        #
        function AddGtCount() {
            local INPUT_VCF_GZ=$1

            bcftools query --format '%CHROM\t%POS\t%ID\t[%GT]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($4=="0/0" || $4=="0|0" || $4=="./."  || $4==".|." || $4=="./0" || $4==".|0" || $4=="0/." || $4=="0|." || $4=="0" || $4==".") GT_COUNT=0; \
                else if ($4=="0/1" || $4=="0|1" || $4=="1/0" || $4=="1|0" || $4=="./1" || $4==".|1" || $4=="1/." || $4=="1|." || $4=="1") GT_COUNT=1; \
                else if ($4=="1/1" || $4=="1|1") GT_COUNT=2; \
                printf("%s\t%d\t%s\t%d\n",$1,$2,$3,GT_COUNT); \
            }' | bgzip -c > annotations.tsv.gz
            tabix -f -s1 -b2 -e2 annotations.tsv.gz
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="Original GT converted to an integer in {0,1,2}.">' > header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/GT_COUNT'
            bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns ${COLUMNS} --output-type z ${INPUT_VCF_GZ} --output annotated.vcf.gz
            rm -f annotations.tsv.gz header.txt ${INPUT_VCF_GZ}
            mv annotated.vcf.gz ${INPUT_VCF_GZ}
            bcftools index -f -t ${INPUT_VCF_GZ}
        }


        # Copies score fields from INFO to FORMAT so that they are retained as
        # sample fields by downstream merges. This does not filter records.
        #
        function CopyScoresToFormat() {
            local INPUT_VCF_GZ=$1

            echo '##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="Score according to the XGBoost model">' > score_header.txt
            echo '##FORMAT=<ID=CALIBRATION_SENSITIVITY,Number=1,Type=Float,Description="Calibration sensitivity according to the model applied by ScoreVariantAnnotations">' >> score_header.txt
            bcftools query --format '%CHROM\t%POS\t%ID\t%SCORE\t%CALIBRATION_SENSITIVITY\n' ${INPUT_VCF_GZ} | bgzip -c > score_format.tsv.gz
            tabix -f -s1 -b2 -e2 score_format.tsv.gz
            bcftools annotate --header-lines score_header.txt --annotations score_format.tsv.gz --columns CHROM,POS,~ID,FORMAT/SCORE,FORMAT/CALIBRATION_SENSITIVITY --output-type z ${INPUT_VCF_GZ} --output score_format.vcf.gz
            rm -f score_header.txt score_format.tsv.gz score_format.tsv.gz.tbi ${INPUT_VCF_GZ} ${INPUT_VCF_GZ}.tbi
            mv score_format.vcf.gz ${INPUT_VCF_GZ}
            bcftools index -f -t ${INPUT_VCF_GZ}
        }




        # ---------------------------- Main program ----------------------------

        # 1. Ensuring that the input VCFs have the correct format
        if [ ~{svtype} != "bnd" ]; then
            bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ~{input_vcf_gz} --output input_cleaned.vcf.gz
            bcftools index -f -t input_cleaned.vcf.gz
            rm -f ~{input_vcf_gz}
            bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ~{resource_vcf_gz} --output resource_cleaned.vcf.gz
            bcftools index -f -t resource_cleaned.vcf.gz
            rm -f ~{resource_vcf_gz}
        else
            # We do not run the command above, since it seems to destroy 
            # BND's ALTs (example: N]chr5:181473415] -> GNcNNNNNNNNNNNNNN ):
            mv ~{input_vcf_gz} input_cleaned.vcf.gz
            mv ~{input_vcf_gz_tbi} input_cleaned.vcf.gz.tbi
            mv ~{resource_vcf_gz} resource_cleaned.vcf.gz
            mv ~{resource_vcf_gz_tbi} resource_cleaned.vcf.gz.tbi
        fi
        N_RECORDS_INPUT="$(bcftools index --nrecords input_cleaned.vcf.gz)"
        if [ ~{annotations_have_gt_count} -eq 1 ]; then
            AddGtCount input_cleaned.vcf.gz
        fi
        N_RECORDS_RESOURCE="$(bcftools index --nrecords resource_cleaned.vcf.gz)"
        echo "Total records: ${N_RECORDS_INPUT}  Marked as true: ${N_RECORDS_RESOURCE}" 1>&2

        # 2. Scoring
        if ~{defined(training_resource_bed)}
        then
            BED_FLAG="-L ~{training_resource_bed}"
        else
            BED_FLAG=""
        fi
        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ExtractVariantAnnotations -V input_cleaned.vcf.gz ~{exclude_chromosomes_string} -O extract -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true resource_cleaned.vcf.gz --maximum-number-of-unlabeled-variants 1000000000 --mode INDEL --mnp-type INDEL ${BED_FLAG}
        ls -laht 1>&2
        # Output:
        # extract.annot.hdf5
        # extract.unlabeled.annot.hdf5
        # extract.vcf.gz
        # extract.vcf.gz.tbi
        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" TrainVariantAnnotationsModel --annotations-hdf5 extract.annot.hdf5 --unlabeled-annotations-hdf5 extract.unlabeled.annot.hdf5 --model-backend PYTHON_SCRIPT --python-script ~{training_python_script} --hyperparameters-json ~{hyperparameters_json} -O train.train --mode INDEL --verbosity DEBUG
        ls -laht 1>&2
        # Output: 
        # train.train.indel.unlabeledScores.hdf5
        # train.train.indel.calibrationScores.hdf5
        # train.train.indel.trainingScores.hdf5
        # train.train.indel.scorer.pkl
        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ScoreVariantAnnotations -V input_cleaned.vcf.gz -O score -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true resource_cleaned.vcf.gz --resource:extracted,extracted=true extract.vcf.gz --model-prefix train.train --model-backend PYTHON_SCRIPT --python-script ~{scoring_python_script} --mode INDEL --mnp-type INDEL --ignore-all-filters --verbosity DEBUG
        ls -laht 1>&2
        # Output:
        # score.vcf.gz
        # score.vcf.gz.tbi
        # score.annot.hdf5
        # score.scores.hdf5
        CopyScoresToFormat score.vcf.gz
        gsutil -m mv score.vcf.gz ~{remote_outdir}/~{id}_score.vcf.gz
        gsutil -m mv score.vcf.gz.tbi ~{remote_outdir}/~{id}_score.vcf.gz.tbi
    >>>
    
    output {
        File extract_annot_hdf5 = "extract.annot.hdf5"
        File extract_unlabeled_annot_hdf5 = "extract.unlabeled.annot.hdf5"

        File train_indel_unlabeled_scores_hdf5 = "train.train.indel.unlabeledScores.hdf5"
        File train_indel_calibration_scores_hdf5 = "train.train.indel.calibrationScores.hdf5"
        File train_indel_training_scores_hdf5 = "train.train.indel.trainingScores.hdf5"
        File train_indel_scorer_pkl = "train.train.indel.scorer.pkl"

        File score_annot_hdf5 = "score.annot.hdf5"
        File score_scores_hdf5 = "score.scores.hdf5"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
