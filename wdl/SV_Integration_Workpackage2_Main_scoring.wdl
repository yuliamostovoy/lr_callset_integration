version 1.0


# Given a single-sample VCF, the program scores it with XGBoost, retains all
# variants, and then splits it into ~100 pieces, in order to run bcftools merge
# over all samples on parallel chunks.
#
workflow SV_Integration_Workpackage2 {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        String filter_string = "none"
        
        String remote_indir
        String remote_outdir
        
        Array[String] annotations = ["KS_1","KS_2","SQ","GQ","DP","AD_NON_ALT","AD_ALL","GT_COUNT","SUPP_PAV","SUPP_SNIFFLES","SUPP_PBSV","SVLEN"]
        
        File training_resource_bed
        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        
        String docker_image = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
    }
    parameter_meta {
        split_for_bcftools_merge_csv: "A partition that covers all chromosomes. Every line is a 0-based, half-open, consecutive chunk of a chromosome. Lines are assumed to be sorted."
        filter_string: "Deprecated; variants are scored but no score filter is applied."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        training_resource_bed: "The same BED used in `SV_Integration_Workpackage1`."
        hyperparameters_json: "Parameters for `gatk TrainVariantAnnotationsModel`."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            filter_string = filter_string,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            annotations = annotations,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
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
# ExtractVariantAnnotations           900 MB
# TrainVariantAnnotationsModel        200 MB
# ScoreVariantAnnotations               1 GB
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        String filter_string
        
        String remote_indir
        String remote_outdir
        
        File training_resource_bed

        Array[String] annotations
        File training_python_script
        File scoring_python_script
        File hyperparameters_json

        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 3
        Int disk_size_gb = 20
        Int preemptible_number = 4
        String upstream_signal = ""
    }
    parameter_meta {
        upstream_signal: "Ordering-only handshake for orchestrator workflows: set to the `done` output of the upstream per-sample step so this task starts only after it. Ignored by standalone runs."
    }
    
    String docker_dir = "/root"
    
    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        GSUTIL_DELAY_S="600"
        export GATK_LOCAL_JAR="/root/gatk.jar"
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local REMOTE_DIR=$2
            
            gsutil cp ${REMOTE_DIR}/${SAMPLE_ID}_kanpig.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_training.vcf.'gz*' .
        }
        
        
        # Deletes all files and directories related to the sample.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -rf ./${SAMPLE_ID}_*
        }
        
        
        # Remark: the procedure's input and output are indexed `.vcf.gz`.
        #
        function JointVcfFiltering() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local RESOURCE_VCF_GZ=$3

            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ExtractVariantAnnotations -V ${INPUT_VCF_GZ} -O ${SAMPLE_ID}_extract -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true ${RESOURCE_VCF_GZ} --maximum-number-of-unlabeled-variants 10000000 --mode INDEL --mnp-type INDEL -L ~{training_resource_bed}
            ls -laht
            # Output:
            # ${SAMPLE_ID}_extract.annot.hdf5
            # ${SAMPLE_ID}_extract.unlabeled.annot.hdf5
            # ${SAMPLE_ID}_extract.vcf.gz
            # ${SAMPLE_ID}_extract.vcf.gz.tbi
            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" TrainVariantAnnotationsModel --annotations-hdf5 ${SAMPLE_ID}_extract.annot.hdf5 --unlabeled-annotations-hdf5 ${SAMPLE_ID}_extract.unlabeled.annot.hdf5 --model-backend PYTHON_SCRIPT --python-script ~{training_python_script} --hyperparameters-json ~{hyperparameters_json} -O ${SAMPLE_ID}.train --mode INDEL --verbosity DEBUG
            ls -laht
            # Output: 
            # ${SAMPLE_ID}.train.*
            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ScoreVariantAnnotations -V ${INPUT_VCF_GZ} -O ${SAMPLE_ID}_score -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true ${RESOURCE_VCF_GZ} --resource:extracted,extracted=true ${SAMPLE_ID}_extract.vcf.gz --model-prefix ${SAMPLE_ID}.train --model-backend PYTHON_SCRIPT --python-script ~{scoring_python_script} --mode INDEL --mnp-type INDEL --ignore-all-filters --verbosity DEBUG
            ls -laht
            # Output:
            # ${SAMPLE_ID}_score.vcf.gz
            # ${SAMPLE_ID}_score.vcf.gz.tbi
            # ${SAMPLE_ID}_score.annot.hdf5
            # ${SAMPLE_ID}_score.scores.hdf5
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_extract.annot.hdf5 ${SAMPLE_ID}_extract.unlabeled.annot.hdf5 ${SAMPLE_ID}_extract.vcf.gz* ${SAMPLE_ID}.train.* ${SAMPLE_ID}_score.annot.hdf5 ${SAMPLE_ID}_score.scores.hdf5
        }
        
        
        # Copies the following fields from INFO to FORMAT, so that they are
        # preserved by the inter-sample merge downstream:
        #
        # SUPP_*, SCORE, CALIBRATION_SENSITIVITY
        #
        # Remark: the procedure outputs an indexed `.bcf`.
        #
        # @param 2 A VCF where all IDs are distinct. This is guaranteed by
        # workpackages upstream.
        #
        function CopyInfoToFormat() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            echo '##FORMAT=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="Score according to the XGBoost model">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=CALIBRATION_SENSITIVITY,Number=1,Type=Float,Description="Calibration sensitivity according to the model applied by ScoreVariantAnnotations">' >> ${SAMPLE_ID}_header.txt
            # REF,ALT are emitted (cols 4,5) and added to --columns so the `~ID`
            # match actually engages. Without REF,ALT present bcftools ignores
            # `~ID` and matches on CHROM,POS only, mis-assigning FORMAT values
            # across records sharing a start coordinate. IDs are distinct (see
            # @param note above).
            bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SUPP_PBSV\t%SUPP_SNIFFLES\t%SUPP_PAV\t%SCORE\t%CALIBRATION_SENSITIVITY\n' ${INPUT_VCF_GZ} | bgzip -c > ${SAMPLE_ID}_format.tsv.gz
            tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_format.tsv.gz
            bcftools annotate --threads ${N_THREADS} --header-lines ${SAMPLE_ID}_header.txt --annotations ${SAMPLE_ID}_format.tsv.gz --columns CHROM,POS,~ID,REF,ALT,FORMAT/SUPP_PBSV,FORMAT/SUPP_SNIFFLES,FORMAT/SUPP_PAV,FORMAT/SCORE,FORMAT/CALIBRATION_SENSITIVITY --output-type b ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_scored.bcf
            bcftools index --threads ${N_THREADS} ${SAMPLE_ID}_scored.bcf
            (bcftools view --no-header ${SAMPLE_ID}_scored.bcf | head -n 1 || echo "0") 1>&2
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_header.txt ${SAMPLE_ID}_format.tsv.gz*
        }
        
        
        # Assumes that `CopyInfoToFormat()` has already been executed.
        #
        function PrintDebugInformation() {
            local SAMPLE_ID=$1
            local INPUT_BCF=$2
            
            rm -rf ${SAMPLE_ID}_xgboost.csv
            local N_RECORDS_BEFORE_FILTERING=$(bcftools index --nrecords ${INPUT_BCF})
            for THRESHOLD in 0.7 0.8 0.9 0.95 ; do
                local N_RECORDS_AFTER_FILTERING=$( bcftools query --format '%ID' --include "FORMAT/CALIBRATION_SENSITIVITY<=${THRESHOLD}" ${INPUT_BCF} | wc -l )
                local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_FILTERING} / ${N_RECORDS_BEFORE_FILTERING}" | bc )
                echo "${N_RECORDS_AFTER_FILTERING},${N_RECORDS_BEFORE_FILTERING},${PERCENT},Number of records with CALIBRATION_SENSITIVITY<=${THRESHOLD}" >> ${SAMPLE_ID}_xgboost.csv
            done
            if [ "~{filter_string}" != "none" ]; then
                local N_RECORDS_AFTER_FILTERING=$( bcftools query --format '%ID' --include "~{filter_string}" ${INPUT_BCF} | wc -l )
                local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_FILTERING} / ${N_RECORDS_BEFORE_FILTERING}" | bc )
                echo "${N_RECORDS_AFTER_FILTERING},${N_RECORDS_BEFORE_FILTERING},${PERCENT},Number of records that pass the specified filter" >> ${SAMPLE_ID}_xgboost.csv
            fi
        }
        
        
        # Remark: the procedure's input and output are indexed `.bcf`.
        # 
        function FilterChunkUpload() {
            local SAMPLE_ID=$1
            local INPUT_BCF=$2
            
            i="0"
            local INTERVAL
            while read -u 4 INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                # Remark: we use `targets` rather than `regions` because
                # the former considers just the POS coordinate for overlaps.
                bcftools view --threads ${N_THREADS} --targets-file ${SAMPLE_ID}.bed --output-type b ${INPUT_BCF} --output ${SAMPLE_ID}_chunk_${i}.bcf
                bcftools index --threads ${N_THREADS} ${SAMPLE_ID}_chunk_${i}.bcf
                gsutil mv ${SAMPLE_ID}_chunk_${i}.bcf ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf
                gsutil mv ${SAMPLE_ID}_chunk_${i}.bcf.csi ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf.csi
                i=$(( ${i} + 1 ))
            done 4< ~{split_for_bcftools_merge_csv}
            touch ${SAMPLE_ID}.done
            gsutil mv ${SAMPLE_ID}.done ~{remote_outdir}/ && echo 0 || echo 1
        }

        

        
        # ---------------------------- Main program ----------------------------
        
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        N_OUTPUT_CHUNKS=$(wc -l < ~{split_for_bcftools_merge_csv})
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Filtering
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}
            JointVcfFiltering ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz ${SAMPLE_ID}_training.vcf.gz
            CopyInfoToFormat ${SAMPLE_ID} ${SAMPLE_ID}_score.vcf.gz
            PrintDebugInformation ${SAMPLE_ID} ${SAMPLE_ID}_scored.bcf
            FilterChunkUpload ${SAMPLE_ID} ${SAMPLE_ID}_scored.bcf
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done 3< chunk.csv

        # Batch-completion signal for orchestrator ordering. Ignored standalone.
        echo "done" > wp2.signal
    >>>

    output {
        String done = read_string("wp2.signal")
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
