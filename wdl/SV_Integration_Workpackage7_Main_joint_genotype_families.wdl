version 1.0


# Builds family-specific candidate VCFs from the WP8 truvari-collapsed cohort
# VCF, re-genotypes each family member with kanpig, and splits each sample's
# re-genotyped VCF into chunks for a downstream bcftools merge by ID.
#
workflow SV_Integration_Workpackage9_families {
    input {
        Array[String] family_ids
        File ped
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bais
        Array[String] aligned_bams
        File split_for_bcftools_merge_csv
        
        String remote_indir
        String remote_outdir
        String requester_pays_project = ""
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed
        
        String kanpig_params_cohort = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        family_ids: "Family IDs to process. Each ID must match column 1 of ped."
        ped: "Standard 6-column PED: family_id, sample_id, paternal_id, maternal_id, sex, phenotype."
        sample_ids: "Sample IDs. Must be in the same order as sample_sexes, aligned_bais, and aligned_bams."
        sample_sexes: "Sample sexes. Must be in the same order as sample_ids, aligned_bais, and aligned_bams."
        aligned_bais: "Remote aligned BAM index paths. Must be in the same order as sample_ids, sample_sexes, and aligned_bams."
        aligned_bams: "Remote aligned BAM paths. Must be in the same order as sample_ids, sample_sexes, and aligned_bais."
        split_for_bcftools_merge_csv: "A partition that covers all chromosomes. Every line is a 0-based, half-open, consecutive chunk of a chromosome. Lines are assumed to be sorted."
        remote_indir: "Without final slash. Contains WP8 genome-wide truvari_collapsed.bcf and truvari_collapsed.bcf.csi."
        remote_outdir: "Without final slash. Output sample BCF chunks are written under chunk_N/sample.bcf."
        requester_pays_project: "Google Cloud project to bill for requester-pays BAM/BAI buckets. Leave blank for non-requester-pays buckets."
    }
    
    call Impl {
        input:
            family_ids = family_ids,
            ped = ped,
            sample_ids = sample_ids,
            sample_sexes = sample_sexes,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            requester_pays_project = requester_pays_project,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            autosomes_bed = autosomes_bed,
            kanpig_params_cohort = kanpig_params_cohort,
            docker_image = docker_image
    }
    
    output {
    }
}


task Impl {
    input {
        Array[String] family_ids
        File ped
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bais
        Array[String] aligned_bams
        File split_for_bcftools_merge_csv
        
        String remote_indir
        String remote_outdir
        String requester_pays_project
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed
        
        String kanpig_params_cohort
        String docker_image
        
        Int n_cpu = 6
        Int ram_size_gb = 8
        Int disk_size_gb = 100
        Int preemptible_number = 4
    }
    parameter_meta {
        disk_size_gb: "Increase for large family batches or large BAMs."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        export RUST_BACKTRACE="full"
        GCLOUD_STORAGE_BILLING_FLAGS=""
        if [ -n "~{requester_pays_project}" ]; then
            GCLOUD_STORAGE_BILLING_FLAGS="--billing-project=~{requester_pays_project}"
        fi
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            local ALIGNED_BAI=$(echo "${LINE}" | cut -f 3)
            local ALIGNED_BAM=$(echo "${LINE}" | cut -f 4)
            
            local AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
            AVAILABLE_GB=${AVAILABLE_GB%G}
            AVAILABLE_GB=${AVAILABLE_GB%.*}
            local BAM_BYTES=$(gcloud storage ls -l ${GCLOUD_STORAGE_BILLING_FLAGS} "${ALIGNED_BAM}" | awk '$1 ~ /^[0-9]+$/ { print $1; exit }')
            if [ -z "${BAM_BYTES}" ]; then
                echo "ERROR: could not determine BAM size for ${ALIGNED_BAM}."
                exit 1
            fi
            local BAM_GB=${BAM_BYTES}
            BAM_GB=$(( (${BAM_GB} + 1073741823) / 1073741824 + 5 ))
            if [ ${BAM_GB} -gt ${AVAILABLE_GB} ]; then
                echo "ERROR: the BAM is larger than the available disk space. BAM size + slack: ${BAM_GB}GB. Available disk: ${AVAILABLE_GB}GB."
                exit 1
            fi
            
            date 1>&2
            gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} "${ALIGNED_BAM}" ./${SAMPLE_ID}_aligned.bam
            gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} "${ALIGNED_BAI}" ./${SAMPLE_ID}_aligned.bam.bai
            date 1>&2
            touch ${SAMPLE_ID}_aligned.bam.bai
        }
        
        
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_*
        }
        
        
        function BuildFamilyCandidateVcf() {
            local FAMILY_ID=$1
            
            awk -v family="${FAMILY_ID}" 'BEGIN { FS="[ \t]+" } $1==family && $2!="0" && $2!="." { print $2 }' ~{ped} | sort -u > ${FAMILY_ID}.samples.txt
            local N_FAMILY_SAMPLES=$(wc -l < ${FAMILY_ID}.samples.txt)
            if [ ${N_FAMILY_SAMPLES} -eq 0 ]; then
                echo "ERROR: family ${FAMILY_ID} has no samples in the PED."
                exit 1
            fi
            
            local SAMPLE_ID
            while read -u 4 SAMPLE_ID; do
                if ! awk -v sample="${SAMPLE_ID}" 'BEGIN { FS="\t"; found=0 } $1==sample { found=1 } END { exit(found ? 0 : 1) }' sample_metadata.tsv; then
                    echo "ERROR: sample ${SAMPLE_ID} from family ${FAMILY_ID} is missing from sample_ids."
                    exit 1
                fi
            done 4< ${FAMILY_ID}.samples.txt
            
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file ${FAMILY_ID}.samples.txt --output-type z cohort.bcf --output ${FAMILY_ID}_all.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${FAMILY_ID}_all.vcf.gz
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z ${FAMILY_ID}_all.vcf.gz --output ${FAMILY_ID}_present.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${FAMILY_ID}_present.vcf.gz
            rm -f ${FAMILY_ID}_all.vcf.gz*
            
            local N_RECORDS=$(bcftools index --nrecords ${FAMILY_ID}_present.vcf.gz)
            echo "${FAMILY_ID},${N_FAMILY_SAMPLES},${N_RECORDS},Number of family samples and family-present records" > ${FAMILY_ID}_family.csv
        }
        
        
        function Kanpig() {
            local FAMILY_ID=$1
            local SAMPLE_ID=$2
            local SEX=$3
            
            local PLOIDY_BED
            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi
            
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${SAMPLE_ID} --output-type z ${FAMILY_ID}_present.vcf.gz --output ${SAMPLE_ID}_personalized.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_personalized.vcf.gz
            
            # Remark: kanpig needs --sizemin >= --kmer.
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --sizemin 10 --sizemax ${INFINITY} ~{kanpig_params_cohort} --reference ~{reference_fa} --ploidy-bed ${PLOIDY_BED} --input ${SAMPLE_ID}_personalized.vcf.gz --reads ${SAMPLE_ID}_aligned.bam --out ${SAMPLE_ID}_out.vcf --sample ${SAMPLE_ID}
            rm -f ${SAMPLE_ID}_personalized.vcf.gz*
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_out.vcf --output ${SAMPLE_ID}_kanpig.vcf.gz
            rm -f ${SAMPLE_ID}_out.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_kanpig.vcf.gz
            
            local N_RECORDS=$(bcftools index --nrecords ${SAMPLE_ID}_kanpig.vcf.gz)
            local N_PRESENT_RECORDS=$(bcftools query --format '%ID\n' --include 'GT="alt"' ${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            local PERCENT="0"
            if [ ${N_RECORDS} -gt 0 ]; then
                PERCENT=$(echo "scale=2; 100 * ${N_PRESENT_RECORDS} / ${N_RECORDS}" | bc)
            fi
            echo "${N_PRESENT_RECORDS},${N_RECORDS},${PERCENT},Number of records that are marked as ALT by kanpig" >> ${SAMPLE_ID}_kanpig.csv
            
            local N_HETS_IN_AUTOSOMES=$(bcftools query --format '%ID\n' --include 'GT="het"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            local N_PRESENT_RECORDS_IN_AUTOSOMES=$(bcftools query --format '%ID\n' --include 'GT="alt"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            PERCENT="0"
            if [ ${N_PRESENT_RECORDS_IN_AUTOSOMES} -gt 0 ]; then
                PERCENT=$(echo "scale=2; 100 * ${N_HETS_IN_AUTOSOMES} / ${N_PRESENT_RECORDS_IN_AUTOSOMES}" | bc)
            fi
            echo "${N_HETS_IN_AUTOSOMES},${N_PRESENT_RECORDS_IN_AUTOSOMES},${PERCENT},Number of records in autosomes that are marked as HET by kanpig" >> ${SAMPLE_ID}_kanpig.csv
            ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigWindows ${SAMPLE_ID}_kanpig.vcf.gz | bgzip > ${SAMPLE_ID}_kanpig.bed.gz
        }
        
        
        function ChunkAndUpload() {
            local SAMPLE_ID=$1
            
            local i="0"
            local INTERVAL
            while read -u 5 INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type b ${SAMPLE_ID}_kanpig.vcf.gz --output ${SAMPLE_ID}_chunk_${i}.bcf
                bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_chunk_${i}.bcf
                gcloud storage cp ${SAMPLE_ID}_chunk_${i}.bcf ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf
                gcloud storage cp ${SAMPLE_ID}_chunk_${i}.bcf.csi ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf.csi
                i=$(( ${i} + 1 ))
            done 5< ~{split_for_bcftools_merge_csv}
            gcloud storage cp ${SAMPLE_ID}_kanpig.bed.gz ${SAMPLE_ID}_kanpig.csv ~{remote_outdir}/
        }
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        ~{docker_dir}/kanpig --version 1>&2
        
        cat > sample_ids.txt <<'EOF_SAMPLE_IDS'
~{sep="\n" sample_ids}
EOF_SAMPLE_IDS
        cat > sample_sexes.txt <<'EOF_SAMPLE_SEXES'
~{sep="\n" sample_sexes}
EOF_SAMPLE_SEXES
        cat > aligned_bais.txt <<'EOF_ALIGNED_BAIS'
~{sep="\n" aligned_bais}
EOF_ALIGNED_BAIS
        cat > aligned_bams.txt <<'EOF_ALIGNED_BAMS'
~{sep="\n" aligned_bams}
EOF_ALIGNED_BAMS
        cat > family_ids.txt <<'EOF_FAMILY_IDS'
~{sep="\n" family_ids}
EOF_FAMILY_IDS
        grep -v '^[[:space:]]*$' family_ids.txt | sort -u > family_ids.unique.txt
        mv family_ids.unique.txt family_ids.txt
        
        N_SAMPLE_IDS=$(wc -l < sample_ids.txt)
        N_SAMPLE_SEXES=$(wc -l < sample_sexes.txt)
        N_ALIGNED_BAIS=$(wc -l < aligned_bais.txt)
        N_ALIGNED_BAMS=$(wc -l < aligned_bams.txt)
        if [ ${N_SAMPLE_IDS} -ne ${N_SAMPLE_SEXES} ] || [ ${N_SAMPLE_IDS} -ne ${N_ALIGNED_BAIS} ] || [ ${N_SAMPLE_IDS} -ne ${N_ALIGNED_BAMS} ]; then
            echo "ERROR: sample_ids, sample_sexes, aligned_bais, and aligned_bams must have the same length."
            echo "sample_ids=${N_SAMPLE_IDS}, sample_sexes=${N_SAMPLE_SEXES}, aligned_bais=${N_ALIGNED_BAIS}, aligned_bams=${N_ALIGNED_BAMS}"
            exit 1
        fi
        paste sample_ids.txt sample_sexes.txt aligned_bais.txt aligned_bams.txt > sample_metadata.tsv
        
        # Localizing the WP8 cohort VCF.
        ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/truvari_collapsed.'bcf*' .
        mv truvari_collapsed.bcf cohort.bcf
        mv truvari_collapsed.bcf.csi cohort.bcf.csi
        
        while read -u 3 FAMILY_ID; do
            if [ -z "${FAMILY_ID}" ] || [[ "${FAMILY_ID}" == \#* ]]; then
                continue
            fi
            
            BuildFamilyCandidateVcf ${FAMILY_ID}
            
            while read -u 4 SAMPLE_ID; do
                TEST=$(gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0")
                if [ ${TEST} != "0" ]; then
                    continue
                fi
                
                LINE=$(awk -v sample="${SAMPLE_ID}" 'BEGIN { FS="\t" } $1==sample { print; exit }' sample_metadata.tsv)
                SEX=$(echo "${LINE}" | cut -f 2)
                LocalizeSample ${SAMPLE_ID} "${LINE}"
                Kanpig ${FAMILY_ID} ${SAMPLE_ID} ${SEX}
                ChunkAndUpload ${SAMPLE_ID}
                
                touch ${SAMPLE_ID}.done
                gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
                DelocalizeSample ${SAMPLE_ID}
                ls -laht 1>&2
            done 4< ${FAMILY_ID}.samples.txt
            
            gcloud storage cp ${FAMILY_ID}_family.csv ~{remote_outdir}/
            rm -f ${FAMILY_ID}.samples.txt ${FAMILY_ID}_present.vcf.gz* ${FAMILY_ID}_family.csv
        done 3< family_ids.txt
    >>>
    
    output {
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
