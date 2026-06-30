version 1.0


# Builds trio-specific candidate VCFs from the WP8 truvari-collapsed cohort
# VCF, re-genotypes each complete mother/father/proband trio with kanpig trio,
# and splits selected sample outputs into chunks for downstream bcftools merge
# by ID. This is intended as a test-mode sibling of Workpackage9_families.
#
workflow SV_Integration_Workpackage9_trios {
    input {
        File family_ids
        File ped
        File sv_integration_sample_tsv
        File split_for_bcftools_merge_csv
        
        String remote_indir
        String remote_outdir
        String requester_pays_project = ""
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed
        
        String kanpig_params_trio = "--neighdist 500 --gpenalty 0.04"
        Boolean upload_parents = false
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        family_ids: "One family ID per line. Each ID must match column 1 of ped. Complete trios are inferred from child rows with nonzero paternal_id and maternal_id."
        ped: "Standard 6-column PED: family_id, sample_id, paternal_id, maternal_id, sex, phenotype."
        sv_integration_sample_tsv: "Sample metadata TSV. First four columns are sample_id, sex, aligned_bai, aligned_bam."
        split_for_bcftools_merge_csv: "A partition that covers all chromosomes. Every line is a 0-based, half-open, consecutive chunk of a chromosome. Lines are assumed to be sorted."
        remote_indir: "Without final slash. Contains WP8 genome-wide truvari_collapsed.bcf and truvari_collapsed.bcf.csi."
        remote_outdir: "Without final slash. Output sample BCF chunks are written under chunk_N/sample.bcf. By default only probands are uploaded; set upload_parents=true to upload parents too."
        requester_pays_project: "Google Cloud project to bill for requester-pays BAM/BAI buckets. Leave blank for non-requester-pays buckets."
    }
    
    call Impl {
        input:
            family_ids = family_ids,
            ped = ped,
            sv_integration_sample_tsv = sv_integration_sample_tsv,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            requester_pays_project = requester_pays_project,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            autosomes_bed = autosomes_bed,
            kanpig_params_trio = kanpig_params_trio,
            upload_parents = upload_parents,
            docker_image = docker_image
    }
    
    output {
    }
}


task Impl {
    input {
        File family_ids
        File ped
        File sv_integration_sample_tsv
        File split_for_bcftools_merge_csv
        
        String remote_indir
        String remote_outdir
        String requester_pays_project
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed
        
        String kanpig_params_trio
        Boolean upload_parents
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
        
        
        function BuildFamilyTrios() {
            local FAMILY_ID=$1
            
            awk -v family="${FAMILY_ID}" 'BEGIN { FS="[ \t]+"; OFS="\t" } $1==family && $2!="0" && $2!="." && $3!="0" && $3!="." && $4!="0" && $4!="." { print $2,$3,$4 }' ~{ped} | sort -u > ${FAMILY_ID}.trios.tsv
            local N_TRIOS=$(wc -l < ${FAMILY_ID}.trios.tsv)
            echo "${FAMILY_ID},${N_TRIOS},Number of complete trios found in PED" > ${FAMILY_ID}_trios.csv
            if [ ${N_TRIOS} -eq 0 ]; then
                echo "WARNING: family ${FAMILY_ID} has no complete trios in the PED."
                return
            fi
            
            local PROBAND_ID
            local FATHER_ID
            local MOTHER_ID
            while read -u 4 PROBAND_ID FATHER_ID MOTHER_ID; do
                ValidateSampleInMetadata ${FAMILY_ID} ${PROBAND_ID}
                ValidateSampleInMetadata ${FAMILY_ID} ${FATHER_ID}
                ValidateSampleInMetadata ${FAMILY_ID} ${MOTHER_ID}
            done 4< ${FAMILY_ID}.trios.tsv
        }
        
        
        function ValidateSampleInMetadata() {
            local FAMILY_ID=$1
            local SAMPLE_ID=$2
            
            if ! awk -v sample="${SAMPLE_ID}" 'BEGIN { FS="\t"; found=0 } $1==sample { found=1 } END { exit(found ? 0 : 1) }' sample_metadata.tsv; then
                echo "ERROR: sample ${SAMPLE_ID} from family ${FAMILY_ID} is missing from sv_integration_sample_tsv."
                exit 1
            fi
        }
        
        
        function BuildTrioCandidateVcf() {
            local PROBAND_ID=$1
            local FATHER_ID=$2
            local MOTHER_ID=$3
            
            printf "%s\n%s\n%s\n" ${PROBAND_ID} ${FATHER_ID} ${MOTHER_ID} > ${PROBAND_ID}.trio.samples.txt
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file ${PROBAND_ID}.trio.samples.txt --output-type z cohort.bcf --output ${PROBAND_ID}_trio_all.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${PROBAND_ID}_trio_all.vcf.gz
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z ${PROBAND_ID}_trio_all.vcf.gz --output ${PROBAND_ID}_trio_present.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${PROBAND_ID}_trio_present.vcf.gz
            rm -f ${PROBAND_ID}_trio_all.vcf.gz*
        }
        
        
        function KanpigTrio() {
            local PROBAND_ID=$1
            local FATHER_ID=$2
            local MOTHER_ID=$3
            local PROBAND_SEX=$4
            
            local KARYOTYPE
            if [ ${PROBAND_SEX} == "M" ]; then
                KARYOTYPE="XY"
            else
                KARYOTYPE="XX"
            fi
            
            # Remark: kanpig needs --sizemin >= --kmer.
            ${TIME_COMMAND} ~{docker_dir}/kanpig trio --threads $(( ${N_THREADS} - 1)) --sizemin 10 --sizemax ${INFINITY} ~{kanpig_params_trio} --reference ~{reference_fa} --XYploidy-bed ~{ploidy_bed_male} --XXploidy-bed ~{ploidy_bed_female} --karyotype ${KARYOTYPE} --input ${PROBAND_ID}_trio_present.vcf.gz --proband ${PROBAND_ID}_aligned.bam --father ${FATHER_ID}_aligned.bam --mother ${MOTHER_ID}_aligned.bam --out ${PROBAND_ID}_trio_out.vcf --proband-sample ${PROBAND_ID} --father-sample ${FATHER_ID} --mother-sample ${MOTHER_ID}
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${PROBAND_ID}_trio_out.vcf --output ${PROBAND_ID}_trio_kanpig.vcf.gz
            rm -f ${PROBAND_ID}_trio_out.vcf
            bcftools index --threads ${N_THREADS} -f -t ${PROBAND_ID}_trio_kanpig.vcf.gz
            
            local N_RECORDS=$(bcftools index --nrecords ${PROBAND_ID}_trio_kanpig.vcf.gz)
            echo "${PROBAND_ID},${FATHER_ID},${MOTHER_ID},${N_RECORDS},Trio records genotyped by kanpig trio" > ${PROBAND_ID}_trio_kanpig.csv
        }
        
        
        function SplitTrioSample() {
            local PROBAND_ID=$1
            local SAMPLE_ID=$2
            
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${SAMPLE_ID} --output-type z ${PROBAND_ID}_trio_kanpig.vcf.gz --output ${PROBAND_ID}_${SAMPLE_ID}_kanpig.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${PROBAND_ID}_${SAMPLE_ID}_kanpig.vcf.gz
            
            local N_RECORDS=$(bcftools index --nrecords ${PROBAND_ID}_${SAMPLE_ID}_kanpig.vcf.gz)
            local N_PRESENT_RECORDS=$(bcftools query --format '%ID\n' --include 'GT="alt"' ${PROBAND_ID}_${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            local PERCENT="0"
            if [ ${N_RECORDS} -gt 0 ]; then
                PERCENT=$(echo "scale=2; 100 * ${N_PRESENT_RECORDS} / ${N_RECORDS}" | bc)
            fi
            echo "${PROBAND_ID},${SAMPLE_ID},${N_PRESENT_RECORDS},${N_RECORDS},${PERCENT},Number of records that are marked as ALT by kanpig trio" >> ${PROBAND_ID}_${SAMPLE_ID}_kanpig.csv
            
            local N_HETS_IN_AUTOSOMES=$(bcftools query --format '%ID\n' --include 'GT="het"' --regions-file ~{autosomes_bed} --regions-overlap pos ${PROBAND_ID}_${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            local N_PRESENT_RECORDS_IN_AUTOSOMES=$(bcftools query --format '%ID\n' --include 'GT="alt"' --regions-file ~{autosomes_bed} --regions-overlap pos ${PROBAND_ID}_${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            PERCENT="0"
            if [ ${N_PRESENT_RECORDS_IN_AUTOSOMES} -gt 0 ]; then
                PERCENT=$(echo "scale=2; 100 * ${N_HETS_IN_AUTOSOMES} / ${N_PRESENT_RECORDS_IN_AUTOSOMES}" | bc)
            fi
            echo "${PROBAND_ID},${SAMPLE_ID},${N_HETS_IN_AUTOSOMES},${N_PRESENT_RECORDS_IN_AUTOSOMES},${PERCENT},Number of records in autosomes that are marked as HET by kanpig trio" >> ${PROBAND_ID}_${SAMPLE_ID}_kanpig.csv
            ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigWindows ${PROBAND_ID}_${SAMPLE_ID}_kanpig.vcf.gz | bgzip > ${PROBAND_ID}_${SAMPLE_ID}_kanpig.bed.gz
        }
        
        
        function ChunkAndUpload() {
            local PROBAND_ID=$1
            local SAMPLE_ID=$2
            
            local i="0"
            local INTERVAL
            while read -u 5 INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${PROBAND_ID}_${SAMPLE_ID}.bed
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file ${PROBAND_ID}_${SAMPLE_ID}.bed --regions-overlap pos --output-type b ${PROBAND_ID}_${SAMPLE_ID}_kanpig.vcf.gz --output ${PROBAND_ID}_${SAMPLE_ID}_chunk_${i}.bcf
                bcftools index --threads ${N_THREADS} -f ${PROBAND_ID}_${SAMPLE_ID}_chunk_${i}.bcf
                gcloud storage cp ${PROBAND_ID}_${SAMPLE_ID}_chunk_${i}.bcf ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf
                gcloud storage cp ${PROBAND_ID}_${SAMPLE_ID}_chunk_${i}.bcf.csi ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf.csi
                i=$(( ${i} + 1 ))
            done 5< ~{split_for_bcftools_merge_csv}
            gcloud storage cp ${PROBAND_ID}_${SAMPLE_ID}_kanpig.bed.gz ${PROBAND_ID}_${SAMPLE_ID}_kanpig.csv ~{remote_outdir}/
        }
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        ~{docker_dir}/kanpig --version 1>&2
        cp ~{sv_integration_sample_tsv} sample_metadata.tsv
        
        # Localizing the WP8 cohort VCF.
        ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/truvari_collapsed.'bcf*' .
        mv truvari_collapsed.bcf cohort.bcf
        mv truvari_collapsed.bcf.csi cohort.bcf.csi
        
        while read -u 3 FAMILY_ID; do
            if [ -z "${FAMILY_ID}" ] || [[ "${FAMILY_ID}" == \#* ]]; then
                continue
            fi
            
            BuildFamilyTrios ${FAMILY_ID}
            
            while read -u 4 PROBAND_ID FATHER_ID MOTHER_ID; do
                TEST=$(gcloud storage ls ~{remote_outdir}/${PROBAND_ID}.done || echo "0")
                if [ ${TEST} != "0" ]; then
                    continue
                fi
                
                PROBAND_LINE=$(awk -v sample="${PROBAND_ID}" 'BEGIN { FS="\t" } $1==sample { print; exit }' sample_metadata.tsv)
                FATHER_LINE=$(awk -v sample="${FATHER_ID}" 'BEGIN { FS="\t" } $1==sample { print; exit }' sample_metadata.tsv)
                MOTHER_LINE=$(awk -v sample="${MOTHER_ID}" 'BEGIN { FS="\t" } $1==sample { print; exit }' sample_metadata.tsv)
                PROBAND_SEX=$(echo "${PROBAND_LINE}" | cut -f 2)
                
                BuildTrioCandidateVcf ${PROBAND_ID} ${FATHER_ID} ${MOTHER_ID}
                LocalizeSample ${PROBAND_ID} "${PROBAND_LINE}"
                LocalizeSample ${FATHER_ID} "${FATHER_LINE}"
                LocalizeSample ${MOTHER_ID} "${MOTHER_LINE}"
                KanpigTrio ${PROBAND_ID} ${FATHER_ID} ${MOTHER_ID} ${PROBAND_SEX}
                
                SplitTrioSample ${PROBAND_ID} ${PROBAND_ID}
                ChunkAndUpload ${PROBAND_ID} ${PROBAND_ID}
                if ~{upload_parents}; then
                    SplitTrioSample ${PROBAND_ID} ${FATHER_ID}
                    ChunkAndUpload ${PROBAND_ID} ${FATHER_ID}
                    SplitTrioSample ${PROBAND_ID} ${MOTHER_ID}
                    ChunkAndUpload ${PROBAND_ID} ${MOTHER_ID}
                fi
                gcloud storage cp ${PROBAND_ID}_trio_kanpig.csv ~{remote_outdir}/
                
                touch ${PROBAND_ID}.done
                gcloud storage mv ${PROBAND_ID}.done ~{remote_outdir}/
                DelocalizeSample ${PROBAND_ID}
                DelocalizeSample ${FATHER_ID}
                DelocalizeSample ${MOTHER_ID}
                ls -laht 1>&2
            done 4< ${FAMILY_ID}.trios.tsv
            
            gcloud storage cp ${FAMILY_ID}_trios.csv ~{remote_outdir}/
            rm -f ${FAMILY_ID}.trios.tsv ${FAMILY_ID}_trios.csv
        done 3< ~{family_ids}
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
