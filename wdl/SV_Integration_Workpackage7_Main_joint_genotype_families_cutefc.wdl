version 1.0

# Builds family-specific candidate VCFs from the WP6 truvari-collapsed cohort
# VCF, re-genotypes each family member with cuteFC, and splits
# each sample's re-genotyped VCF into chunks for a downstream bcftools merge by
# ID.
#
# This is a drop-in sibling of SV_Integration_Workpackage7_Main_joint_genotype_families,
# which uses kanpig. This one is appropriate for variants of any size. 
# The family-candidate construction, per-sample chunking, and
# upload logic are identical to the kanpig workflow, so the outputs plug into
# the same downstream bcftools-merge-by-ID and WP8 concat steps.
#
# Differences vs the kanpig workflow:
#   - cuteFC does not take a ploidy BED: it genotypes purely from read support,
#     so per-sample sex/ploidy inputs are not used and are dropped. Sex-aware
#     genotyping on the sex chromosomes is therefore NOT applied here; this is a
#     known behavioral difference to keep in mind when comparing to kanpig.

workflow SV_Integration_Workpackage7_families_cutefc {
    input {
        Array[String] family_ids
        File ped
        Array[String] sample_ids
        Array[String] aligned_bais
        Array[String] aligned_bams
        File split_for_bcftools_merge_csv

        String remote_indir
        String remote_outdir
        String requester_pays_project = ""

        File reference_fa
        File reference_fai
        File autosomes_bed

        String cutefc_params_cohort = "--max_size -1 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        family_ids: "Family IDs to process. Each ID must match column 1 of ped."
        ped: "Standard 6-column PED: family_id, sample_id, paternal_id, maternal_id, sex, phenotype."
        sample_ids: "Sample IDs. Must be in the same order as aligned_bais and aligned_bams."
        aligned_bais: "Remote aligned BAM index paths. Must be in the same order as sample_ids and aligned_bams."
        aligned_bams: "Remote aligned BAM paths. Must be in the same order as sample_ids and aligned_bais."
        split_for_bcftools_merge_csv: "A partition that covers all chromosomes. Every line is a 0-based, half-open, consecutive chunk of a chromosome. Lines are assumed to be sorted."
        remote_indir: "Without final slash. Contains WP6 genome-wide truvari_collapsed.bcf and truvari_collapsed.bcf.csi."
        remote_outdir: "Without final slash. Output sample BCF chunks are written under chunk_N/sample.bcf."
        requester_pays_project: "Google Cloud project to bill for requester-pays BAM/BAI buckets. Leave blank for non-requester-pays buckets."
        cutefc_params_cohort: "Extra cuteFC force-calling parameters. --genotype and --sample are added automatically. --max_size -1 is what makes cuteFC report full-length (large) SVs."
    }

    call Impl {
        input:
            family_ids = family_ids,
            ped = ped,
            sample_ids = sample_ids,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            requester_pays_project = requester_pays_project,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            autosomes_bed = autosomes_bed,
            cutefc_params_cohort = cutefc_params_cohort,
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
        Array[String] aligned_bais
        Array[String] aligned_bams
        File split_for_bcftools_merge_csv

        String remote_indir
        String remote_outdir
        String requester_pays_project

        File reference_fa
        File reference_fai
        File autosomes_bed

        String cutefc_params_cohort
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
        N_PHYSICAL_CORES=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        # cuteFC does not run faster when given more threads than physical cores
        # (hyperthreading does not help it), so we cap it at the physical count.
        CUTEFC_N_THREADS=${N_PHYSICAL_CORES}
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        export RUST_BACKTRACE="full"
        GCLOUD_STORAGE_BILLING_FLAGS=""
        if [ -n "~{requester_pays_project}" ]; then
            GCLOUD_STORAGE_BILLING_FLAGS="--billing-project=~{requester_pays_project}"
        fi

        # cuteFC (via pysam) needs the .fai next to the .fa. WDL localizes them
        # to separate paths, so we co-locate them here with symlinks.
        ln -sf ~{reference_fa} ./reference.fa
        ln -sf ~{reference_fai} ./reference.fa.fai
        REFERENCE_FA="./reference.fa"


        # ----------------------- Steps of the pipeline ------------------------

        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2

            local ALIGNED_BAI=$(echo "${LINE}" | cut -f 2)
            local ALIGNED_BAM=$(echo "${LINE}" | cut -f 3)

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


        function Cutefc() {
            local FAMILY_ID=$1
            local SAMPLE_ID=$2

            # Single-sample view of the family-present sites. cuteFC force-calls
            # (re-genotypes) exactly these sites against the sample's reads,
            # regardless of the sites' original genotypes. --sample names the
            # output sample column, and IDs are preserved so the downstream
            # bcftools merge by ID lines records up across family members.
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${SAMPLE_ID} --output-type v ${FAMILY_ID}_present.vcf.gz --output ${SAMPLE_ID}_personalized.vcf

            mkdir -p ./${SAMPLE_ID}_cutefc_dir
            ${TIME_COMMAND} cuteFC --threads ${CUTEFC_N_THREADS} --genotype --sample ${SAMPLE_ID} ~{cutefc_params_cohort} -Ivcf ${SAMPLE_ID}_personalized.vcf ${SAMPLE_ID}_aligned.bam ${REFERENCE_FA} ${SAMPLE_ID}_out.vcf ./${SAMPLE_ID}_cutefc_dir
            rm -rf ./${SAMPLE_ID}_cutefc_dir ${SAMPLE_ID}_personalized.vcf

            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_out.vcf --output ${SAMPLE_ID}_cutefc.vcf.gz
            rm -f ${SAMPLE_ID}_out.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_cutefc.vcf.gz

            local N_RECORDS=$(bcftools index --nrecords ${SAMPLE_ID}_cutefc.vcf.gz)
            local N_PRESENT_RECORDS=$(bcftools query --format '%ID\n' --include 'GT="alt"' ${SAMPLE_ID}_cutefc.vcf.gz | wc -l)
            local PERCENT="0"
            if [ ${N_RECORDS} -gt 0 ]; then
                PERCENT=$(echo "scale=2; 100 * ${N_PRESENT_RECORDS} / ${N_RECORDS}" | bc)
            fi
            echo "${N_PRESENT_RECORDS},${N_RECORDS},${PERCENT},Number of records that are marked as ALT by cuteFC" >> ${SAMPLE_ID}_cutefc.csv

            local N_HETS_IN_AUTOSOMES=$(bcftools query --format '%ID\n' --include 'GT="het"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_cutefc.vcf.gz | wc -l)
            local N_PRESENT_RECORDS_IN_AUTOSOMES=$(bcftools query --format '%ID\n' --include 'GT="alt"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_cutefc.vcf.gz | wc -l)
            PERCENT="0"
            if [ ${N_PRESENT_RECORDS_IN_AUTOSOMES} -gt 0 ]; then
                PERCENT=$(echo "scale=2; 100 * ${N_HETS_IN_AUTOSOMES} / ${N_PRESENT_RECORDS_IN_AUTOSOMES}" | bc)
            fi
            echo "${N_HETS_IN_AUTOSOMES},${N_PRESENT_RECORDS_IN_AUTOSOMES},${PERCENT},Number of records in autosomes that are marked as HET by cuteFC" >> ${SAMPLE_ID}_cutefc.csv
        }


        function ChunkAndUpload() {
            local SAMPLE_ID=$1

            local i="0"
            local INTERVAL
            while read -u 5 INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type b ${SAMPLE_ID}_cutefc.vcf.gz --output ${SAMPLE_ID}_chunk_${i}.bcf
                bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_chunk_${i}.bcf
                gcloud storage cp ${SAMPLE_ID}_chunk_${i}.bcf ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf
                gcloud storage cp ${SAMPLE_ID}_chunk_${i}.bcf.csi ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf.csi
                i=$(( ${i} + 1 ))
            done 5< ~{split_for_bcftools_merge_csv}
            gcloud storage cp ${SAMPLE_ID}_cutefc.csv ~{remote_outdir}/
        }


        # ---------------------------- Main program ----------------------------

        cuteFC --version 1>&2 || echo "cuteFC returns an error code when called with --version"

        cat > sample_ids.txt <<'EOF_SAMPLE_IDS'
~{sep="\n" sample_ids}
EOF_SAMPLE_IDS
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
        N_ALIGNED_BAIS=$(wc -l < aligned_bais.txt)
        N_ALIGNED_BAMS=$(wc -l < aligned_bams.txt)
        if [ ${N_SAMPLE_IDS} -ne ${N_ALIGNED_BAIS} ] || [ ${N_SAMPLE_IDS} -ne ${N_ALIGNED_BAMS} ]; then
            echo "ERROR: sample_ids, aligned_bais, and aligned_bams must have the same length."
            echo "sample_ids=${N_SAMPLE_IDS}, aligned_bais=${N_ALIGNED_BAIS}, aligned_bams=${N_ALIGNED_BAMS}"
            exit 1
        fi
        paste sample_ids.txt aligned_bais.txt aligned_bams.txt > sample_metadata.tsv

        # Localizing the WP6 cohort VCF.
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
                LocalizeSample ${SAMPLE_ID} "${LINE}"
                Cutefc ${FAMILY_ID} ${SAMPLE_ID}
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
