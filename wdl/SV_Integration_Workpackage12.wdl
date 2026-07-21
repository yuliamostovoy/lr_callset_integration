version 1.0


# Performs a simple bcftools merge of all the ultralong or BND intra-sample
# VCFs. Stores in output per-chromosome VCFs that should then be split for
# parallel truvari collapse.
#
workflow SV_Integration_Workpackage12 {
    input {
        File sample_ids
        String suffix = "ultralong"
        String chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        
        String remote_indir
        Int n_expected_samples
        
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sample_ids: "Specifies the order of the samples to use in bcftools merge."
        remote_indir: "Without final slash. Must contain one sample-specific BCF and CSI pair per sample in sample_ids, named SAMPLE_ID_suffix.bcf and SAMPLE_ID_suffix.bcf.csi."
        remote_outdir: "Without final slash"
        suffix: "Denoting the type of intra-sample VCFs we want to merge: 'ultralong' or 'bnd'."
    }
    
    call Impl {
        input:
            sample_ids = sample_ids,
            suffix = suffix,
            chromosomes = chromosomes,
            
            remote_indir = remote_indir,
            n_expected_samples = n_expected_samples,
            
            remote_outdir = remote_outdir,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Historical performance on 12'680 samples, 15x, GRCh38, HDD, ultralong VCFs:
#
# TOOL                           CPU     RAM     TIME
# gcloud storage cp                                3m            // Whole genome
# bcftools merge level 1        300%    600M      10s            // Whole genome
# bcftools norm level 1         300%    300M      10s            // Whole genome
# bcftools merge level 2        200%    3.5G       4m            // Per chr
# bcftools norm level 2         250%      4G       2m            // Per chr
#
# Peak disk usage (all input files): 10G
#
task Impl {
    input {
        File sample_ids
        String suffix
        String chromosomes
        
        String remote_indir
        Int n_expected_samples
        
        Int n_files_per_merge = 100
        String remote_outdir
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number = 4
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeFiles() {
            local N_SAMPLES=$(wc -l < ~{sample_ids})
            if [ ${N_SAMPLES} -ne ~{n_expected_samples} ]; then
                echo "ERROR: sample_ids has ${N_SAMPLES} samples != ~{n_expected_samples}"
                exit 1
            fi

            rm -f all_remote_files.txt uri_list.txt
            while read -u 5 SAMPLE_ID; do
                echo "~{remote_indir}/${SAMPLE_ID}_~{suffix}.bcf" >> uri_list.txt
                echo "~{remote_indir}/${SAMPLE_ID}_~{suffix}.bcf.csi" >> uri_list.txt
            done 5< ~{sample_ids}
            while read -u 6 URI; do
                gcloud storage ls -l "${URI}" | grep -E '^[[:space:]]*[0-9]+' | sed 's/^[ ]*//' >> all_remote_files.txt
            done 6< uri_list.txt
        
            # Failing immediately if the files are too large WRT the available
            # disk. Otherwise the VM may get stuck forever, and this gets worse
            # with preemption.
            local AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
            AVAILABLE_GB=${AVAILABLE_GB%G}
            AVAILABLE_GB=${AVAILABLE_GB%.*}
            local REMOTE_GB=$(java -cp ~{docker_dir} SumFileSizes all_remote_files.txt)
            local SLACK_GB="5"
            REMOTE_GB=$(( ${REMOTE_GB} + ${SLACK_GB} ))
            if [ ${REMOTE_GB} -gt ${AVAILABLE_GB} ]; then
                echo "ERROR: the remote files are larger than the available disk space. Remote files + slack: ${REMOTE_GB}GB. Available disk: ${AVAILABLE_GB}GB."
                exit 1
            fi
        
            # Localizing all the single-sample VCFs.
            date 1>&2
            mkdir ./input_files/
            ${TIME_COMMAND} gcloud storage cp -I ./input_files/ < uri_list.txt
            date 1>&2
            local N_DOWNLOADED_SAMPLES=$(ls ./input_files/*.bcf | wc -l)
            if [ ${N_DOWNLOADED_SAMPLES} -lt ${N_SAMPLES} ]; then
                echo "ERROR: The number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is smaller than the number of samples specified (${N_SAMPLES})."
                exit 1
            elif [ ${N_DOWNLOADED_SAMPLES} -gt ${N_SAMPLES} ]; then
                echo "ERROR: The number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is larger than the number of samples specified (${N_SAMPLES})."
                exit 1
            fi
            df -h 1>&2
        }
        
        
        cat << 'END' > chunk_by_chr.sh
#!/bin/bash
INPUT_BCF=$1
CHROMOSOME=$2
mkdir -p ./${CHROMOSOME}/
bcftools view --output-type b ${INPUT_BCF} ${CHROMOSOME} --output ./${CHROMOSOME}/${INPUT_BCF}
bcftools index -f ./${CHROMOSOME}/${INPUT_BCF}
END
        chmod +x chunk_by_chr.sh
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        echo ~{chromosomes} | tr ',' '\n' > chromosomes.txt
        LocalizeFiles
        
        # Trivial "hierarchical" bcftools merge with just two steps.
        # Step 1: merging a few samples at a time over the whole genome.
        # Reheader each per-sample BCF to its sample_ids name before merging so
        # the cohort uses the same sample names as the main branch (WP3).
        rm -f list.txt
        while read -u 3 SAMPLE_ID; do
            echo ${SAMPLE_ID} > ${SAMPLE_ID}.sample_name.txt
            ${TIME_COMMAND} bcftools reheader --samples ${SAMPLE_ID}.sample_name.txt --output ./input_files/${SAMPLE_ID}_~{suffix}.reheader.bcf ./input_files/${SAMPLE_ID}_~{suffix}.bcf
            mv ./input_files/${SAMPLE_ID}_~{suffix}.reheader.bcf ./input_files/${SAMPLE_ID}_~{suffix}.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./input_files/${SAMPLE_ID}_~{suffix}.bcf
            echo ./input_files/${SAMPLE_ID}_~{suffix}.bcf >> list.txt
        done 3< ~{sample_ids}
        split -l ~{n_files_per_merge} -d -a 4 list.txt list_
        N_LIST_FILES=$(ls list_* | wc -l)
        for LIST_FILE in $(ls list_* | sort -V); do
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list ${LIST_FILE} --output-type b --output ${LIST_FILE}_merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_merged.bcf
            xargs --arg-file=${LIST_FILE} --max-lines=1 --max-procs=${N_THREADS} rm -f
            rm -f ${LIST_FILE}
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ${LIST_FILE}_merged.bcf --output ${LIST_FILE}_normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_normed.bcf
            rm -f ${LIST_FILE}_merged.bcf*
            ${TIME_COMMAND} xargs --arg-file=chromosomes.txt --max-lines=1 --max-procs=${N_THREADS} ./chunk_by_chr.sh ./${LIST_FILE}_normed.bcf
            rm -f ${LIST_FILE}_normed.bcf*
        done
        rm -rf ./input_files/
        
        # Step 2: merging all samples over each chromosome.
        rm -f files_list.txt
        while read -u 4 CHROMOSOME; do
            ls ./${CHROMOSOME}/*.bcf | sort -V > list.txt
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type b --output ./${CHROMOSOME}/merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./${CHROMOSOME}/merged.bcf
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ./${CHROMOSOME}/merged.bcf --output ./${CHROMOSOME}/normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./${CHROMOSOME}/normed.bcf
            mv ./${CHROMOSOME}/normed.bcf ${CHROMOSOME}.bcf
            mv ./${CHROMOSOME}/normed.bcf.csi ${CHROMOSOME}.bcf.csi
            echo "${CHROMOSOME}.bcf" >> files_list.txt
            echo "${CHROMOSOME}.bcf.csi" >> files_list.txt
            rm -rf ./${CHROMOSOME}/
        done 4< chromosomes.txt
        df -h 1>&2
        ls -laht 1>&2
        
        # Uploading
        date 1>&2
        cat files_list.txt | gcloud storage mv -I ~{remote_outdir}/
        date 1>&2

        # Completion signal for orchestrator ordering. Ignored standalone.
        echo "done" > wp12.signal
    >>>

    output {
        String done = read_string("wp12.signal")
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
