version 1.0


# Performs a bcftools merge of all the VCFs of a given chromosome chunk.
#
workflow SV_Integration_Workpackage3 {
    input {
        Int chunk_id
        File sample_ids
        
        String remote_indir
        
        Int merge_mode
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sample_ids: "Speficies the order of the samples to use in bcftools merge."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        merge_mode: "1: standard bcftools merge (CHROM,POS,REF,ALT). 2: bcftools merge by ID only."
    }
    
    call Impl {
        input:
            chunk_id = chunk_id,
            sample_ids = sample_ids,
            remote_indir = remote_indir,
            
            merge_mode = merge_mode,
            remote_outdir = remote_outdir,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, first 30 MB chunk of chr6:
#
# TOOL                          CPU     RAM     TIME
# gcloud storage ls                             30 s
# gcloud storage cp             280%    100 M    1 m
#
# Merge by CHROM,POS,REF,ALT:
# bcftools merge level 1        100%    300 M    3 s          // 100 files
# bcftools norm level 1         300%     50 M    1 s
# bcftools merge level 2        170%      1 G   20 m          // 127 files
# bcftools norm level 2         300%     11 G    6 m          // 16.5G in chr1
#
# Peak disk usage (all input files of chunk 0): 2 GB
#
# Merge by ID:
# bcftools merge level 1        400%    1.5G     5 s          // 100 files
# bcftools merge level 2        300%    2.5G    30 s          // 126 files
#
# Peak disk usage (all input files of chunk 0): 74G
#
task Impl {
    input {
        Int chunk_id
        File sample_ids
        String remote_indir
        
        Int merge_mode
        Int n_files_per_merge = 100
        String remote_outdir
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 16
        Int disk_size_gb = 100
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
        
        #
        function LocalizeChunkFiles() {
            local N_FILES
            local N_SAMPLES=$(cat ~{sample_ids} | wc -l)
            local SLACK_GB="5"
            
            # Ensuring that the input dataset has the expected number of
            # samples in the chunk.
            date 1>&2
            gcloud storage ls -l ~{remote_indir}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > remote_files.txt
            N_FILES=$(wc -l < remote_files.txt)
            N_FILES=$(( ${N_FILES} - 1 ))
            if [ ${N_FILES} -ne ${N_SAMPLES} ]; then
                echo "ERROR: input chunk has ${N_FILES} files != ${N_SAMPLES} samples in sample_ids"
                exit 1
            fi
            head -n ${N_FILES} remote_files.txt > all_remote_files.txt
            date 1>&2
        
            # Failing immediately if the files are too large WRT the available
            # disk. Otherwise the VM may get stuck forever, and this gets worse
            # with preemption.
            local AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
            AVAILABLE_GB=${AVAILABLE_GB%G}
            AVAILABLE_GB=${AVAILABLE_GB%.*}
            local REMOTE_GB=$(java -cp ~{docker_dir} SumFileSizes all_remote_files.txt)
            REMOTE_GB=$(( ${REMOTE_GB} + ${SLACK_GB} ))
            if [ ${REMOTE_GB} -gt ${AVAILABLE_GB} ]; then
                echo "ERROR: the remote files are larger than the available disk space. Remote files + slack: ${REMOTE_GB}GB. Available disk: ${AVAILABLE_GB}GB."
                exit 1
            fi
            rm -f remote_files.txt
        
            # Localizing all samples for the given chunk.
            date 1>&2
            mkdir ./input_files/
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/chunk_~{chunk_id}/'*' ./input_files/
            date 1>&2
            local N_DOWNLOADED_SAMPLES=$(ls ./input_files/*.bcf | wc -l)
            if [ ${N_DOWNLOADED_SAMPLES} -ne ${N_SAMPLES} ]; then
                echo "ERROR: The number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is different from the number of samples specified (${N_SAMPLES})."
                exit 1
            fi
            local SAMPLE_ID
            while read -u 3 SAMPLE_ID; do
                if [ ! -f ./input_files/${SAMPLE_ID}.bcf ]; then
                    echo "ERROR: Missing input BCF for sample ${SAMPLE_ID}."
                    exit 1
                fi
                echo ${SAMPLE_ID} > ${SAMPLE_ID}.sample_name.txt
                ${TIME_COMMAND} bcftools reheader --samples ${SAMPLE_ID}.sample_name.txt --output ./input_files/${SAMPLE_ID}.reheader.bcf ./input_files/${SAMPLE_ID}.bcf
                mv ./input_files/${SAMPLE_ID}.reheader.bcf ./input_files/${SAMPLE_ID}.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./input_files/${SAMPLE_ID}.bcf
            done 3< ~{sample_ids}
            df -h
        }
        
        
        # Trivial "hierarchical" merge with just two steps.
        #
        function MergeChunkFiles() {
            local MERGE_FLAG
            if [ ~{merge_mode} -eq 1 ]; then
                MERGE_FLAG="none"
            elif [ ~{merge_mode} -eq 2 ]; then
                MERGE_FLAG="id"
            else
                echo "ERROR: Merge mode unknown."
                exit 1
            fi
            
            # Step 1
            rm -f list.txt
            local SAMPLE_ID
            while read -u 4 SAMPLE_ID; do
                echo ./input_files/${SAMPLE_ID}.bcf >> list.txt
            done 4< ~{sample_ids}
            split -l ~{n_files_per_merge} -d -a 4 list.txt list_
            local N_LIST_FILES=$(ls list_* | wc -l)
            local LIST_FILE
            for LIST_FILE in $(ls list_* | sort -V); do
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge ${MERGE_FLAG} --file-list ${LIST_FILE} --output-type b --output ${LIST_FILE}_merged.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_merged.bcf
                df -h 1>&2
                xargs --arg-file=${LIST_FILE} --max-lines=1 --max-procs=${N_THREADS} rm -f
                if [ ~{merge_mode} -eq 1 ]; then
                    ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ${LIST_FILE}_merged.bcf --output ${LIST_FILE}_normed.bcf
                    ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_normed.bcf
                    df -h 1>&2
                    rm -f ${LIST_FILE}_merged.bcf* ; mv ${LIST_FILE}_normed.bcf ${LIST_FILE}_merged.bcf ; mv ${LIST_FILE}_normed.bcf.csi ${LIST_FILE}_merged.bcf.csi
                fi
            done
            
            # Step 2
            ls list_*.bcf | sort -V > list.txt
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge ${MERGE_FLAG} --file-list list.txt --output-type b --output ~{chunk_id}_merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ~{chunk_id}_merged.bcf
            df -h 1>&2
            xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} rm -f
            ls -laht 1>&2
            
            # Making sure no multiallelic record is passed downstream. This is
            # not needed when merging by ID, by construction.
            if [ ~{merge_mode} -eq 1 ]; then
                ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ~{chunk_id}_merged.bcf --output ~{chunk_id}_normed.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ~{chunk_id}_normed.bcf
                rm -f ~{chunk_id}_merged.bcf* ; mv ~{chunk_id}_normed.bcf ~{chunk_id}_merged.bcf ; mv ~{chunk_id}_normed.bcf.csi ~{chunk_id}_merged.bcf.csi
            fi
            
            # Removing records that are REF in all samples. This is not needed
            # in the standard merge, since at that step of the pipeline every
            # input record is ALT in some sample by construction.
            if [ ~{merge_mode} -eq 2 ]; then
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type b ~{chunk_id}_merged.bcf --output ~{chunk_id}_cleaned.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ~{chunk_id}_cleaned.bcf
                local N_RECORDS=$(bcftools index --nrecords ~{chunk_id}_merged.bcf)
                local N_ALT_RECORDS=$(bcftools index --nrecords ~{chunk_id}_cleaned.bcf)
                local PERCENT=$( echo "scale=2; 100 * ${N_ALT_RECORDS} / ${N_RECORDS}" | bc )
                echo "${N_ALT_RECORDS},${N_RECORDS},${PERCENT},Number of cohort-VCF records that are marked as ALT in >=1 sample by kanpig" > ~{chunk_id}_n_alt.csv
                rm -f ~{chunk_id}_merged.bcf* ; mv ~{chunk_id}_cleaned.bcf ~{chunk_id}_merged.bcf ; mv ~{chunk_id}_cleaned.bcf.csi ~{chunk_id}_merged.bcf.csi
            fi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        LocalizeChunkFiles
        MergeChunkFiles
        gcloud storage mv ~{chunk_id}_merged.bcf ~{remote_outdir}/chunk_~{chunk_id}.bcf
        gcloud storage mv ~{chunk_id}_merged.bcf.csi ~{remote_outdir}/chunk_~{chunk_id}.bcf.csi
        if [ ~{merge_mode} -eq 2 ]; then
            gcloud storage mv ~{chunk_id}_n_alt.csv ~{remote_outdir}/
        fi
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
