version 1.0


# Computes `truvari collapse` chunks from all the `bcftools merge` chunks of a
# chromosome.
#
workflow SV_Integration_Workpackage6 {
    input {
        String chromosome_id
        String bcftools_chunks
        
        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        
        String remote_indir
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        bcftools_chunks: "Comma-separated and sorted integers. Chunks are assumed to be sorted by POS."
        truvari_chunk_min_records: "Min number of records per output chunk"
        truvari_collapse_refdist: "The actual collapse downstream will run `truvari collapse --refdist X`, where X is this value."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chromosome_id = chromosome_id,
            bcftools_chunks = bcftools_chunks,
            truvari_chunk_min_records = truvari_chunk_min_records,
            truvari_collapse_refdist = truvari_collapse_refdist,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999:
#
# TOOL                                        CPU     RAM     TIME
# bcftools concat                             12%     20M     1h
#
# truvari divide --buffer 500 --min 1769     100%     12G     6h   // 562 chunks
# truvari divide --min 88468                 100%  >=256G
#
# First speedup strategy:
# TruvariDivide.java 2000                    100%    600M     2h   // 524 chunks
# Conversion GZ->BGZ                         800%     10M    16m
#
# Second speedup strategy:
# bcftools query                             100%      6G     7m
# TruvariDivide2.java 1000 2000              100%    300M     5s   // 524 chunks
# bcftools view                              800%      4G    30m
#
# Output of truvari divide (each chunk is ~5 MB):
# count      562.000000
# mean      3148.354093
# std       2623.478656
# min         68.000000
# 25%       1849.000000
# 50%       2070.500000
# 75%       3062.250000
# max      19218.000000
#
# Remark: as a third speedup strategy, we could create truvari collapse chunks
# from each bcftools merge chunk in parallel, and then merge the first/last
# truvari collapse chunks of consecutive bcftools merge chunks. This might be
# even faster.
# 
task Impl {
    input {
        String chromosome_id
        String bcftools_chunks
        
        Int truvari_chunk_min_records
        Int truvari_collapse_refdist
        
        String remote_indir
        String remote_outdir
        Int consistency_checks = 1
        
        String docker_image
        Int n_cpu = 8
        Int ram_size_gb = 12
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        
        
        # Localizing all the bcftools merge chunks of the chromosome
        rm -f uri_list.txt file_list.txt
        for CHUNK in $(echo ~{bcftools_chunks} | tr ',' ' '); do
            echo ~{remote_indir}/chunk_${CHUNK}.bcf >> uri_list.txt
            echo ~{remote_indir}/chunk_${CHUNK}.bcf.csi >> uri_list.txt
            echo chunk_${CHUNK}.bcf >> file_list.txt
        done
        date 1>&2
        cat uri_list.txt | gcloud storage cp -I .
        date 1>&2
        df -h 1>&2
        rm -f uri_list.txt
        
        # Concatenating all the bcftools merge chunks to build a whole-
        # chromosome VCF. This is because a truvari collapse chunk may straddle
        # multiple bcftools merge chunks.
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list.txt --output-type b --output ~{chromosome_id}.bcf
        bcftools index --threads ${N_THREADS} -f ~{chromosome_id}.bcf
        df -h 1>&2
        rm -f chunk_*.bcf* file_list.txt
        N_RECORDS=$(bcftools index --nrecords ~{chromosome_id}.bcf.csi)
        if [ ~{consistency_checks} -eq 1 ]; then
            ${TIME_COMMAND} bcftools query --format '%ID\n' ~{chromosome_id}.bcf > ids_truth.txt
        fi
        
        # Chunking the chromosome for truvari collapse
        ${TIME_COMMAND} bcftools query --format '%POS\t%REF\t%ALT\n' ~{chromosome_id}.bcf > pos_ref_alt.tsv
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G TruvariDivide2 pos_ref_alt.tsv ~{truvari_collapse_refdist} ~{truvari_chunk_min_records} ~{chromosome_id} ${N_RECORDS} > regions.txt
        rm -f pos_ref_alt.tsv
        cat << 'END' > chunk_by_region.sh
#!/bin/bash
INPUT_BCF=$1
REGION=$2
CHUNK_ID=$3
bcftools view --regions ${REGION} --regions-overlap pos --output-type b ${INPUT_BCF} --output chunk_${CHUNK_ID}.bcf
bcftools index -f chunk_${CHUNK_ID}.bcf
df -h 1>&2
END
        chmod +x chunk_by_region.sh
        ${TIME_COMMAND} xargs --arg-file=regions.txt --max-lines=1 --max-procs=${N_THREADS} ./chunk_by_region.sh ~{chromosome_id}.bcf
        ls -laht 1>&2
        df -h  1>&2
        rm -f ~{chromosome_id}.bcf*
        
        # Simple consistency checks
        if [ ~{consistency_checks} -eq 1 ]; then
            N_RECORDS_CHUNKED="0"
            for FILE in $(ls chunk_*.bcf.csi | sort -V); do
                N=$( bcftools index --nrecords ${FILE} )
                N_RECORDS_CHUNKED=$(( ${N_RECORDS_CHUNKED} + ${N} ))
            done
            if [ ${N_RECORDS_CHUNKED} -ne ${N_RECORDS} ]; then
                echo "ERROR: The truvari collapse chunks contain ${N_RECORDS_CHUNKED} total records, but the chromosome VCF contains ${N_RECORDS} records."
                exit 1
            fi
            rm -f ids_test.txt
            for FILE in $(ls chunk_*.bcf | sort -V); do
                bcftools query --format '%ID\n' ${FILE} >> ids_test.txt
            done
            diff --brief ids_test.txt ids_truth.txt
        fi
        
        # Uploading
        ls chunk_*.bcf* > file_list.txt
        cat file_list.txt | gcloud storage cp -I ~{remote_outdir}/~{chromosome_id}/
        gcloud storage cp regions.txt ~{remote_outdir}/~{chromosome_id}/
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
