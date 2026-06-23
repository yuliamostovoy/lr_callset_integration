version 1.0


# Splits ultralong and BND whole-chromosome VCFs into chunks for running
# truvari collapse in parallel downstream.
#
workflow SV_Integration_Workpackage13 {
    input {
        String chromosome_id
        String suffix
        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        Int consistency_checks
        
        String remote_indir
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        suffix: "Denoting the type of VCF we want to split: 'ultralong' or 'bnd'."
        truvari_chunk_min_records: "Min number of records per output chunk"
        truvari_collapse_refdist: "The actual collapse downstream will run `truvari collapse --refdist X`, where X is this value."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chromosome_id = chromosome_id,
            suffix = suffix,
            truvari_chunk_min_records = truvari_chunk_min_records,
            truvari_collapse_refdist = truvari_collapse_refdist,
            consistency_checks = consistency_checks,
            
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, HDD, ultralong VCFs:
#
# TOOL                           CPU     RAM     TIME
# TruvariDivide2Ultralong
# xargs bcftools view
#
# Peak disk usage (all input files): 10G
#
# Remark: `TruvariDivide2Ultralong.java` tends to create few chunks in the
# ultralong cohort VCF in practice. It might be more effective to partition the
# VCF by SVTYPE before splitting, to reduce overlaps.
#
task Impl {
    input {
        String chromosome_id
        String suffix
        Int truvari_chunk_min_records
        Int truvari_collapse_refdist
        Int consistency_checks
        
        String remote_indir
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        cat << 'END' > chunk_by_region.sh
#!/bin/bash
REGION=$1
CHUNK_ID=$2
bcftools view --regions ${REGION} --regions-overlap pos --output-type b ~{chromosome_id}.bcf --output chunk_${CHUNK_ID}.bcf
bcftools index -f chunk_${CHUNK_ID}.bcf
END
        chmod +x chunk_by_region.sh
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        gcloud storage cp ~{remote_indir}/~{chromosome_id}.'bcf*' .
        
        # Splitting
        N_RECORDS=$(bcftools index --nrecords ~{chromosome_id}.bcf)
        if [ ~{consistency_checks} -eq 1 ]; then
            ${TIME_COMMAND} bcftools query --format '%ID\n' ~{chromosome_id}.bcf > ids_truth.txt
        fi
        ${TIME_COMMAND} bcftools query --format '%POS\t%REF\t%ALT\n' ~{chromosome_id}.bcf > pos_ref_alt.tsv
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G TruvariDivide2Ultralong pos_ref_alt.tsv ~{truvari_collapse_refdist} ~{truvari_chunk_min_records} ~{chromosome_id} ${N_RECORDS} ~{suffix} > regions.txt
        rm -f pos_ref_alt.tsv
        ${TIME_COMMAND} xargs --arg-file=regions.txt --max-lines=1 --max-procs=${N_THREADS} ./chunk_by_region.sh
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
        ${TIME_COMMAND} gcloud storage mv 'chunk_*.bcf*' ~{remote_outdir}/~{chromosome_id}/
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
