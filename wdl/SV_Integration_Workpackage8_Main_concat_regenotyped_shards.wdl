version 1.0


# Combines the bcftools merge chunks of a re-genotyped chromosome, in the given
# order.
#
workflow SV_Integration_Workpackage11 {
    input {
        String chunk_ids
        
        String remote_indir
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        chunk_ids: "Comma-separated. Chunk order is assumed to reflect POS order, and every chunk is assumed to be sorted by POS."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chunk_ids = chunk_ids,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999:
#
# TOOL                          CPU     RAM     TIME
# bcftools concat               17%     20M     10s
# 
task Impl {
    input {
        String chunk_ids
        
        String remote_indir
        String remote_outdir
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 4
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        
        
        # Localizing
        for CHUNK in $(echo ~{chunk_ids} | tr ',' ' '); do
            echo ~{remote_indir}/chunk_${CHUNK}.bcf >> uri_list.txt
            echo ~{remote_indir}/chunk_${CHUNK}.bcf.csi >> uri_list.txt
            echo chunk_${CHUNK}.bcf >> file_list.txt
        done
        date 1>&2
        cat uri_list.txt | gcloud storage cp -I .
        date 1>&2
        df -h 1>&2
        
        # Concatenating
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list.txt --output-type b --output merged.bcf
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f merged.bcf
        df -h 1>&2
        
        # Uploading
        gcloud storage mv merged.'bcf*' ~{remote_outdir}/
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_number
    }
}
