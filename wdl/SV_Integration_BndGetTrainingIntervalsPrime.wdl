version 1.0


# Like `SV_Integration_BndGetTrainingIntervals.wdl`, but uses breakpoints 
# derived directly from the SAM rather than SVIM-asm's BND calls as a source of 
# truth.
#
# Remark: the TPs in output are canonized, but the input is left untouched so it
# may still be not canonized.
#
workflow SV_Integration_BndGetTrainingIntervalsPrime {
    input {
        File samples_tsv

        Int breakpoint_max_distance = 500
        File reference_agp
        
        String remote_indir_query
        String remote_indir_truth
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        samples_tsv: "Format: ID,..."
        breakpoint_max_distance: "Max distance between a BND and an assembly breakpoint to mark TPs"
        reference_agp: "Reference AGP file."
        remote_indir_query: "Without final slash"
        remote_indir_truth: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            samples_tsv = samples_tsv,

            breakpoint_max_distance = breakpoint_max_distance,
            reference_agp = reference_agp,

            remote_indir_query = remote_indir_query,
            remote_indir_truth = remote_indir_truth,
            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 2-core, 8GB VM:
#
# TOOL                                      CPU%        RAM         TIME
# samtools view                              30%         3G          10m
# BndFilterWithAssemblyBreakpoints          100%        60M           1m
#
task Impl {
    input {
        File samples_tsv

        Int breakpoint_max_distance
        File reference_agp
        
        String remote_indir_query
        String remote_indir_truth
        String remote_outdir
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 4
        Int disk_size_gb = 20
        Int preemptible_number = 3
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{ram_size_gb} * 1024 - 500 ))

        samtools --version 1>&2
        bcftools --version 1>&2
        df -h 1>&2

        cat ~{samples_tsv} | tr '\t' ',' > samples.csv
        while read -u 3 LINE || [ -n "${LINE}" ]; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "1" )
            if [ "${TEST}" != "1" ]; then
                continue
            fi

            # Filtering
            gcloud storage cp ~{remote_indir_query}/${SAMPLE_ID}_bnd.vcf.'gz*' .
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M BndCanonize ${SAMPLE_ID}_bnd.vcf.gz > ${SAMPLE_ID}_bnd_canonized.vcf
            rm -f ${SAMPLE_ID}_bnd.vcf.gz*
            gcloud storage cp ~{remote_indir_truth}/${SAMPLE_ID}_breakpoints.csv .
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M BndFilterWithAssemblyBreakpoints2 ${SAMPLE_ID}_bnd_canonized.vcf ${SAMPLE_ID}_breakpoints.csv ~{breakpoint_max_distance} ~{reference_agp} | bcftools sort - --output-type z > ${SAMPLE_ID}_bnd_training.vcf.gz
            rm -f ${SAMPLE_ID}_bnd_canonized.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_bnd_training.vcf.gz
            bcftools filter --threads ${N_THREADS} --exclude 'INFO/TRUTH_BND_TYPE=4' --output-type z ${SAMPLE_ID}_bnd_training.vcf.gz --output ${SAMPLE_ID}_bnd_training_no_4.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_bnd_training_no_4.vcf.gz
            
            # Uploading and deallocating the sample
            gcloud storage mv ${SAMPLE_ID}_bnd_training'*.vcf.gz*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            rm -rf ${SAMPLE_ID}_*
            ls -laht 1>&2
        done 3< samples.csv
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
    }
}
