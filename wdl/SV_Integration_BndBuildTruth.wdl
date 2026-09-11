version 1.0


# Builds a BND truthset directly from diploid-assembly-to-reference BAMs 
# (without passing through an SV caller). The truth is a CSV file that encodes
# breakpoint locations.
#
workflow SV_Integration_BndBuildTruth {
    input {
        String sample_id
        File hap1_bam
        File hap2_bam

        Int max_adjacency_distance = 1000
        Int min_violation_distance = 100000
        Int chromosome_mode = 0
        Int print_chain_start_end = 1
        Int containment_slack_bp = 1000
        Int min_internal_sv_length = 100000
        Int min_distance_from_contig_end = 1000

        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        max_adjacency_distance: "Max distance (on an assembled contig) between two alignments for them to be considered adjacent. 1kbp seems a good value based on a histogram of nearest-neighbor distances."
        min_violation_distance: "Min distance (on the same reference chr) between two alignments (that are adjacent on some contig) for them to be considered a colinearity violation. Any setting will capture some ultralong DELs."
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            hap1_bam = hap1_bam,
            hap2_bam = hap2_bam,

            max_adjacency_distance = max_adjacency_distance,
            min_violation_distance = min_violation_distance,
            chromosome_mode = chromosome_mode,
            print_chain_start_end = print_chain_start_end,
            containment_slack_bp = containment_slack_bp,
            min_internal_sv_length = min_internal_sv_length,
            min_distance_from_contig_end = min_distance_from_contig_end,

            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 2-core, 8GB VM:
#
# TOOL                                      CPU%        RAM         TIME
# samtools sort                             100%         6G           4m
# AssemblySam2Breakpoints2                   30%         1G           7m
#
task Impl {
    input {
        String sample_id
        File hap1_bam
        File hap2_bam

        Int max_adjacency_distance
        Int min_violation_distance
        Int chromosome_mode
        Int print_chain_start_end
        Int containment_slack_bp
        Int min_internal_sv_length
        Int min_distance_from_contig_end

        String remote_outdir
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 20
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
        EFFECTIVE_RAM_MB=$(( (~{ram_size_gb} - 1) * 1024 ))
        RAM_PER_THREAD_MB=$(( ${EFFECTIVE_RAM_MB} / ${N_THREADS} ))


        samtools --version 1>&2
        df -h 1>&2

        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} -n -O SAM -o hap1.sam ~{hap1_bam}
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} -n -O SAM -o hap2.sam ~{hap2_bam}
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${RAM_PER_THREAD_MB}M AssemblySam2Breakpoints2 hap1.sam ~{max_adjacency_distance} ~{min_violation_distance} ~{chromosome_mode} ~{print_chain_start_end} ~{containment_slack_bp} ~{min_internal_sv_length} ~{min_distance_from_contig_end} > ~{sample_id}_breakpoints1.csv &
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${RAM_PER_THREAD_MB}M AssemblySam2Breakpoints2 hap2.sam ~{max_adjacency_distance} ~{min_violation_distance} ~{chromosome_mode} ~{print_chain_start_end} ~{containment_slack_bp} ~{min_internal_sv_length} ~{min_distance_from_contig_end} > ~{sample_id}_breakpoints2.csv &
        wait
        cat ~{sample_id}_breakpoints1.csv ~{sample_id}_breakpoints2.csv | sort -t , -k1,1 -k2,2n | uniq > ~{sample_id}_breakpoints.csv
        gcloud storage mv ~{sample_id}_breakpoints.csv ~{remote_outdir}/
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
