version 1.0

workflow TRGTMergeCohort {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String prefix = "TRGT_merged"
        File reference_fasta
        File? reference_fai
        String trgt_merge_params = ""
        RuntimeAttr? runtime_attr_override
    }

    call MergeTRGT {
        input:
            input_vcfs = input_vcfs,
            input_vcf_indices = input_vcf_indices,
            prefix = prefix,
            reference_fasta = reference_fasta,
            reference_fai = reference_fai,
            trgt_merge_params = trgt_merge_params,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File merged_vcf = MergeTRGT.merged_vcf
        File merged_vcf_index = MergeTRGT.merged_vcf_index
    }
}

task MergeTRGT {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String prefix
        File reference_fasta
        File? reference_fai
        String trgt_merge_params

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + 3 * ceil(size(input_vcfs, "GiB")) + ceil(size(reference_fasta, "GiB"))

    command <<<
        set -euxo pipefail

        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        if [ ~{length(input_vcfs)} -ne ~{length(input_vcf_indices)} ]; then
            echo "ERROR: input_vcfs and input_vcf_indices must have the same length."
            exit 1
        fi

        cat > input_vcfs.list <<'EOF_VCFS'
~{sep="\n" input_vcfs}
EOF_VCFS

        trgt merge \
            --vcf-list input_vcfs.list \
            --output ~{prefix}.vcf.gz \
            --genome ~{reference_fasta} \
            --threads ${N_THREADS} \
            ~{trgt_merge_params} \
            -Oz

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-trgt:5.0"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
