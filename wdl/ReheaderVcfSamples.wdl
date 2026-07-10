version 1.0

# Renames the samples in a VCF according to a cross-reference table.
#
# The cross-reference is a 2-column TSV: column 1 = original sample name,
# column 2 = new sample name. `bcftools reheader -s` accepts this pair form
# directly: it matches each row by the OLD name and renames in place, so the
# original header order is preserved automatically and extra rows (samples not
# in this VCF) are ignored.
#
# The one sharp edge of `-s` is that a VCF sample with NO row in the table is
# silently left with its old name rather than erroring. To catch that, we first
# assert every sample in the VCF header has a mapping and fail loudly if not.
#
# `bcftools reheader` only rewrites the header line, not the records, so it is
# fast and cheap regardless of VCF size.
workflow ReheaderVcfSamples {
    input {
        File input_vcf
        File crossref_tsv
        String prefix = basename(basename(input_vcf, ".gz"), ".vcf") + ".reheadered"
        RuntimeAttr? runtime_attr_override
    }
    parameter_meta {
        input_vcf: "VCF to reheader (.vcf or .vcf.gz)."
        crossref_tsv: "2-column TSV: col1 = original sample name, col2 = new sample name. Must contain a mapping for every sample present in the VCF; extra rows are ignored."
        prefix: "Output file prefix. Result is written as <prefix>.vcf.gz."
    }

    call Reheader {
        input:
            input_vcf = input_vcf,
            crossref_tsv = crossref_tsv,
            prefix = prefix,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File reheadered_vcf = Reheader.reheadered_vcf
        File reheadered_vcf_index = Reheader.reheadered_vcf_index
    }
}

task Reheader {
    input {
        File input_vcf
        File crossref_tsv
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + 3 * ceil(size(input_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # Sample names currently in the VCF header, in header order.
        bcftools query -l ~{input_vcf} > old_samples.txt
        if [ ! -s old_samples.txt ]; then
            echo "ERROR: no samples found in the VCF header." 1>&2
            exit 1
        fi

        # Normalize the cross-reference: strip carriage returns (CRLF safety) and
        # drop blank / short lines so the map is clean.
        sed 's/\r$//' ~{crossref_tsv} | awk -F'\t' 'NF>=2 && $1!=""' > crossref.clean.tsv

        # Guard against `-s`'s silent pass-through: every VCF sample MUST have a
        # mapping, otherwise it would be left with its old name unnoticed.
        awk -F'\t' '
            NR==FNR { map[$1]=1; next }
            !($1 in map) { print "ERROR: no mapping in cross-reference for sample: " $1 > "/dev/stderr"; err=1 }
            END { if (err) exit 1 }
        ' crossref.clean.tsv old_samples.txt

        # Reject a table that maps two different samples to the same new name,
        # which would yield an invalid (duplicate-sample) VCF header.
        cut -f2 crossref.clean.tsv | sort | uniq -d > dup_new_names.txt
        if [ -s dup_new_names.txt ]; then
            echo "ERROR: cross-reference maps multiple samples to the same new name:" 1>&2
            cat dup_new_names.txt 1>&2
            exit 1
        fi

        # Rename by old->new pairs; header order is preserved, records untouched.
        bcftools reheader --samples crossref.clean.tsv --output ~{prefix}.vcf.gz ~{input_vcf}
        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File reheadered_vcf = "~{prefix}.vcf.gz"
        File reheadered_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-utils-basic:latest"
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
