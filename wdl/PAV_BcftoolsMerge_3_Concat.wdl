version 1.0


# STAGE 3 of the PAV whole-callset `bcftools merge` pipeline.
#
# Concatenates the merged chunks into ONE BCF PER PRIMARY CONTIG (chr1-22, X, Y).
# There is deliberately no genome-wide concat: at 13k samples with SNVs a single
# genome-wide BCF is ~1 TB, which is unwieldy to build and consume. Per-contig
# BCFs are the pipeline's final deliverable; downstream tools that want finer
# granularity can also consume stage 2's per-chunk BCFs directly.
#
# Which chunks belong to which contig (and in what order) is DERIVED from the
# split CSV -- chunk id == 0-based line number, CSV sorted by position -- so no
# chunk list is passed by hand. Point it at stage 2's output dir.
#
# `bcftools concat --naive` is byte-level (no re-compression), so it is I/O
# bound; disk (inputs + output, ~2x the contig's merged size, dominated by the
# densest contig ~chr1) is the binding constraint, not CPU/RAM.
#
workflow PAV_BcftoolsMerge_3_Concat {
    input {
        File split_for_bcftools_merge_csv
        String remote_indir
        String remote_outdir

        Array[String] contigs = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]

        Array[String]? upstream_signal

        Int concat_disk_size_gb = 512

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        split_for_bcftools_merge_csv: "The SAME CSV used in stages 1-2. Maps chunk id (0-based line number) to contig and order."
        remote_indir: "Without final slash. Stage 2's output dir (contains `chunk_<i>.bcf`)."
        remote_outdir: "Without final slash. Receives one `<contig>/<contig>.bcf` per contig."
        contigs: "Primary contigs to emit, one BCF each. Defaults to chr1-22, chrX, chrY (no chrM by default; add it if wanted). A contig with no chunks in the CSV is skipped."
        upstream_signal: "Optional ordering dependency (the orchestrator passes stage 2's `done`). Otherwise unused."
        concat_disk_size_gb: "Disk per contig. Must hold inputs + output (~2x the contig's merged size); the densest contig (~chr1) dominates."
    }

    scatter (contig in contigs) {
        call ConcatContig {
            input:
                contig = contig,
                split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
                remote_indir = remote_indir,
                remote_outdir = remote_outdir,
                disk_size_gb = concat_disk_size_gb,
                docker_image = docker_image
        }
    }

    output {
        Array[String] per_contig_bcfs = ConcatContig.out_bcf
    }
}


task ConcatContig {
    input {
        String contig
        File split_for_bcftools_merge_csv
        String remote_indir
        String remote_outdir

        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 512
        Int preemptible_number = 4
    }

    command <<<
        set -euxo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # Chunk ids for this contig, in CSV (position) order. chunk id ==
        # 0-based CSV line number.
        awk -F, -v chr="~{contig}" '$1==chr {print NR-1}' ~{split_for_bcftools_merge_csv} > chunk_ids.txt
        N_CHUNKS=$(wc -l < chunk_ids.txt)
        if [ ${N_CHUNKS} -eq 0 ]; then
            echo "~{contig} has no chunks in the split CSV; skipping." 1>&2
            echo "NONE" > out_bcf.txt
            exit 0
        fi

        # Idempotency.
        TEST=$( gcloud storage ls ~{remote_outdir}/~{contig}/~{contig}.done || echo "0" )
        if [ ${TEST} != "0" ]; then
            echo "~{contig} already done; skipping." 1>&2
            echo "~{remote_outdir}/~{contig}/~{contig}.bcf" > out_bcf.txt
            exit 0
        fi

        # Localize this contig's chunks and build an ordered file list.
        rm -f uri_list.txt file_list.txt
        while read -u 3 ID; do
            echo ~{remote_indir}/chunk_${ID}.bcf >> uri_list.txt
            echo ~{remote_indir}/chunk_${ID}.bcf.csi >> uri_list.txt
            echo chunk_${ID}.bcf >> file_list.txt
        done 3< chunk_ids.txt
        ${TIME_COMMAND} cat uri_list.txt | gcloud storage cp -I .
        df -h 1>&2

        # Concatenate chunks (already position-ordered).
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list.txt --output-type b --output ~{contig}.bcf
        bcftools index --threads ${N_THREADS} -f ~{contig}.bcf
        rm -f chunk_*.bcf*
        df -h 1>&2

        gcloud storage cp ~{contig}.bcf ~{remote_outdir}/~{contig}/~{contig}.bcf
        gcloud storage cp ~{contig}.bcf.csi ~{remote_outdir}/~{contig}/~{contig}.bcf.csi
        touch ~{contig}.done
        gcloud storage mv ~{contig}.done ~{remote_outdir}/~{contig}/

        echo "~{remote_outdir}/~{contig}/~{contig}.bcf" > out_bcf.txt
    >>>

    output {
        String out_bcf = read_string("out_bcf.txt")
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
