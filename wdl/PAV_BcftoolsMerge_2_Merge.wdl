version 1.0


# STAGE 2 of the PAV whole-callset `bcftools merge` pipeline.
#
# Cohort-level `bcftools merge` for every genome chunk produced by stage 1. The
# number of chunks is DERIVED from the split CSV (chunk id == 0-based line
# number), and the workflow scatters over `range(n_chunks)` in a single
# submission -- no hand-listed chunk ids and no data table required here.
#
# Per chunk (mirrors SV_Integration_Workpackage3), on one VM:
#   1. Verify the chunk dir has exactly one BCF per sample.
#   2. Fail fast if the inputs won't fit on disk (`SumFileSizes`).
#   3. Localize all per-sample BCFs (already reheadered + indexed by stage 1, so
#      no reheader/re-index here).
#   4. Two-level hierarchical `bcftools merge --merge none` (batches of
#      `n_files_per_merge`, then a merge-of-merges).
#   5. Split any residual multiallelics with `bcftools norm --do-not-normalize
#      -m -any` to GUARANTEE biallelic-only output. `--do-not-normalize` means
#      NO left-alignment here: positions must not move, or a record could shift
#      out of its chunk. This step is phase-safe (split only, never join).
#   6. Upload `chunk_<i>.bcf` and a `chunk_<i>.done` marker (re-running the
#      submission then only redoes chunks that failed).
#
# `--merge none` copies each sample's GT verbatim, so PAV phasing is preserved.
#
workflow PAV_BcftoolsMerge_2_Merge {
    input {
        File split_for_bcftools_merge_csv
        String remote_indir
        String remote_outdir

        File? sample_ids_file

        Int n_files_per_merge = 100

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
        Int merge_disk_size_gb = 256
        Int merge_ram_size_gb = 32
        Int merge_n_cpu = 8
    }
    parameter_meta {
        split_for_bcftools_merge_csv: "The SAME CSV used in stage 1. Its line count is the number of chunks; chunk id == 0-based line number."
        remote_indir: "Without final slash. Stage 1's output dir (contains `chunk_<i>/<sample>.bcf` and `sample_ids.txt`)."
        remote_outdir: "Without final slash. Receives `chunk_<i>.bcf` and `chunk_<i>.done`."
        sample_ids_file: "Optional. Merge column order. If omitted, `sample_ids.txt` is read from `remote_indir` (produced by stage 1). The orchestrator passes stage 1's output directly to create a dependency and skip the extra download."
        n_files_per_merge: "Batch size for level-1 of the hierarchical merge."
    }

    Int n_chunks = length(read_lines(split_for_bcftools_merge_csv))

    if (!defined(sample_ids_file)) {
        call GetSampleList {
            input:
                remote_indir = remote_indir,
                docker_image = docker_image
        }
    }
    File sample_ids = select_first([sample_ids_file, GetSampleList.sample_ids_file])

    scatter (chunk_id in range(n_chunks)) {
        call MergeChunk {
            input:
                chunk_id = chunk_id,
                sample_ids = sample_ids,
                remote_indir = remote_indir,
                remote_outdir = remote_outdir,
                n_files_per_merge = n_files_per_merge,
                docker_image = docker_image,
                disk_size_gb = merge_disk_size_gb,
                ram_size_gb = merge_ram_size_gb,
                n_cpu = merge_n_cpu
        }
    }

    output {
        Array[String] done = MergeChunk.done
    }
}


task GetSampleList {
    input {
        String remote_indir
        String docker_image
    }
    command <<<
        set -euxo pipefail
        gcloud storage cp ~{remote_indir}/sample_ids.txt ./sample_ids.txt
        wc -l ./sample_ids.txt 1>&2
    >>>
    output {
        File sample_ids_file = "sample_ids.txt"
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}


# Performance reference for the SV-only pipeline this borrows from (12,680
# samples, 15x, GRCh38, first 30 MB chunk of chr6, merge by CHROM,POS,REF,ALT):
# peak input disk ~2 GB; whole normed chr1 ~16.5 GB. PAV-with-SNVs is ~2 orders
# of magnitude heavier per bp, hence smaller chunks (~5 Mbp) and larger disk.
# Pilot one dense + one sparse chunk and read the emitted `df -h` / `time`
# output before sizing the full fan-out.
#
task MergeChunk {
    input {
        Int chunk_id
        File sample_ids
        String remote_indir
        String remote_outdir

        Int n_files_per_merge
        Int slack_gb = 20

        String docker_image
        Int n_cpu = 8
        Int ram_size_gb = 32
        Int disk_size_gb = 256
        Int preemptible_number = 4
    }

    String docker_dir = "/callset_integration"

    command <<<
        set -euxo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # Idempotency: skip if already merged.
        TEST=$( gcloud storage ls ~{remote_outdir}/chunk_~{chunk_id}.done || echo "0" )
        if [ ${TEST} != "0" ]; then
            echo "Chunk ~{chunk_id} already done; skipping." 1>&2
            exit 0
        fi

        # ------------------------- Localize chunk files -----------------------
        N_SAMPLES=$(cat ~{sample_ids} | wc -l)

        date 1>&2
        gcloud storage ls -l ~{remote_indir}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > remote_files.txt
        N_FILES=$(wc -l < remote_files.txt)
        N_FILES=$(( ${N_FILES} - 1 ))   # drop the trailing TOTAL line
        if [ ${N_FILES} -ne ${N_SAMPLES} ]; then
            echo "ERROR: chunk ~{chunk_id} has ${N_FILES} BCFs != ${N_SAMPLES} samples in sample_ids."
            exit 1
        fi
        head -n ${N_FILES} remote_files.txt > all_remote_files.txt

        # Fail fast if the inputs won't fit on disk.
        AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
        AVAILABLE_GB=${AVAILABLE_GB%G}
        AVAILABLE_GB=${AVAILABLE_GB%.*}
        REMOTE_GB=$(java -cp ~{docker_dir} SumFileSizes all_remote_files.txt)
        REMOTE_GB=$(( ${REMOTE_GB} + ~{slack_gb} ))
        if [ ${REMOTE_GB} -gt ${AVAILABLE_GB} ]; then
            echo "ERROR: chunk ~{chunk_id} inputs + slack (${REMOTE_GB}GB) exceed available disk (${AVAILABLE_GB}GB). Increase merge_disk_size_gb or use a finer split CSV."
            exit 1
        fi
        rm -f remote_files.txt

        date 1>&2
        mkdir ./input_files/
        ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/chunk_~{chunk_id}/'*' ./input_files/
        date 1>&2
        N_DOWNLOADED=$(ls ./input_files/*.bcf | wc -l)
        if [ ${N_DOWNLOADED} -ne ${N_SAMPLES} ]; then
            echo "ERROR: downloaded ${N_DOWNLOADED} BCFs != ${N_SAMPLES} samples."
            exit 1
        fi
        # Stage 1 already reheadered to sample_id and uploaded indexes, so no
        # reheader / re-index is needed here.
        while read -u 3 SAMPLE_ID; do
            if [ ! -f ./input_files/${SAMPLE_ID}.bcf ]; then
                echo "ERROR: missing input BCF for sample ${SAMPLE_ID}."
                exit 1
            fi
        done 3< ~{sample_ids}
        df -h 1>&2

        # --------------------- Two-level hierarchical merge -------------------
        # Level 1: merge in batches of n_files_per_merge.
        rm -f list.txt
        while read -u 4 SAMPLE_ID; do
            echo ./input_files/${SAMPLE_ID}.bcf >> list.txt
        done 4< ~{sample_ids}
        split -l ~{n_files_per_merge} -d -a 4 list.txt list_
        for LIST_FILE in $(ls list_* | sort -V); do
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list ${LIST_FILE} --output-type b --output ${LIST_FILE}_merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_merged.bcf
            xargs --arg-file=${LIST_FILE} --max-lines=1 --max-procs=${N_THREADS} rm -f
            df -h 1>&2
        done

        # Level 2: merge the level-1 outputs.
        ls list_*.bcf | sort -V > list.txt
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type b --output merged.bcf
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f merged.bcf
        xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} rm -f
        rm -rf ./input_files/
        df -h 1>&2

        # Guarantee biallelic-only output. `--do-not-normalize` => split only, NO
        # left-alignment, so no record can move out of this chunk. Phase-safe.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b merged.bcf --output normed.bcf
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f normed.bcf
        rm -f merged.bcf*
        df -h 1>&2

        # Upload.
        gcloud storage mv normed.bcf ~{remote_outdir}/chunk_~{chunk_id}.bcf
        gcloud storage mv normed.bcf.csi ~{remote_outdir}/chunk_~{chunk_id}.bcf.csi
        touch chunk_~{chunk_id}.done
        gcloud storage mv chunk_~{chunk_id}.done ~{remote_outdir}/

        echo "chunk_~{chunk_id}" > chunk.signal
    >>>

    output {
        String done = read_string("chunk.signal")
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
