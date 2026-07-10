version 1.0


# STAGE 1 of the PAV whole-callset (SNV+indel+SV) `bcftools merge` pipeline.
#
# Given a Terra data table of samples (sample_id + PAV VCF URI), this workflow
# splits the cohort into batches of `batch_size` and, per batch, runs bcftools norm 
# and slices the normalized BCFs into chunks defined by a CSV.

workflow PAV_BcftoolsMerge_1_Split {
    input {
        Array[String] sample_ids
        Array[String] pav_vcfs

        File reference_fa
        File reference_fai
        File split_for_bcftools_merge_csv

        Int batch_size = 100

        String norm_remote_dir
        String remote_outdir

        Boolean left_align_and_split = true
        String check_ref = "w"

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.samples.sample_id`. Parallel to `pav_vcfs`."
        pav_vcfs: "From the data table, e.g. `this.samples.pav_vcf`. MUST be Array[String] (gs:// URIs) so files are streamed in-task rather than localized up front. Parallel to `sample_ids`."
        reference_fa: "Reference PAV was called against. Used for left-alignment."
        split_for_bcftools_merge_csv: "A partition covering all chromosomes: every line is a 0-based, half-open chunk `chrom,start,end`, sorted. Chunk id == 0-based line number. Same file the merge/concat stages use."
        batch_size: "Number of samples processed sequentially per VM."
        norm_remote_dir: "Without final slash. STABLE location for the normalized per-sample BCFs (`<norm_remote_dir>/<sample>.bcf`). Reuse the SAME path across runs so normalization is done once and re-chunking never re-normalizes."
        remote_outdir: "Without final slash. Receives `chunk_<i>/<sample>.bcf`, per-sample `<sample>.done`, and `sample_ids.txt`."
        left_align_and_split: "true: `bcftools norm -f ref -m -any` (left-align + split). false: `bcftools norm --do-not-normalize -m -any` (split multiallelics only, positions untouched)."
        check_ref: "`bcftools norm --check-ref` mode when left-aligning. 'w' warns and continues on REF mismatch; 'e' exits."
    }

    call MakeBatches {
        input:
            sample_ids = sample_ids,
            pav_vcfs = pav_vcfs,
            batch_size = batch_size,
            docker_image = docker_image
    }

    scatter (manifest in MakeBatches.batches) {
        call NormalizeBatch {
            input:
                manifest = manifest,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                norm_remote_dir = norm_remote_dir,
                left_align_and_split = left_align_and_split,
                check_ref = check_ref,
                docker_image = docker_image
        }
        call ChunkBatch {
            input:
                manifest = manifest,
                split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
                norm_remote_dir = norm_remote_dir,
                remote_outdir = remote_outdir,
                upstream_signal = NormalizeBatch.done,
                docker_image = docker_image
        }
    }

    call WriteSampleList {
        input:
            sample_ids = sample_ids,
            batch_signals = ChunkBatch.done,
            remote_outdir = remote_outdir,
            docker_image = docker_image
    }

    output {
        File sample_ids_file = WriteSampleList.sample_ids_file
        String done = WriteSampleList.done
    }
}


# Partitions the parallel (sample_id, pav_vcf) arrays into batches of
# `batch_size` rows, emitting one manifest TSV per batch. Independent of the
# split CSV, so batch composition is stable across chunk-size changes.
#
task MakeBatches {
    input {
        Array[String] sample_ids
        Array[String] pav_vcfs
        Int batch_size
        String docker_image
    }

    command <<<
        set -euxo pipefail

        N_IDS=$(wc -l < ~{write_lines(sample_ids)})
        N_VCFS=$(wc -l < ~{write_lines(pav_vcfs)})
        if [ ${N_IDS} -ne ${N_VCFS} ]; then
            echo "ERROR: sample_ids (${N_IDS}) and pav_vcfs (${N_VCFS}) have different lengths."
            exit 1
        fi

        # 2-column manifest: sample_id <TAB> pav_vcf_uri
        paste ~{write_lines(sample_ids)} ~{write_lines(pav_vcfs)} > all.tsv
        split --lines=~{batch_size} --numeric-suffixes=0 --suffix-length=6 all.tsv batch_
        ls batch_* 1>&2
    >>>

    output {
        Array[File] batches = glob("batch_*")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}


# Normalizes one batch of samples: localize -> reheader to sample_id ->
# left-align + split multiallelics. Idempotent via existence-check on
# `norm_remote_dir`. Does NOT reference the split CSV, so its result is reused
# unchanged across chunk-size retunes and pilot->full runs.
#
task NormalizeBatch {
    input {
        File manifest
        File reference_fa
        File reference_fai
        String norm_remote_dir

        Boolean left_align_and_split
        String check_ref

        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 64
        Int preemptible_number = 4
    }

    command <<<
        set -euxo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        if [ ~{true="1" false="0" left_align_and_split} -eq 1 ]; then
            NORM_ARGS="--fasta-ref ~{reference_fa} --check-ref ~{check_ref} --multiallelics -any"
        else
            NORM_ARGS="--do-not-normalize --multiallelics -any"
        fi

        function NormalizeSample() {
            local SAMPLE_ID=$1
            local PAV_URI=$2

            # Reuse: skip if the normalized BCF already exists. The .bcf is
            # uploaded LAST (after its index), so its presence implies both
            # objects are complete.
            local TEST=$( gcloud storage ls ~{norm_remote_dir}/${SAMPLE_ID}.bcf || echo "0" )
            if [ ${TEST} != "0" ]; then
                echo "Sample ${SAMPLE_ID} already normalized; skipping." 1>&2
                return 0
            fi

            ${TIME_COMMAND} gcloud storage cp ${PAV_URI} ./${SAMPLE_ID}.input.vcf.gz

            local N_IN_SAMPLES=$(bcftools query --list-samples ./${SAMPLE_ID}.input.vcf.gz | wc -l)
            if [ ${N_IN_SAMPLES} -ne 1 ]; then
                echo "ERROR: ${SAMPLE_ID} PAV VCF has ${N_IN_SAMPLES} samples (expected 1)."
                exit 1
            fi
            echo ${SAMPLE_ID} > ${SAMPLE_ID}.sample_name.txt
            bcftools reheader --samples ${SAMPLE_ID}.sample_name.txt --output ./${SAMPLE_ID}.reheader.vcf.gz ./${SAMPLE_ID}.input.vcf.gz
            rm -f ./${SAMPLE_ID}.input.vcf.gz

            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} ${NORM_ARGS} --output-type b ./${SAMPLE_ID}.reheader.vcf.gz --output ./${SAMPLE_ID}.bcf
            bcftools index --threads ${N_THREADS} ./${SAMPLE_ID}.bcf
            rm -f ./${SAMPLE_ID}.reheader.vcf.gz

            gcloud storage mv ./${SAMPLE_ID}.bcf.csi ~{norm_remote_dir}/${SAMPLE_ID}.bcf.csi
            gcloud storage mv ./${SAMPLE_ID}.bcf ~{norm_remote_dir}/${SAMPLE_ID}.bcf
            rm -f ./${SAMPLE_ID}.*
            df -h 1>&2
        }

        while read -u 3 -r SAMPLE_ID PAV_URI; do
            [ -z "${SAMPLE_ID}" ] && continue
            NormalizeSample "${SAMPLE_ID}" "${PAV_URI}"
        done 3< ~{manifest}

        echo "normalize_done" > batch.signal
    >>>

    output {
        String done = read_string("batch.signal")
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


# Slices one batch of already-normalized BCFs into chunks. Cheap (no norm, no
# reference). Idempotent per sample via `.done`. Depends on NormalizeBatch via
# `upstream_signal`.
#
task ChunkBatch {
    input {
        File manifest
        File split_for_bcftools_merge_csv
        String norm_remote_dir
        String remote_outdir
        String upstream_signal

        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 64
        Int preemptible_number = 4
    }

    command <<<
        set -euxo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        function ChunkSample() {
            local SAMPLE_ID=$1

            # Idempotency: skip if already chunked.
            local TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                echo "Sample ${SAMPLE_ID} already chunked; skipping." 1>&2
                return 0
            fi

            gcloud storage cp ~{norm_remote_dir}/${SAMPLE_ID}.bcf ./${SAMPLE_ID}.norm.bcf

            # Split into chunks. `--targets-file` uses POS-only overlap, so every
            # record is assigned to exactly one chunk by its final POS. The
            # interval file must end in `.bed` for 0-based interpretation,
            # matching the 0-based half-open CSV.
            local I=0
            local INTERVAL
            while read -u 4 INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > interval.bed
                bcftools view --threads ${N_THREADS} --targets-file interval.bed --output-type b ./${SAMPLE_ID}.norm.bcf --output ./${SAMPLE_ID}_chunk_${I}.bcf
                bcftools index --threads ${N_THREADS} ./${SAMPLE_ID}_chunk_${I}.bcf
                gcloud storage mv ./${SAMPLE_ID}_chunk_${I}.bcf ~{remote_outdir}/chunk_${I}/${SAMPLE_ID}.bcf
                gcloud storage mv ./${SAMPLE_ID}_chunk_${I}.bcf.csi ~{remote_outdir}/chunk_${I}/${SAMPLE_ID}.bcf.csi
                I=$(( ${I} + 1 ))
            done 4< ~{split_for_bcftools_merge_csv}
            rm -f ./${SAMPLE_ID}.norm.bcf*

            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            df -h 1>&2
        }

        while read -u 3 -r SAMPLE_ID PAV_URI; do
            [ -z "${SAMPLE_ID}" ] && continue
            ChunkSample "${SAMPLE_ID}"
        done 3< ~{manifest}

        echo "chunk_done" > batch.signal
    >>>

    output {
        String done = read_string("batch.signal")
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


# Writes the canonical, sorted sample list to the chunk output dir. Runs after
# all batches (via `batch_signals`) so the merge stage can rely on it existing.
#
task WriteSampleList {
    input {
        Array[String] sample_ids
        Array[String] batch_signals
        String remote_outdir
        String docker_image
    }

    command <<<
        set -euxo pipefail
        sort -u ~{write_lines(sample_ids)} > sample_ids.txt
        wc -l sample_ids.txt 1>&2
        gcloud storage cp sample_ids.txt ~{remote_outdir}/sample_ids.txt
        echo "split_done" > split.signal
    >>>

    output {
        File sample_ids_file = "sample_ids.txt"
        String done = read_string("split.signal")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}
