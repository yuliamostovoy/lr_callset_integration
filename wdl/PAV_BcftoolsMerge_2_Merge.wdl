version 1.0


# STAGE 2 of the PAV whole-callset `bcftools merge` pipeline.
#
# Cohort-level `bcftools merge` for every genome chunk produced by stage 1. The
# number of chunks is DERIVED from the split CSV (chunk id == 0-based line
# number), and the workflow scatters over `range(n_chunks)` in a single
# submission -- no hand-listed chunk ids and no data table required here.
#
# Sourced from up to six independently-run stage-1 outputs, one per center/
# control set (bi, ha, bcm, uw, controls_15x, controls_30x) -- this mirrors
# SV_Integration_Workpackage5's InterCenter merge pattern exactly, including
# its overlap-resolution rule: a sample present in more than one source is
# resolved by LOCALIZATION ORDER (bi, then ha, then uw, then bcm, then the
# controls -- each later source silently overwrites an earlier source's file
# for the same sample_id), with one explicit override: `bi_samples_to_prefer_over_ha`
# lists sample_ids for which bi's copy is re-localized after ha, flipping the
# winner back to bi for exactly those samples. There is no equivalent override
# for uw/bcm/controls; whichever source is localized last always wins for a
# sample it shares with an earlier source.
#
# Per chunk (mirrors SV_Integration_Workpackage5's Impl task), on one VM:
#   1. Verify each source's chunk dir has exactly the expected number of BCFs
#      (n_expected_samples_<source>; controls may exceed their expected count,
#      which only warns).
#   2. Fail fast if the inputs won't fit on disk (`SumFileSizes`).
#   3. Localize all per-sample BCFs in source order (already reheadered +
#      indexed by stage 1, so no reheader/re-index here), applying the
#      bi_samples_to_prefer_over_ha override last.
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

        String remote_indir_bi
        String remote_indir_ha
        String remote_indir_bcm
        String remote_indir_uw
        String remote_indir_controls_15x
        String remote_indir_controls_30x

        Array[String] bi_samples_to_prefer_over_ha

        Int n_expected_samples_bi
        Int n_expected_samples_ha
        Int n_expected_samples_bcm
        Int n_expected_samples_uw
        Int n_expected_samples_controls_15x
        Int n_expected_samples_controls_30x

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
        remote_indir_bi: "Without final slash. bi's stage-1 output dir (contains `chunk_<i>/<sample>.bcf` and `sample_ids.txt`). Set n_expected_samples_bi=0 and this to anything if bi is unused."
        remote_indir_ha: "Without final slash. ha's stage-1 output dir. Set n_expected_samples_ha=0 if unused."
        remote_indir_bcm: "Without final slash. bcm's stage-1 output dir. Set n_expected_samples_bcm=0 if unused."
        remote_indir_uw: "Without final slash. uw's stage-1 output dir. Set n_expected_samples_uw=0 if unused."
        remote_indir_controls_15x: "Without final slash. 15x-controls stage-1 output dir. Set n_expected_samples_controls_15x=0 if unused."
        remote_indir_controls_30x: "Without final slash. 30x-controls stage-1 output dir. Set n_expected_samples_controls_30x=0 if unused."
        bi_samples_to_prefer_over_ha: "sample_ids present in BOTH bi and ha for which bi's copy should win instead of the default (ha wins on overlap). No effect on samples not present in both."
        n_expected_samples_bi: "Exact expected file count in bi's chunk dirs. 0 disables bi entirely."
        n_expected_samples_ha: "Exact expected file count in ha's chunk dirs. 0 disables ha entirely."
        n_expected_samples_bcm: "Exact expected file count in bcm's chunk dirs. 0 disables bcm entirely."
        n_expected_samples_uw: "Exact expected file count in uw's chunk dirs. 0 disables uw entirely."
        n_expected_samples_controls_15x: "Minimum expected file count in the 15x-controls chunk dirs (more than expected only warns). 0 disables this source entirely."
        n_expected_samples_controls_30x: "Minimum expected file count in the 30x-controls chunk dirs (more than expected only warns). 0 disables this source entirely."
        remote_outdir: "Without final slash. Receives `chunk_<i>.bcf` and `chunk_<i>.done`."
        sample_ids_file: "Optional. Merge column order (the UNION of sample_ids across all enabled sources). If omitted, it is derived by unioning `sample_ids.txt` from every source with a nonzero n_expected_samples_<source> (each produced by stage 1). The orchestrator passes stage 1's union directly to create a dependency and skip the extra download."
        n_files_per_merge: "Batch size for level-1 of the hierarchical merge."
    }

    Int n_chunks = length(read_lines(split_for_bcftools_merge_csv))

    if (!defined(sample_ids_file)) {
        call GetSampleList {
            input:
                remote_indir_bi = remote_indir_bi,
                remote_indir_ha = remote_indir_ha,
                remote_indir_bcm = remote_indir_bcm,
                remote_indir_uw = remote_indir_uw,
                remote_indir_controls_15x = remote_indir_controls_15x,
                remote_indir_controls_30x = remote_indir_controls_30x,
                n_expected_samples_bi = n_expected_samples_bi,
                n_expected_samples_ha = n_expected_samples_ha,
                n_expected_samples_bcm = n_expected_samples_bcm,
                n_expected_samples_uw = n_expected_samples_uw,
                n_expected_samples_controls_15x = n_expected_samples_controls_15x,
                n_expected_samples_controls_30x = n_expected_samples_controls_30x,
                docker_image = docker_image
        }
    }
    File sample_ids = select_first([sample_ids_file, GetSampleList.sample_ids_file])

    scatter (chunk_id in range(n_chunks)) {
        call MergeChunk {
            input:
                chunk_id = chunk_id,
                sample_ids = sample_ids,
                remote_indir_bi = remote_indir_bi,
                remote_indir_ha = remote_indir_ha,
                remote_indir_bcm = remote_indir_bcm,
                remote_indir_uw = remote_indir_uw,
                remote_indir_controls_15x = remote_indir_controls_15x,
                remote_indir_controls_30x = remote_indir_controls_30x,
                bi_samples_to_prefer_over_ha = bi_samples_to_prefer_over_ha,
                n_expected_samples_bi = n_expected_samples_bi,
                n_expected_samples_ha = n_expected_samples_ha,
                n_expected_samples_bcm = n_expected_samples_bcm,
                n_expected_samples_uw = n_expected_samples_uw,
                n_expected_samples_controls_15x = n_expected_samples_controls_15x,
                n_expected_samples_controls_30x = n_expected_samples_controls_30x,
                remote_outdir = remote_outdir,
                n_files_per_merge = n_files_per_merge,
                docker_image = docker_image,
                disk_size_gb = merge_disk_size_gb,
                ram_size_gb = merge_ram_size_gb,
                n_cpu = merge_n_cpu
        }
    }

    output {
        File merged_sample_ids_file = sample_ids
        Array[String] done = MergeChunk.done
    }
}


# Derives the canonical, deduplicated sample list by unioning `sample_ids.txt`
# (written by stage 1) from every source with a nonzero expected count. Only
# used when the caller doesn't already have this (e.g. from the orchestrator).
#
task GetSampleList {
    input {
        String remote_indir_bi
        String remote_indir_ha
        String remote_indir_bcm
        String remote_indir_uw
        String remote_indir_controls_15x
        String remote_indir_controls_30x

        Int n_expected_samples_bi
        Int n_expected_samples_ha
        Int n_expected_samples_bcm
        Int n_expected_samples_uw
        Int n_expected_samples_controls_15x
        Int n_expected_samples_controls_30x

        String docker_image
    }
    command <<<
        set -euxo pipefail
        rm -f all_sample_ids.txt
        touch all_sample_ids.txt

        if [ ~{n_expected_samples_bi} -gt 0 ]; then
            gcloud storage cp ~{remote_indir_bi}/sample_ids.txt bi_sample_ids.txt
            cat bi_sample_ids.txt >> all_sample_ids.txt
        fi
        if [ ~{n_expected_samples_ha} -gt 0 ]; then
            gcloud storage cp ~{remote_indir_ha}/sample_ids.txt ha_sample_ids.txt
            cat ha_sample_ids.txt >> all_sample_ids.txt
        fi
        if [ ~{n_expected_samples_bcm} -gt 0 ]; then
            gcloud storage cp ~{remote_indir_bcm}/sample_ids.txt bcm_sample_ids.txt
            cat bcm_sample_ids.txt >> all_sample_ids.txt
        fi
        if [ ~{n_expected_samples_uw} -gt 0 ]; then
            gcloud storage cp ~{remote_indir_uw}/sample_ids.txt uw_sample_ids.txt
            cat uw_sample_ids.txt >> all_sample_ids.txt
        fi
        if [ ~{n_expected_samples_controls_15x} -gt 0 ]; then
            gcloud storage cp ~{remote_indir_controls_15x}/sample_ids.txt controls_15x_sample_ids.txt
            cat controls_15x_sample_ids.txt >> all_sample_ids.txt
        fi
        if [ ~{n_expected_samples_controls_30x} -gt 0 ]; then
            gcloud storage cp ~{remote_indir_controls_30x}/sample_ids.txt controls_30x_sample_ids.txt
            cat controls_30x_sample_ids.txt >> all_sample_ids.txt
        fi

        sort -u all_sample_ids.txt > sample_ids.txt
        wc -l sample_ids.txt 1>&2
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

        String remote_indir_bi
        String remote_indir_ha
        String remote_indir_bcm
        String remote_indir_uw
        String remote_indir_controls_15x
        String remote_indir_controls_30x

        Array[String] bi_samples_to_prefer_over_ha

        Int n_expected_samples_bi
        Int n_expected_samples_ha
        Int n_expected_samples_bcm
        Int n_expected_samples_uw
        Int n_expected_samples_controls_15x
        Int n_expected_samples_controls_30x

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

        N_SAMPLES=$(cat ~{sample_ids} | wc -l)

        # ------------------- List + validate every enabled source --------------
        date 1>&2
        rm -f all_remote_files.txt
        touch all_remote_files.txt

        if [ ~{n_expected_samples_bi} -gt 0 ]; then
            gcloud storage ls -l ~{remote_indir_bi}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > bi_files.txt
            N_FILES=$(wc -l < bi_files.txt)
            N_FILES=$(( ${N_FILES} - 1 ))
            if [ ${N_FILES} -ne ~{n_expected_samples_bi} ]; then
                echo "ERROR: bi has ${N_FILES} files != ~{n_expected_samples_bi}"
                exit 1
            fi
            head -n ${N_FILES} bi_files.txt >> all_remote_files.txt
        fi

        if [ ~{n_expected_samples_ha} -gt 0 ]; then
            gcloud storage ls -l ~{remote_indir_ha}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > ha_files.txt
            N_FILES=$(wc -l < ha_files.txt)
            N_FILES=$(( ${N_FILES} - 1 ))
            if [ ${N_FILES} -ne ~{n_expected_samples_ha} ]; then
                echo "ERROR: ha has ${N_FILES} files != ~{n_expected_samples_ha}"
                exit 1
            fi
            head -n ${N_FILES} ha_files.txt >> all_remote_files.txt
        fi

        if [ ~{n_expected_samples_bcm} -gt 0 ]; then
            gcloud storage ls -l ~{remote_indir_bcm}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > bcm_files.txt
            N_FILES=$(wc -l < bcm_files.txt)
            N_FILES=$(( ${N_FILES} - 1 ))
            if [ ${N_FILES} -ne ~{n_expected_samples_bcm} ]; then
                echo "ERROR: bcm has ${N_FILES} files != ~{n_expected_samples_bcm}"
                exit 1
            fi
            head -n ${N_FILES} bcm_files.txt >> all_remote_files.txt
        fi

        if [ ~{n_expected_samples_uw} -gt 0 ]; then
            gcloud storage ls -l ~{remote_indir_uw}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > uw_files.txt
            N_FILES=$(wc -l < uw_files.txt)
            N_FILES=$(( ${N_FILES} - 1 ))
            if [ ${N_FILES} -ne ~{n_expected_samples_uw} ]; then
                echo "ERROR: uw has ${N_FILES} files != ~{n_expected_samples_uw}"
                exit 1
            fi
            head -n ${N_FILES} uw_files.txt >> all_remote_files.txt
        fi

        if [ ~{n_expected_samples_controls_15x} -gt 0 ]; then
            gcloud storage ls -l ~{remote_indir_controls_15x}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > control_15x_files.txt
            N_FILES=$(wc -l < control_15x_files.txt)
            N_FILES=$(( ${N_FILES} - 1 ))
            if [ ${N_FILES} -lt ~{n_expected_samples_controls_15x} ]; then
                echo "ERROR: controls_15x has ${N_FILES} files < ~{n_expected_samples_controls_15x}"
                exit 1
            elif [ ${N_FILES} -gt ~{n_expected_samples_controls_15x} ]; then
                echo "WARNING: controls_15x has ${N_FILES} files > ~{n_expected_samples_controls_15x}"
            fi
            head -n ${N_FILES} control_15x_files.txt >> all_remote_files.txt
        fi

        if [ ~{n_expected_samples_controls_30x} -gt 0 ]; then
            gcloud storage ls -l ~{remote_indir_controls_30x}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > control_30x_files.txt
            N_FILES=$(wc -l < control_30x_files.txt)
            N_FILES=$(( ${N_FILES} - 1 ))
            if [ ${N_FILES} -lt ~{n_expected_samples_controls_30x} ]; then
                echo "ERROR: controls_30x has ${N_FILES} files < ~{n_expected_samples_controls_30x}"
                exit 1
            elif [ ${N_FILES} -gt ~{n_expected_samples_controls_30x} ]; then
                echo "WARNING: controls_30x has ${N_FILES} files > ~{n_expected_samples_controls_30x}"
            fi
            head -n ${N_FILES} control_30x_files.txt >> all_remote_files.txt
        fi
        date 1>&2

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
        rm -f *_files.txt

        # ------------------- Localize, resolving cross-source overlap ----------
        # Sources are localized into the SAME flat ./input_files/ dir, keyed by
        # `<sample_id>.bcf`. A later source silently overwrites an earlier
        # source's file for any sample_id present in both. Default order below
        # (bi, ha, uw, bcm, controls) means ha wins over bi, and uw/bcm/controls
        # win over everything before them, on overlap. bi_samples_to_prefer_over_ha
        # is the one explicit override: for exactly the listed sample_ids, bi's
        # copy is re-localized after ha, flipping the winner back to bi. There is
        # no equivalent override for uw/bcm/controls.
        date 1>&2
        mkdir ./input_files/
        if [ ~{n_expected_samples_bi} -gt 0 ]; then
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir_bi}/chunk_~{chunk_id}/'*' ./input_files/
        fi
        if [ ~{n_expected_samples_ha} -gt 0 ]; then
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir_ha}/chunk_~{chunk_id}/'*' ./input_files/
        fi
        if [ ~{n_expected_samples_bi} -gt 0 -a ~{n_expected_samples_ha} -gt 0 ]; then
            echo ~{sep="," bi_samples_to_prefer_over_ha} | tr ',' '\n' > bi_samples_to_prefer_over_ha.txt
            rm -f bi_override_files.txt
            touch bi_override_files.txt
            while read -u 5 SAMPLE_ID; do
                [ -z "${SAMPLE_ID}" ] && continue
                echo "~{remote_indir_bi}/chunk_~{chunk_id}/${SAMPLE_ID}.bcf" >> bi_override_files.txt
                echo "~{remote_indir_bi}/chunk_~{chunk_id}/${SAMPLE_ID}.bcf.csi" >> bi_override_files.txt
            done 5< bi_samples_to_prefer_over_ha.txt
            if [ -s bi_override_files.txt ]; then
                cat bi_override_files.txt | gcloud storage cp -I ./input_files/
            fi
            rm -f bi_samples_to_prefer_over_ha.txt bi_override_files.txt
        fi
        if [ ~{n_expected_samples_uw} -gt 0 ]; then
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir_uw}/chunk_~{chunk_id}/'*' ./input_files/
        fi
        if [ ~{n_expected_samples_bcm} -gt 0 ]; then
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir_bcm}/chunk_~{chunk_id}/'*' ./input_files/
        fi
        if [ ~{n_expected_samples_controls_15x} -gt 0 ]; then
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir_controls_15x}/chunk_~{chunk_id}/'*' ./input_files/
        fi
        if [ ~{n_expected_samples_controls_30x} -gt 0 ]; then
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir_controls_30x}/chunk_~{chunk_id}/'*' ./input_files/
        fi
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
