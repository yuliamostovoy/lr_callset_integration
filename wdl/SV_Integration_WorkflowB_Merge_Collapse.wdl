version 1.0

import "SV_Integration_Workpackage3_Main_bcftools_merge.wdl" as wp3
import "SV_Integration_Workpackage4_Main_shard.wdl" as wp4
import "SV_Integration_Workpackage5_Main_truvari_collapse.wdl" as wp5
import "SV_Integration_Workpackage6_Main_concat_shards.wdl" as wp6


# CONSOLIDATED STEP B of the SV integration pipeline: cohort bcftools merge
# (WP3) -> shard into truvari-collapse chunks (WP4) -> truvari collapse (WP5) ->
# concat into the genome-wide callset (WP6), in a SINGLE Terra submission with
# NO data table and NO hand-built inter-step files.
#
# The cohort sample list, the per-chromosome `bcftools_chunks` strings, and the
# per-chromosome truvari chunk-id lists are all derived in-graph (from the same
# interval CSV plus a listing of Workflow A's output dir), replacing the manual
# sample_ids file, the comma-string chunks, and make_workpackage7_chunk_id_files.sh.
#
workflow SV_Integration_WorkflowB_Merge_Collapse {
    input {
        File split_for_bcftools_merge_csv

        String remote_indir
        String remote_outdir

        Int merge_mode = 1
        File? sample_ids_file

        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        Int consistency_checks = 1

        String truvari_matching_parameters = "--refdist 500 --pctseq 0.95 --pctsize 0.95 --pctovl 0.0"
        Boolean use_bed = false
        Int chunk_ids_per_file = 100

        Int concat_all_naive = 1

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        split_for_bcftools_merge_csv: "The interval partition CSV; chunk id == 0-based line number. Same file Workflow A used."
        remote_indir: "Workflow A's remote_outdir (this workflow reads its /02_scoring subdir automatically)."
        remote_outdir: "Stage outputs go to /03_merge, /04_shard, /05_collapse, /06_concat under here; the genome-wide callset is /06_concat/truvari_collapsed.bcf."
        merge_mode: "1: bcftools merge by CHROM,POS,REF,ALT (default for the main chain). 2: merge by ID."
        sample_ids_file: "OPTIONAL. If omitted, the cohort sample list (and thus the bcftools-merge column order) is auto-derived from the <sample>.done markers in remote_indir."
        chunk_ids_per_file: "Truvari chunk ids per WP5 VM (replaces make_workpackage7_chunk_id_files.sh split size)."
    }

    # Workflow A writes its scored per-sample chunks to the fixed /02_scoring
    # subdir of its remote_outdir, so the user passes A's remote_outdir here.
    String indir = sub(remote_indir, "/+$", "") + "/02_scoring"
    String outdir = sub(remote_outdir, "/+$", "")
    String merge_dir = outdir + "/03_merge"
    String shard_dir = outdir + "/04_shard"
    String collapse_dir = outdir + "/05_collapse"
    String concat_dir = outdir + "/06_concat"

    Int n_chunks = length(read_lines(split_for_bcftools_merge_csv))

    # Cohort sample list (bcftools merge column order), auto-derived unless given.
    if (!defined(sample_ids_file)) {
        call WriteSampleList {
            input:
                remote_indir = indir,
                docker_image = docker_image
        }
    }
    File sample_ids = select_first([sample_ids_file, WriteSampleList.sample_ids_file])

    # --- WP3: cohort bcftools merge, one VM per chunk ---
    scatter (chunk_id in range(n_chunks)) {
        call wp3.Impl as Merge {
            input:
                chunk_id = chunk_id,
                sample_ids = sample_ids,
                remote_indir = indir,
                merge_mode = merge_mode,
                remote_outdir = merge_dir,
                docker_image = docker_image
        }
    }

    # --- WP4 -> DeriveChunkIds -> WP5, per chromosome ---
    scatter (chr in chromosomes) {
        call ChunksForChromosome {
            input:
                split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
                chromosome_id = chr,
                docker_image = docker_image
        }
        call wp4.Impl as Shard {
            input:
                chromosome_id = chr,
                bcftools_chunks = ChunksForChromosome.chunks,
                truvari_chunk_min_records = truvari_chunk_min_records,
                truvari_collapse_refdist = truvari_collapse_refdist,
                consistency_checks = consistency_checks,
                remote_indir = merge_dir,
                remote_outdir = shard_dir,
                docker_image = docker_image,
                upstream_signal = Merge.done
        }
        call DeriveChunkIds {
            input:
                regions_txt = Shard.regions_txt,
                chunk_ids_per_file = chunk_ids_per_file,
                docker_image = docker_image
        }
        scatter (idfile in DeriveChunkIds.chunk_id_files) {
            call wp5.Impl as Collapse {
                input:
                    remote_indir = shard_dir,
                    chromosome_id = chr,
                    chunks_ids = idfile,
                    remote_outdir = collapse_dir,
                    truvari_matching_parameters = truvari_matching_parameters,
                    use_bed = use_bed,
                    docker_image = docker_image
            }
        }
        call wp6.SingleChromosome {
            input:
                chromosome = chr,
                remote_indir = collapse_dir,
                remote_outdir = concat_dir,
                docker_image = docker_image,
                upstream_signal = Collapse.done
        }
    }

    # --- WP6: genome-wide concat ---
    call wp6.AllChromosomes {
        input:
            chromosomes = chromosomes,
            out_txt = SingleChromosome.out_txt,
            remote_outdir = concat_dir,
            naive = concat_all_naive,
            docker_image = docker_image
    }

    output {
    }
}


# Cohort sample list from a listing of Workflow A's <sample>.done markers. This
# is both the merge column order and the per-chunk file-count expectation.
#
task WriteSampleList {
    input {
        String remote_indir
        String docker_image
    }

    command <<<
        set -euxo pipefail
        gcloud storage ls ~{remote_indir}/'*.done' | sed 's#.*/##; s#\.done$##' | sort -u > sample_ids.txt
        wc -l sample_ids.txt 1>&2
        if [ ! -s sample_ids.txt ]; then
            echo "ERROR: no <sample>.done markers found under ~{remote_indir}."
            exit 1
        fi
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


# Comma-separated, POS-sorted bcftools-chunk ids for one chromosome, taken from
# the interval CSV (chunk id == 0-based line number). Replaces the hand-typed
# `bcftools_chunks` string.
#
task ChunksForChromosome {
    input {
        File split_for_bcftools_merge_csv
        String chromosome_id
        String docker_image
    }

    command <<<
        set -euxo pipefail
        awk -F, -v chr=~{chromosome_id} '$1==chr {print NR-1}' ~{split_for_bcftools_merge_csv} | paste -sd, - > chunks.txt
        if [ ! -s chunks.txt ]; then
            echo "ERROR: chromosome ~{chromosome_id} not present in the interval CSV."
            exit 1
        fi
        cat chunks.txt 1>&2
    >>>

    output {
        String chunks = read_string("chunks.txt")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}


# Per-chromosome truvari chunk-id lists from WP4's regions.txt (col 2 = chunk
# id), split into files of `chunk_ids_per_file` ids. Replaces
# make_workpackage7_chunk_id_files.sh; the resulting array drives the WP5 scatter.
#
task DeriveChunkIds {
    input {
        File regions_txt
        Int chunk_ids_per_file
        String docker_image
    }

    command <<<
        set -euxo pipefail
        awk 'NF>=2 {print $2}' ~{regions_txt} > ids.txt
        if [ ! -s ids.txt ]; then
            echo "ERROR: no chunk ids found in regions.txt."
            exit 1
        fi
        split -l ~{chunk_ids_per_file} -d -a 4 ids.txt chunk_ids_
        ls chunk_ids_* 1>&2
    >>>

    output {
        Array[File] chunk_id_files = glob("chunk_ids_*")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}
