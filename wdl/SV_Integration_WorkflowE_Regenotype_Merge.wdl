version 1.0

import "SV_Integration_Workpackage3_Main_bcftools_merge.wdl" as wp3
import "SV_Integration_Workpackage8_Main_concat_regenotyped_shards.wdl" as wp8


# CONSOLIDATED STEP E of the SV integration pipeline: the cross-family GATHER
# that turns per-family regenotyped chunks into ONE cohort VCF.
#
# The regenotyping workflows (Workflow C for the main/kanpig branch, and the
# cuteFC test workflow for the ultralong branch) scatter one family per Terra
# instance and write per-sample chunks (chunk_<i>/<sample>.bcf + <sample>.done)
# into a shared dir. No single instance sees all samples, so the cohort-wide
# merge happens here, once, after they all finish:
#   WriteSampleList (all samples across all families, from the .done markers)
#   -> WP3 merge_mode=2 (bcftools merge by ID) per chunk
#   -> WP8 concat -> <remote_outdir>/concat/merged.bcf
#
# Merge-by-ID is variant-type-agnostic, so ONE definition serves both branches:
# run it once per branch, pointed at that branch's regenotyping dir. Only the
# type-specific truvari-collapse parameters differ between branches, and those
# live in the cohort BUILD (Workflows B/D), not here.
#
workflow SV_Integration_WorkflowE_Regenotype_Merge {
    input {
        File split_for_bcftools_merge_csv

        # remote_indir = the SHARED regenotyping dir the per-family instances wrote
        # to (Workflow C's remote_outdir, or the cuteFC test's regenotyping dir).
        String remote_indir
        String remote_outdir

        Int merge_mode = 2
        File? sample_ids_file

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"

        # Runtime knobs (ultralong records are bulkier -> may want more).
        Int merge_n_cpu = 4
        Int merge_ram_size_gb = 16
        Int merge_disk_size_gb = 100
    }
    parameter_meta {
        split_for_bcftools_merge_csv: "The interval partition CSV the regenotyping used; chunk id == 0-based line number. Defines the chunk scatter width."
        remote_indir: "Shared dir the per-family regenotyping instances wrote to: chunk_<i>/<sample>.bcf + <sample>.done for every sample across every family."
        remote_outdir: "Merged chunks go to /merge; the final cohort callset is /concat/merged.bcf."
        merge_mode: "bcftools merge mode. 2 (default) = merge by ID, correct for regenotyped calls sharing cohort IDs. Both branches use 2."
        sample_ids_file: "OPTIONAL. If omitted, the cohort sample list (merge column order) is auto-derived from the <sample>.done markers in remote_indir."
    }

    String indir = sub(remote_indir, "/+$", "")
    String outdir = sub(remote_outdir, "/+$", "")
    String merge_dir = outdir + "/merge"
    String concat_dir = outdir + "/concat"

    Int n_chunks = length(read_lines(split_for_bcftools_merge_csv))

    # Cohort sample list = every sample across every family that was regenotyped.
    if (!defined(sample_ids_file)) {
        call WriteSampleList {
            input:
                remote_indir = indir,
                docker_image = docker_image
        }
    }
    File sample_ids = select_first([sample_ids_file, WriteSampleList.sample_ids_file])

    call MakeChunkIdsCsv {
        input:
            n_chunks = n_chunks,
            docker_image = docker_image
    }

    scatter (chunk_id in range(n_chunks)) {
        call wp3.Impl as MergeChunk {
            input:
                chunk_id = chunk_id,
                sample_ids = sample_ids,
                remote_indir = indir,
                merge_mode = merge_mode,
                remote_outdir = merge_dir,
                docker_image = docker_image,
                n_cpu = merge_n_cpu,
                ram_size_gb = merge_ram_size_gb,
                disk_size_gb = merge_disk_size_gb
        }
    }

    call wp8.Impl as Concat {
        input:
            chunk_ids = MakeChunkIdsCsv.csv,
            remote_indir = merge_dir,
            remote_outdir = concat_dir,
            docker_image = docker_image,
            upstream_signal = MergeChunk.done
    }

    output {
    }
}


# Cohort sample list (bcftools merge column order) from the <sample>.done markers
# the regenotyping instances left in the shared dir.
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


# Comma-separated 0..n_chunks-1 for the final WP8 concat.
#
task MakeChunkIdsCsv {
    input {
        Int n_chunks
        String docker_image
    }

    command <<<
        set -euxo pipefail
        seq 0 $(( ~{n_chunks} - 1 )) | paste -sd, - > csv.txt
    >>>

    output {
        String csv = read_string("csv.txt")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}
