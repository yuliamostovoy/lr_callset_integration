version 1.0

import "SV_Integration_Workpackage12.wdl" as wp12
import "SV_Integration_Workpackage13.wdl" as wp13
import "SV_Integration_Workpackage14.wdl" as wp14
import "SV_Integration_Workpackage15.wdl" as wp15


# One suffix (`ultralong` or `bnd`) of the non-scored cohort integration:
# cohort merge (WP12) -> shard (WP13) -> truvari collapse (WP14) -> concat (WP15),
# chained in one submission with no data table and no hand-built inter-step
# files. Called once per suffix by SV_Integration_WorkflowD_Ultralong_Bnd.
#
workflow SV_Integration_WorkflowD_PerSuffix {
    input {
        String suffix
        String remote_indir
        String remote_outdir_suffix

        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        Int? n_expected_samples

        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        Int consistency_checks = 1

        File reference_fa
        File reference_fai
        String truvari_matching_parameters = "--refdist 500 --pctseq 0.95 --pctsize 0.95 --pctovl 0.0"
        Int max_resolve = 100000
        Boolean use_bed = false
        Int chunk_ids_per_file = 100

        Int concat_all_naive = 1

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        suffix: "'ultralong' or 'bnd'."
        remote_indir: "Workflow A intra-sample output dir holding <sample>_<suffix>.bcf per sample."
        remote_outdir_suffix: "Stage outputs go to /12_merge, /13_shard, /14_collapse, /15_concat under here; the genome-wide callset is /15_concat/truvari_collapsed.bcf."
        n_expected_samples: "OPTIONAL. Auto-derived from the remote_indir listing when omitted."
    }

    String chromosomes_str = "~{sep=',' chromosomes}"
    String merge_dir = remote_outdir_suffix + "/12_merge"
    String shard_dir = remote_outdir_suffix + "/13_shard"
    String collapse_dir = remote_outdir_suffix + "/14_collapse"
    String concat_dir = remote_outdir_suffix + "/15_concat"

    call WriteSampleListSuffix {
        input:
            remote_indir = remote_indir,
            suffix = suffix,
            docker_image = docker_image
    }
    Int n_samples = select_first([n_expected_samples, WriteSampleListSuffix.n_samples])

    call wp12.Impl as Merge {
        input:
            sample_ids = WriteSampleListSuffix.sample_ids_file,
            suffix = suffix,
            chromosomes = chromosomes_str,
            remote_indir = remote_indir,
            n_expected_samples = n_samples,
            remote_outdir = merge_dir,
            docker_image = docker_image
    }

    scatter (chr in chromosomes) {
        call wp13.Impl as Shard {
            input:
                chromosome_id = chr,
                suffix = suffix,
                truvari_chunk_min_records = truvari_chunk_min_records,
                truvari_collapse_refdist = truvari_collapse_refdist,
                consistency_checks = consistency_checks,
                remote_indir = merge_dir,
                remote_outdir = shard_dir,
                docker_image = docker_image,
                upstream_signal = [Merge.done]
        }
        call DeriveChunkIds {
            input:
                regions_txt = Shard.regions_txt,
                chunk_ids_per_file = chunk_ids_per_file,
                docker_image = docker_image
        }
        scatter (idfile in DeriveChunkIds.chunk_id_files) {
            call wp14.Impl as Collapse {
                input:
                    remote_indir = shard_dir,
                    chromosome_id = chr,
                    chunks_ids = idfile,
                    remote_outdir = collapse_dir,
                    reference_fa = reference_fa,
                    reference_fai = reference_fai,
                    truvari_matching_parameters = truvari_matching_parameters,
                    max_resolve = max_resolve,
                    use_bed = use_bed,
                    docker_image = docker_image
            }
        }
        call wp15.SingleChromosome {
            input:
                chromosome = chr,
                remote_indir = collapse_dir,
                remote_outdir = concat_dir,
                docker_image = docker_image,
                upstream_signal = Collapse.done
        }
    }

    call wp15.AllChromosomes {
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


# Cohort sample list + count from a listing of the per-sample <sample>_<suffix>.bcf
# files in remote_indir. Replaces the hand-built sample_ids file and the manual
# n_expected_samples value.
#
task WriteSampleListSuffix {
    input {
        String remote_indir
        String suffix
        String docker_image
    }

    command <<<
        set -euxo pipefail
        gcloud storage ls ~{remote_indir}/'*_~{suffix}.bcf' | sed 's#.*/##; s#_~{suffix}\.bcf$##' | sort -u > sample_ids.txt
        if [ ! -s sample_ids.txt ]; then
            echo "ERROR: no <sample>_~{suffix}.bcf files found under ~{remote_indir}."
            exit 1
        fi
        wc -l < sample_ids.txt > n.txt
        cat sample_ids.txt 1>&2
    >>>

    output {
        File sample_ids_file = "sample_ids.txt"
        Int n_samples = read_int("n.txt")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}


# Per-chromosome truvari chunk-id lists from WP13's regions.txt (col 2 = chunk
# id), split into files of `chunk_ids_per_file` ids. Replaces
# make_workpackage7_chunk_id_files.sh; the resulting array drives the WP14 scatter.
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
