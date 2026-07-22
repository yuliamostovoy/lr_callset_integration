version 1.0

import "SV_Integration_Workpackage7_Main_joint_genotype_families.wdl" as wp7
import "SV_Integration_Workpackage3_Main_bcftools_merge.wdl" as wp3
import "SV_Integration_Workpackage8_Main_concat_regenotyped_shards.wdl" as wp8


# CONSOLIDATED STEP C of the SV integration pipeline: family joint regenotyping
# (WP7) -> per-chunk cross-member merge by ID (WP3, merge_mode=2) -> genome-wide
# concat of the merged chunks (WP8), in a SINGLE Terra submission driven by a
# Terra SAMPLE SET (one family per workflow instance; Terra runs N parallel
# copies over the selected sets).
#
# The family's members are bound from `this.samples.*`; `MakeFamilyInputs`
# synthesizes the PED and the merge sample-order file, and derives the
# genome-ordered chunk-id list from the interval CSV. This replaces the
# hand-aligned parallel arrays, make_terra_family_set_membership.sh, and the
# hand-typed chunk_ids string.
#
# The WP3 per-chunk merge between WP7 and WP8 is REQUIRED: WP7 writes per-sample
# chunk_<i>/<sample>.bcf, but WP8 reads a single merged chunk_<i>.bcf.
#
workflow SV_Integration_WorkflowC_Regenotype {
    input {
        # --- One sample_set per instance; members from the set ---
        String family_id
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bais
        Array[String] aligned_bams

        File split_for_bcftools_merge_csv

        # --- GCS dirs (no final slash) ---
        String wp6_remote_indir
        String remote_outdir

        String requester_pays_project = ""
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed
        String kanpig_params_cohort = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        family_id: "e.g. `this.sample_set_id`. If Terra does not expose the set id, pass a literal per submission."
        sample_ids: "e.g. `this.samples.sample_id`. Members of this set. Parallel to the other member arrays."
        sample_sexes: "e.g. `this.samples.sample_sex`."
        aligned_bais: "e.g. `this.samples.aligned_bai` (gs:// URI)."
        aligned_bams: "e.g. `this.samples.aligned_bam` (gs:// URI)."
        split_for_bcftools_merge_csv: "The interval partition CSV; chunk id == 0-based line number. Same file the rest of the pipeline used."
        wp6_remote_indir: "Workflow B's /06_concat dir holding the genome-wide truvari_collapsed.bcf(.csi)."
        remote_outdir: "Outputs go to /07_regenotype (per-sample chunks), /07b_merge (per-chunk merged), /08_concat (final merged.bcf) under here."
    }

    String wp6_indir = sub(wp6_remote_indir, "/+$", "")
    String outdir = sub(remote_outdir, "/+$", "")
    String regeno_dir = outdir + "/07_regenotype"
    String merge_dir = outdir + "/07b_merge"
    String concat_dir = outdir + "/08_concat"

    Int n_chunks = length(read_lines(split_for_bcftools_merge_csv))

    call MakeFamilyInputs {
        input:
            family_id = family_id,
            sample_ids = sample_ids,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            docker_image = docker_image
    }

    call wp7.Impl as JointGenotype {
        input:
            family_ids = [family_id],
            ped = MakeFamilyInputs.ped,
            sample_ids = sample_ids,
            sample_sexes = sample_sexes,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = wp6_indir,
            remote_outdir = regeno_dir,
            requester_pays_project = requester_pays_project,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            autosomes_bed = autosomes_bed,
            kanpig_params_cohort = kanpig_params_cohort,
            docker_image = docker_image
    }

    scatter (chunk_id in range(n_chunks)) {
        call wp3.Impl as MergeChunk {
            input:
                chunk_id = chunk_id,
                sample_ids = MakeFamilyInputs.sample_ids_file,
                remote_indir = regeno_dir,
                merge_mode = 2,
                remote_outdir = merge_dir,
                docker_image = docker_image,
                upstream_signal = [JointGenotype.done]
        }
    }

    call wp8.Impl as Concat {
        input:
            chunk_ids = MakeFamilyInputs.chunk_ids_csv,
            remote_indir = merge_dir,
            remote_outdir = concat_dir,
            docker_image = docker_image,
            upstream_signal = MergeChunk.done
    }

    output {
    }
}


# Synthesizes, for one sample set:
#   - a 6-column PED (family_id, sample_id, 0, 0, 0, 0) that WP7 reads (it only
#     uses columns 1-2), replacing make_terra_family_set_membership.sh;
#   - the sample-order file WP3 needs as its `File sample_ids` merge order;
#   - the genome-ordered comma-separated chunk-id string WP8 needs, derived from
#     the interval CSV (chunk id == 0-based line number).
#
task MakeFamilyInputs {
    input {
        String family_id
        Array[String] sample_ids
        File split_for_bcftools_merge_csv
        String docker_image
    }

    command <<<
        set -euxo pipefail

        cat > sample_ids.txt <<'EOF_SAMPLE_IDS'
~{sep="\n" sample_ids}
EOF_SAMPLE_IDS
        grep -v '^[[:space:]]*$' sample_ids.txt | sort -u > sample_ids.clean.txt
        mv sample_ids.clean.txt sample_ids.txt
        if [ ! -s sample_ids.txt ]; then
            echo "ERROR: sample set ~{family_id} has no members."
            exit 1
        fi

        awk -v fam="~{family_id}" 'BEGIN { OFS="\t" } { print fam, $1, 0, 0, 0, 0 }' sample_ids.txt > family.ped

        N=$(wc -l < ~{split_for_bcftools_merge_csv})
        seq 0 $(( ${N} - 1 )) | paste -sd, - > chunk_ids_csv.txt
    >>>

    output {
        File ped = "family.ped"
        File sample_ids_file = "sample_ids.txt"
        String chunk_ids_csv = read_string("chunk_ids_csv.txt")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}
