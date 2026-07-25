version 1.0

import "SV_Integration_Workpackage7_Main_joint_genotype_families.wdl" as wp7


# CONSOLIDATED STEP C of the SV integration pipeline: family joint regenotyping
# (WP7 kanpig) ONLY. Driven by a Terra SAMPLE SET (one sample_set == one family;
# Terra runs N parallel copies over the selected sets). Each instance force-calls
# its family's members and writes per-sample chunks to a SHARED `remote_outdir`.
#
# The cross-family cohort merge is deliberately NOT here: because the run scatters
# one family per instance, no single instance sees all samples. Run Workflow E
# (SV_Integration_WorkflowE_Regenotype_Merge) once afterward, pointed at this same
# `remote_outdir`, to merge all families' per-sample chunks by ID and concat them
# into the single cohort VCF. (Per-sample chunk filenames are unique across
# families, so concurrent instances sharing one `remote_outdir` do not collide.)
#
# The family's members are bound from `this.samples.*`. `MakeFamilyInputs`
# validates that the set's members EXACTLY match the family defined in the
# provided PED (catches a mis-built Terra set); the PED is passed through to WP7,
# which derives family membership from it, as in the standalone runs.
#
workflow SV_Integration_WorkflowC_Regenotype {
    input {
        # --- One sample_set per instance; members from the set ---
        String family_id
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bais
        Array[String] aligned_bams

        File ped
        File split_for_bcftools_merge_csv

        # --- GCS dirs (no final slash) ---
        String remote_indir
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
        ped: "Standard 6-column PED. The family named by `family_id` (column 1) must list exactly this set's members in column 2 — the workflow errors otherwise. Passed through to WP7, which derives family membership from it. Columns 3-6 are unused by WP7_families."
        split_for_bcftools_merge_csv: "The interval partition CSV; chunk id == 0-based line number. Same file the rest of the pipeline used."
        remote_indir: "Workflow B's remote_outdir (this workflow reads its /06_concat subdir, holding the genome-wide truvari_collapsed.bcf, automatically)."
        remote_outdir: "SHARED across all family-sets of the cohort. Per-sample regenotyped chunks (chunk_<i>/<sample>.bcf) + <sample>.done land directly here. Feed this same dir to Workflow E to produce the cohort VCF."
    }

    # Workflow B writes the genome-wide cohort callset to the fixed /06_concat
    # subdir of its remote_outdir, so the user passes B's remote_outdir here.
    String cohort_indir = sub(remote_indir, "/+$", "") + "/06_concat"
    String regeno_dir = sub(remote_outdir, "/+$", "")

    # Validate the set against the PED before spending compute on regenotyping.
    call MakeFamilyInputs {
        input:
            family_id = family_id,
            sample_ids = sample_ids,
            ped = ped,
            docker_image = docker_image
    }

    call wp7.Impl as JointGenotype {
        input:
            family_ids = [family_id],
            ped = ped,
            sample_ids = sample_ids,
            sample_sexes = sample_sexes,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = cohort_indir,
            remote_outdir = regeno_dir,
            requester_pays_project = requester_pays_project,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            autosomes_bed = autosomes_bed,
            kanpig_params_cohort = kanpig_params_cohort,
            docker_image = docker_image,
            upstream_signal = [MakeFamilyInputs.done]
    }

    output {
    }
}


# Validates that the Terra set's members EXACTLY match the family defined in the
# PED (col 1 == family_id, col 2 == sample_id -- the same extraction WP7 uses).
# Fails the run if they disagree (mis-built set) or if `family_id` is absent from
# PED column 1. Emits a `done` signal that gates WP7 so validation fails fast,
# before any regenotyping compute.
#
task MakeFamilyInputs {
    input {
        String family_id
        Array[String] sample_ids
        File ped
        String docker_image
    }

    command <<<
        set -euxo pipefail

        # Members from the Terra sample_set.
        cat > set_members.raw <<'EOF_SAMPLE_IDS'
~{sep="\n" sample_ids}
EOF_SAMPLE_IDS
        grep -v '^[[:space:]]*$' set_members.raw | sort -u > set_members.txt
        if [ ! -s set_members.txt ]; then
            echo "ERROR: sample set ~{family_id} has no members."
            exit 1
        fi

        # Members of this family per the PED.
        awk -v fam="~{family_id}" 'BEGIN { FS="[ \t]+" } $1==fam && $2!="0" && $2!="." { print $2 }' ~{ped} | sort -u > ped_members.txt
        if [ ! -s ped_members.txt ]; then
            echo "ERROR: family_id '~{family_id}' not found in PED column 1 (or it lists no members). The Terra sample_set id must match the PED family id."
            exit 1
        fi

        # The set must EXACTLY match the PED family -- guards against a mis-built set.
        if ! diff -q set_members.txt ped_members.txt >/dev/null; then
            echo "ERROR: sample set '~{family_id}' does not match its PED family." 1>&2
            echo "  In the set but NOT in the PED family:" 1>&2; comm -23 set_members.txt ped_members.txt 1>&2
            echo "  In the PED family but NOT in the set:" 1>&2; comm -13 set_members.txt ped_members.txt 1>&2
            exit 1
        fi

        echo "validated" > validated.txt
    >>>

    output {
        String done = read_string("validated.txt")
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "2GB"
        disks: "local-disk 16 HDD"
        preemptible: 3
    }
}
