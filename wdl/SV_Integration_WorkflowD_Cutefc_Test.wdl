version 1.0

import "SV_Integration_Workpackage7_Main_joint_genotype_families_cutefc.wdl" as wp7cutefc


# TEST-ONLY workflow (not production): per-family cuteFC regenotyping of the
# cohort ULTRALONG callset, for testing a cuteFC fork.
#
# Mirrors Workflow C's structure exactly, but with the cuteFC WP7 variant and
# reading the ultralong cohort VCF:
#   - Root the submission on a Terra SAMPLE SET (one set == one family); Terra
#     scatters one instance per set.
#   - Validate the set's members against the PED, then force-call the family's
#     members with cuteFC against the family-present sites of the cohort ultralong
#     callset, writing per-sample chunks to a SHARED remote_outdir.
#
# The cohort ultralong callset is built ONCE by a separate Workflow D (ultralong)
# submission and read here as a prebuilt input -- exactly as Workflow C reads
# Workflow B's cohort VCF. This workflow does NOT build it (WP12-15), and does
# NOT merge: run Workflow E once afterward on this same remote_outdir to produce
# the single cohort VCF.
#
workflow SV_Integration_WorkflowD_Cutefc_Test {
    input {
        # --- One sample_set per instance; members from the set ---
        String family_id
        Array[String] sample_ids
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
        File autosomes_bed
        String cutefc_params_cohort = "--max_size -1 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"
        String cutefc_docker_image = "quay.io/ymostovoy/lr-ultralong:latest"
        String wp_docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        family_id: "e.g. `this.sample_set_id`. Must match a family id in PED column 1."
        sample_ids: "e.g. `this.samples.sample_id`. Members of this set; must exactly match the PED family."
        aligned_bais: "e.g. `this.samples.aligned_bai` (gs:// URI)."
        aligned_bams: "e.g. `this.samples.aligned_bam` (gs:// URI)."
        ped: "6-column PED. The family named by `family_id` (col 1) must list exactly this set's members (col 2); the workflow errors otherwise. Passed through to the cuteFC WP7 variant."
        split_for_bcftools_merge_csv: "The interval partition CSV; chunk id == 0-based line number."
        remote_indir: "Workflow D's remote_outdir (this reads its /ultralong/15_concat subdir, holding the cohort ultralong truvari_collapsed.bcf, automatically)."
        remote_outdir: "SHARED across all family-sets. Per-sample cuteFC chunks + <sample>.done land directly here. Feed this dir to Workflow E for the cohort VCF."
        cutefc_docker_image: "The forked-cuteFC ultralong image under test."
    }

    # Workflow D writes the cohort ultralong callset to /ultralong/15_concat.
    String cohort_indir = sub(remote_indir, "/+$", "") + "/ultralong/15_concat"
    String regeno_dir = sub(remote_outdir, "/+$", "")

    # Validate the set against the PED before spending compute.
    call ValidateFamilySet {
        input:
            family_id = family_id,
            sample_ids = sample_ids,
            ped = ped,
            docker_image = wp_docker_image
    }

    call wp7cutefc.Impl as Regenotype {
        input:
            family_ids = [family_id],
            ped = ped,
            sample_ids = sample_ids,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = cohort_indir,
            remote_outdir = regeno_dir,
            requester_pays_project = requester_pays_project,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            autosomes_bed = autosomes_bed,
            cutefc_params_cohort = cutefc_params_cohort,
            docker_image = cutefc_docker_image,
            upstream_signal = [ValidateFamilySet.done]
    }

    output {
    }
}


# Validates that the Terra set's members EXACTLY match the family defined in the
# PED (col 1 == family_id, col 2 == sample_id). Fails fast, gating the cuteFC
# stage. Identical to Workflow C's MakeFamilyInputs validation.
#
task ValidateFamilySet {
    input {
        String family_id
        Array[String] sample_ids
        File ped
        String docker_image
    }

    command <<<
        set -euxo pipefail

        cat > set_members.raw <<'EOF_SAMPLE_IDS'
~{sep="\n" sample_ids}
EOF_SAMPLE_IDS
        grep -v '^[[:space:]]*$' set_members.raw | sort -u > set_members.txt
        if [ ! -s set_members.txt ]; then
            echo "ERROR: sample set ~{family_id} has no members."
            exit 1
        fi

        awk -v fam="~{family_id}" 'BEGIN { FS="[ \t]+" } $1==fam && $2!="0" && $2!="." { print $2 }' ~{ped} | sort -u > ped_members.txt
        if [ ! -s ped_members.txt ]; then
            echo "ERROR: family_id '~{family_id}' not found in PED column 1 (or it lists no members). The Terra sample_set id must match the PED family id."
            exit 1
        fi

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
