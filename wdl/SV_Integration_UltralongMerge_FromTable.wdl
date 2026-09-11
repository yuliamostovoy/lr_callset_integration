version 1.0

import "SV_Integration_UltralongMerge.wdl" as base


# DATA-TABLE FRONT-END for `SV_Integration_UltralongMerge`.
#
# Merge is an AGGREGATION: it pulls one per-sample VCF for every sample and
# concatenates them into a single cohort-wide file, so there is no batching /
# scatter — one submission consumes the whole sample list at once.
#
# The upstream workflow ingests the list through a `samples_csv` (one sample ID
# per line; extra CSV fields ignored). This wrapper takes `sample_ids` as a Terra
# data-table column and writes that CSV directly, then makes a single call to the
# unmodified `SV_Integration_UltralongMerge` workflow.
#
# Per sample it pulls `remote_indir/<sample_id>_<svtype><suffix>.vcf.gz` and
# writes `<svtype><suffix>_merged.vcf.gz` to `remote_outdir`. Run it once per
# (svtype, suffix) combination you need:
#   suffix=""          over the Annotate output dir       -> Score `input_vcf_gz`
#   suffix="_training" over the GetTrainingIntervals dir  -> Score `resource_vcf_gz`
#
workflow SV_Integration_UltralongMerge_FromTable {
    input {
        # --- Per-sample column from the Terra data table ---
        Array[String] sample_ids

        # --- Underlying-workflow pass-through ---
        String remote_indir
        String remote_outdir
        String svtype
        String suffix

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        sample_ids: "From the data table, e.g. `this.sample_id`. The whole set is merged in one run."
        remote_indir: "Folder holding the per-sample `<sample_id>_<svtype><suffix>.vcf.gz` to merge. Trailing slashes stripped."
        remote_outdir: "Where `<svtype><suffix>_merged.vcf.gz` is written. Trailing slashes stripped."
        svtype: "One of del, dup, ins, insdup, inv, bnd."
        suffix: "'' for the all-calls merge (Score input); '_training' for the labeled-true merge (Score resource)."
    }

    File samples_csv = write_lines(sample_ids)

    call base.SV_Integration_UltralongMerge as Merge {
        input:
            samples_csv = samples_csv,
            remote_indir = sub(remote_indir, "/+$", ""),
            remote_outdir = sub(remote_outdir, "/+$", ""),
            svtype = svtype,
            suffix = suffix,
            docker_image = docker_image
    }

    output {
        String merged_outdir = sub(remote_outdir, "/+$", "")
    }
}
