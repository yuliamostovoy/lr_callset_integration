#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  make_workpackage7_chunk_id_files.sh <wp6_output_dir> [output_dir] [chunk_ids_per_file]

Arguments:
  wp6_output_dir       WP6 output directory. Supports gs://... or a local path.
                       Expected structure: <wp6_output_dir>/<chrom>/regions.txt
  output_dir           Local directory for WP7 chunk-id files. Default: .
  chunk_ids_per_file   Max IDs per file. Default: 100

Output files are named:
  workpackage7_<chrom>_0000
  workpackage7_<chrom>_0001
  ...

Each output file contains one WP6 chunk ID per line.
EOF
}

if [[ $# -lt 1 || $# -gt 3 ]]; then
    usage
    exit 1
fi

WP6_OUTPUT_DIR="${1%/}"
OUTPUT_DIR="${2:-.}"
CHUNK_IDS_PER_FILE="${3:-100}"

if ! [[ "${CHUNK_IDS_PER_FILE}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: chunk_ids_per_file must be a positive integer: ${CHUNK_IDS_PER_FILE}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

is_gcs_path() {
    [[ "$1" == gs://* ]]
}

list_chromosomes_gcs() {
    gcloud storage ls "${WP6_OUTPUT_DIR}/" |
        awk -F/ 'NF >= 2 { print $(NF-1) }' |
        grep -v '^$' |
        sort -V
}

list_chromosomes_local() {
    find "${WP6_OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -type d -print |
        xargs -n 1 basename |
        sort -V
}

copy_regions_txt() {
    local chrom="$1"
    local local_regions="$2"

    if is_gcs_path "${WP6_OUTPUT_DIR}"; then
        gcloud storage cp "${WP6_OUTPUT_DIR}/${chrom}/regions.txt" "${local_regions}"
    else
        cp "${WP6_OUTPUT_DIR}/${chrom}/regions.txt" "${local_regions}"
    fi
}

CHROMS=()
if is_gcs_path "${WP6_OUTPUT_DIR}"; then
    while IFS= read -r CHROM; do
        CHROMS+=("${CHROM}")
    done < <(list_chromosomes_gcs)
else
    while IFS= read -r CHROM; do
        CHROMS+=("${CHROM}")
    done < <(list_chromosomes_local)
fi

if [[ ${#CHROMS[@]} -eq 0 ]]; then
    echo "ERROR: No chromosome directories found under ${WP6_OUTPUT_DIR}" >&2
    exit 1
fi

for CHROM in "${CHROMS[@]}"; do
    REGIONS_TXT="${TMP_DIR}/${CHROM}.regions.txt"
    IDS_FILE="${TMP_DIR}/${CHROM}.chunk_ids.txt"

    copy_regions_txt "${CHROM}" "${REGIONS_TXT}"

    awk 'NF >= 2 { print $2 }' "${REGIONS_TXT}" > "${IDS_FILE}"
    N_IDS="$(wc -l < "${IDS_FILE}" | tr -d ' ')"
    if [[ "${N_IDS}" -eq 0 ]]; then
        echo "ERROR: No chunk IDs found in ${WP6_OUTPUT_DIR}/${CHROM}/regions.txt" >&2
        exit 1
    fi

    split -l "${CHUNK_IDS_PER_FILE}" -d -a 4 \
        "${IDS_FILE}" "${OUTPUT_DIR}/workpackage7_${CHROM}_"

    N_FILES="$(find "${OUTPUT_DIR}" -maxdepth 1 -type f -name "workpackage7_${CHROM}_[0-9][0-9][0-9][0-9]" | wc -l | tr -d ' ')"
    echo "${CHROM}: wrote ${N_FILES} file(s) for ${N_IDS} chunk ID(s)" >&2
done
