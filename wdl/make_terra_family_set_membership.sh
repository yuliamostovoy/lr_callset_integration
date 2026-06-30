#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 INPUT_SAMPLE_TABLE OUTPUT_MEMBERSHIP_TSV [ENTITY_TYPE]" >&2
    exit 1
fi

INPUT_TABLE=$1
OUTPUT_TSV=$2
ENTITY_TYPE=${3:-}

if [ -z "${ENTITY_TYPE}" ]; then
    HEADER=$(head -n 1 "${INPUT_TABLE}" | cut -f 1)
    ENTITY_TYPE=${HEADER#entity:}
    ENTITY_TYPE=${ENTITY_TYPE%_id}
fi

awk -v entity_type="${ENTITY_TYPE}" '
    BEGIN {
        FS = OFS = "\t"
        print "membership:" entity_type "_set_id", entity_type
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "family_id") {
                family_col = i
            }
        }
        if (family_col == "") {
            print "ERROR: no family_id column found" > "/dev/stderr"
            exit 1
        }
        next
    }
    $family_col != "" && $family_col != "." {
        print $family_col, $1
    }
' "${INPUT_TABLE}" > "${OUTPUT_TSV}"
