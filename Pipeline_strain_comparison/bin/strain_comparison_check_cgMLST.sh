#!/usr/bin/env bash

set -euo pipefail

################################################################################
# Extract selected cgMLST loci from a chewBBACA summary table
#
# Description:
#   Keep all samples from a chewBBACA summary table and retain only the
#   requested cgMLST loci. Report metadata samples from a given comparison
#   missing from the summary. Write an error flag if no sample of the
#   comparison is found in the summary.
#
# Arguments:
#   1 - chewBBACA summary table
#   2 - Metadata TSV
#   3 - Gene list (one locus per line)
#   4 - Comparison value
#   5 - Output TSV
#   6 - Log file
#   7 - Error file
#
# Requirements:
#   - bash
#   - awk (mawk or gawk)
################################################################################

SUMMARY="$1"
METADATA="$2"
GENES="$3"
COMPARISON="$4"
OUTPUT="$5"
LOG="$6"
ERROR_FILE="$7"

################################################################################
# Check input files
################################################################################

for file in "$SUMMARY" "$METADATA" "$GENES"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: File not found: $file" >&2
        exit 1
    fi
done

: > "$LOG"
: > "$ERROR_FILE"

################################################################################
# Extract loci and report missing samples
################################################################################

awk \
-v comparison="$COMPARISON" \
-v genes_file="$GENES" \
-v metadata_file="$METADATA" \
-v summary_file="$SUMMARY" \
-v log_file="$LOG" \
-v error_file="$ERROR_FILE" \
'

BEGIN {
    FS = OFS = "\t"
    printf "" > log_file
    close(log_file)
}

###############################################################################
# Read gene list
###############################################################################

FILENAME == genes_file {
    genes[$1] = 1
    next
}

###############################################################################
# Read metadata
###############################################################################

FILENAME == metadata_file {

    if (FNR == 1) {

        for (i = 1; i <= NF; i++) {
            if ($i == "Comparison")
                comparison_col = i
        }

        next
    }

    if ($comparison_col == comparison)
        metadata_samples[$1] = 1

    next
}

###############################################################################
# Read chewBBACA summary
###############################################################################

FILENAME == summary_file {

    if (FNR == 1) {

        keep[1] = 1

        for (i = 2; i <= NF; i++) {
            if ($i in genes)
                keep[i] = 1
        }

        first = 1

        for (i = 1; i <= NF; i++) {
            if (keep[i]) {
                if (!first)
                    printf OFS
                printf "%s", $i
                first = 0
            }
        }

        printf "\n"
        next
    }

    if ($1 in metadata_samples)
        found_comparison++

    found[$1] = 1

    first = 1

    for (i = 1; i <= NF; i++) {
        if (keep[i]) {
            if (!first)
                printf OFS
            printf "%s", $i
            first = 0
        }
    }

    printf "\n"
}

###############################################################################
# Report missing metadata samples and error flag
###############################################################################

END {

    for (id in metadata_samples) {
        if (!(id in found))
            print comparison, "SAMPLE MISSING FROM PREVIOUS", id >> log_file
    }

    if (found_comparison == 0) {
        print 1 > error_file
        print comparison, "NO ALLELIC PROFILE", "NA" >> log_file
    }
    else {
        print 0 > error_file
    }

}

' "$GENES" "$METADATA" "$SUMMARY" > "$OUTPUT"
