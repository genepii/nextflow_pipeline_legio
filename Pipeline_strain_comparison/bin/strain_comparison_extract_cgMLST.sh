#!/usr/bin/env bash

set -euo pipefail

################################################################################
# Extract selected cgMLST loci from a chewBBACA summary table
#
# Description:
#   Extract samples belonging to a given comparison from a metadata file and
#   keep only the requested cgMLST loci from a chewBBACA summary table.
#
# Arguments:
#   1 - chewBBACA summary table
#   2 - Metadata TSV
#   3 - Gene list (one locus per line)
#   4 - Comparison value
#   5 - Output TSV
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

################################################################################
# Check input files
################################################################################

for file in "$SUMMARY" "$METADATA" "$GENES"; do

    if [[ ! -f "$file" ]]; then
        echo "ERROR: File not found: $file" >&2
        exit 1
    fi

done

################################################################################
# Extract samples
################################################################################

awk \
-v comparison="$COMPARISON" \
-v genes_file="$GENES" \
-v metadata_file="$METADATA" \
-v summary_file="$SUMMARY" \
'

BEGIN {
    FS = OFS = "\t"
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

        comparison_col = 0

        for (i = 1; i <= NF; i++) {

            if ($i == "Comparison") {
                comparison_col = i
                break
            }

        }

        if (comparison_col == 0) {

            print "ERROR: Column \"Comparison\" not found." > "/dev/stderr"
            exit 1

        }

        next
    }

    if ($comparison_col == comparison)
        samples[$1] = 1

    next
}

###############################################################################
# Read chewBBACA summary
###############################################################################

FILENAME == summary_file {

    ###########################################################################
    # Header
    ###########################################################################

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

    ###########################################################################
    # Selected samples
    ###########################################################################

    if (!($1 in samples))
        next

    if ($1 in written)
        next

    written[$1] = 1

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

' "$GENES" "$METADATA" "$SUMMARY" > "$OUTPUT"
