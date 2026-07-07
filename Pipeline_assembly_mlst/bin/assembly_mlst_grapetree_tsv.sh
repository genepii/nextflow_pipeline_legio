#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Convert coma to tabulation (remove ST column for GrapeTree input)
# AND generate metadata (Sample_ID + ST)
# ------------------------------------------------------------

INPUT_CSV="$1"
OUTPUT_TSV="$2"
OUTPUT_META="$3"

awk -F ',' -v OFS='\t' '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "ST") st=i
    }
    print
    next
}
{
    for (i=1; i<=NF; i++) {
        if (i != st) {
            printf "%s%s", $i, (i==NF ? ORS : OFS)
        }
    }
}
' "$INPUT_CSV" > "$OUTPUT_TSV"

awk -F '\t' '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "ST") st=i
    }
    print "Sample_ID\tST"
    next
}
{
    print $1 "\t" $st
}
' "$INPUT_CSV" > "$OUTPUT_META"