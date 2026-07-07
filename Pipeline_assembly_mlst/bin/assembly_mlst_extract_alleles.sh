#!/usr/bin/env bash

set -euo pipefail

mapping="$1"
input_tsv="$2"
output_tsv="$3"

awk -F '\t' -v MAP="${mapping}" '
BEGIN {

    OFS="\t"

    # Parse mapping string:
    # lpLIT0656:lag,lpLIT2833:lpeA,...

    n = split(MAP, arr, ",")

    for(i=1; i<=n; i++) {

        split(arr[i], kv, ":")

        wanted[kv[1]] = kv[2]
    }
}

NR==1 {

    out = $1

    for(i=2; i<=NF; i++) {

        if($i in wanted) {

            keep[++k] = i
            header[i] = wanted[$i]

            out = out OFS wanted[$i]
        }
    }

    print out
    next
}

{

    out = $1

    for(i=1; i<=k; i++) {

        col = keep[i]

        out = out OFS $(col)
    }

    print out
}
' "${input_tsv}" > "${output_tsv}"