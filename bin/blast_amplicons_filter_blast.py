#!/usr/bin/env python3
import sys
from collections import defaultdict

# -------------------------
# Helpers
# -------------------------

def extract_size(qseqid: str) -> int:
    """
    Extract read size from FASTA header (format: ...;size=XXX).
    """
    for part in qseqid.split(";"):
        if part.startswith("size="):
            try:
                return int(part.split("=")[1])
            except ValueError:
                return 1
    return 1


def clean_qseqid(qseqid: str) -> str:
    """
    Remove metadata from query ID (keep only primary identifier).
    """
    return qseqid.split(";")[0]


def update(best, second, taxon_best, size_dict, taxon_set, count,
           q, score, taxon, size):
    """
    Track best and second best bitscore per query,
    and store all observed taxa for ambiguity detection.
    """
    count[q] += 1
    taxon_set[q].add(taxon)

    if q not in best:
        best[q] = score
        second[q] = float("-inf")
        taxon_best[q] = taxon
        size_dict[q] = size
        return

    if score > best[q]:
        second[q] = best[q]
        best[q] = score
        taxon_best[q] = taxon
        size_dict[q] = size

    elif score > second[q]:
        second[q] = score


def resolve_taxon(q, taxon_best, taxon_set, best, second, count, delta):
    """
    Final taxonomic resolution:

    Priority:
    1. Species conflict = ambiguous
    2. Low score separation = strain=multi
    3. Clean assignment
    """

    best_taxon = taxon_best[q]

    # Case 1: different species detected and weak score separation
    if len(taxon_set[q]) > 1 and (best[q] - second[q]) <= delta:
        return "L-ambiguous"

    # Case 2: same species but weak score separation
    if count[q] > 1 and (best[q] - second[q]) <= delta:
        return best_taxon + "|strain=multi"

    # Case 3: confident assignment
    return best_taxon


def write_output(outfile, best, second, taxon_best, taxon_set,
                 size, count, delta):
    """
    Write final filtered BLAST table:
    qseqid, taxon, size
    """
    with open(outfile, "w") as out:
        for q in best:
            final_taxon = resolve_taxon(
                q, taxon_best, taxon_set, best, second, count, delta
            )
            out.write(f"{clean_qseqid(q)}\t{final_taxon}\t{size[q]}\n")


# -------------------------
# Main
# -------------------------

def main():
    if len(sys.argv) != 11:
        sys.exit(
            "Usage: blast_amplicons_filter_blast.py <blast> <pid> <qcov> <minlen> <delta> "
            "<pid_loose> <qcov_loose> <minlen_loose> <strict_out> <loose_out>"
        )

    blast_file = sys.argv[1]

    pid = float(sys.argv[2])
    qcov = float(sys.argv[3])
    minlen = int(sys.argv[4])
    delta = float(sys.argv[5])

    pid_loose = float(sys.argv[6])
    qcov_loose = float(sys.argv[7])
    minlen_loose = int(sys.argv[8])

    strict_out = sys.argv[9]
    loose_out = sys.argv[10]

    # -------------------------
    # Containers (strict)
    # -------------------------
    best_s = {}
    second_s = {}
    tax_s = {}
    size_s = {}
    count_s = defaultdict(int)
    taxa_s = defaultdict(set)

    # -------------------------
    # Containers (loose)
    # -------------------------
    best_l = {}
    second_l = {}
    tax_l = {}
    size_l = {}
    count_l = defaultdict(int)
    taxa_l = defaultdict(set)

    # -------------------------
    # Parse BLAST file
    # -------------------------
    with open(blast_file) as f:
        for line in f:
            cols = line.rstrip().split("\t")

            qseqid = cols[0]
            taxon = cols[1]
            pident = float(cols[2])
            qlen = int(cols[4])
            qcovhsp = float(cols[6])
            bitscore = float(cols[7])

            size = extract_size(qseqid)

            # -------------------------
            # STRICT FILTER
            # -------------------------
            if qlen >= minlen and qcovhsp >= qcov and pident >= pid:
                update(best_s, second_s, tax_s, size_s, taxa_s,
                       count_s, qseqid, bitscore, taxon, size)

            # -------------------------
            # LOOSE FILTER
            # -------------------------
            if (qlen >= minlen_loose and
                qcovhsp >= qcov_loose and
                pident >= pid_loose):

                update(best_l, second_l, tax_l, size_l, taxa_l,
                       count_l, qseqid, bitscore, taxon, size)

    # -------------------------
    # Write outputs
    # -------------------------
    write_output(strict_out, best_s, second_s, tax_s, taxa_s,
                 size_s, count_s, delta)

    write_output(loose_out, best_l, second_l, tax_l, taxa_l,
                 size_l, count_l, delta)


if __name__ == "__main__":
    main()