#!/usr/bin/env python3
import sys
import random
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
           q, score, taxon, size, best_taxa_all):

    count[q] += 1
    taxon_set[q].add(taxon)

    if q not in best:
        best[q] = score
        second[q] = float("-inf")
        taxon_best[q] = taxon
        size_dict[q] = size
        best_taxa_all[q] = {taxon}
        return

    if score > best[q]:
        second[q] = best[q]
        best[q] = score
        taxon_best[q] = taxon
        size_dict[q] = size
        best_taxa_all[q] = {taxon}

    elif score == best[q]:
        best_taxa_all[q].add(taxon)

    elif score > second[q]:
        second[q] = score


def resolve_taxon(q, taxon_best, taxon_set, best, second, count, delta, best_taxa_all):
    """
    Final taxonomic resolution:

    Special mode:
    - delta < 0 → take best hit only
      if multiple best hits (equal score), pick randomly

    Standard mode:
    1. Species conflict = Legionella spp.
    2. Low score separation = strain=multi
    3. Clean assignment
    """

    # -------------------------
    # SPECIAL MODE
    # -------------------------
    if delta < 0:
        # taxons with same bitscore
        candidates = best_taxa_all[q]

        if len(candidates) > 1:
            return random.choice(list(candidates))
        else:
            return next(iter(candidates))

    best_taxon = taxon_best[q]

    # Case 1: different species detected and weak score separation
    if len(taxon_set[q]) > 1 and (best[q] - second[q]) <= delta:
        return "Legionella spp."

    # Case 2: same species but weak score separation
    if count[q] > 1 and (best[q] - second[q]) <= delta:
        return best_taxon + "|strain=multi"

    # Case 3: confident assignment
    return best_taxon


def write_output(outfile, best, second, taxon_best, taxon_set,
                 size, count, delta, best_taxa_all, all_queries):

    written = set()

    with open(outfile, "w") as out:

        # -------------------------
        # With hits
        # -------------------------
        for q in best:
            final_taxon = resolve_taxon(
                q, taxon_best, taxon_set, best, second,
                count, delta, best_taxa_all
            )
            out.write(f"{clean_qseqid(q)}\t{final_taxon}\t{size[q]}\n")
            written.add(q)

        # -------------------------
        # Not_assigned
        # -------------------------
        for q in all_queries:
            if q not in written:
                out.write(f"{clean_qseqid(q)}\tNot_assigned|Not_assigned\t{extract_size(q)}\n")


# -------------------------
# Main
# -------------------------

def main():
    random.seed(42) #Reproducitbility set

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
    best_taxa_all_s = defaultdict(set)
    all_q_s = set()

    # -------------------------
    # Containers (loose)
    # -------------------------
    best_l = {}
    second_l = {}
    tax_l = {}
    size_l = {}
    count_l = defaultdict(int)
    taxa_l = defaultdict(set)
    best_taxa_all_l = defaultdict(set)
    all_q_l = set()

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
            all_q_s.add(qseqid)
            all_q_l.add(qseqid)

            # -------------------------
            # STRICT FILTER
            # -------------------------
            if qlen >= minlen and qcovhsp >= qcov and pident >= pid:
                update(best_s, second_s, tax_s, size_s, taxa_s,
                       count_s, qseqid, bitscore, taxon, size, best_taxa_all_s)

            # -------------------------
            # LOOSE FILTER
            # -------------------------
            if (qlen >= minlen_loose and
                qcovhsp >= qcov_loose and
                pident >= pid_loose):

                update(best_l, second_l, tax_l, size_l, taxa_l,
                       count_l, qseqid, bitscore, taxon, size, best_taxa_all_l)

    # -------------------------
    # Write outputs
    # -------------------------
    write_output(strict_out, best_s, second_s, tax_s, taxa_s,
            size_s, count_s, delta, best_taxa_all_s, all_q_s)

    write_output(loose_out, best_l, second_l, tax_l, taxa_l,
            size_l, count_l, delta, best_taxa_all_l, all_q_l)


if __name__ == "__main__":
    main()