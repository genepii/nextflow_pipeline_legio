#!/usr/bin/env python3

import sys


def compute_ratio(count, total_reads):
    return count / total_reads if total_reads > 0 else 0.0


def process_kraken(sample_id, total_reads, kraken_file, output_tsv,
                   target_species, target_threshold, legio_threshold, legia_threshold):

    total_reads = float(total_reads)
    target_threshold = float(target_threshold)
    legio_threshold = float(legio_threshold)
    legia_threshold = float(legia_threshold)

    target_ratio = 0.0
    legio_ratio = 0.0
    legia_ratio = 0.0

    found_target = False
    found_legio = False
    found_legia = False

    with open(kraken_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 2:
                continue

            taxonomy = cols[0]

            try:
                count = float(cols[1])
            except ValueError:
                continue

            ratio = compute_ratio(count, total_reads)

            # Target species match
            if target_species in taxonomy and not found_target:
                target_ratio = ratio
                found_target = True

            # Legionella genus match
            if "g__Legionella" in taxonomy and not found_legio:
                legio_ratio = ratio
                found_legio = True

            # Legionellaceae family match
            if "f__Legionellaceae" in taxonomy and not found_legia:
                legia_ratio = ratio
                found_legia = True

            if found_target and found_legio and found_legia:
                break

    # Classification rules
    if target_ratio >= target_threshold:
        classification = target_species
    elif legio_ratio >= legio_threshold:
        classification = "Legionella spp."
    elif legia_ratio >= legia_threshold:
        classification = "Legionella spp. (extended)"
    elif target_ratio < target_threshold and legio_ratio < legio_threshold:
        classification = "Contamination"
    else:
        classification = "Contamination"

    with open(output_tsv, "w") as out:
        out.write(
            f"{sample_id}\t{target_ratio}\t{legio_ratio}\t{classification}\n"
        )

    return True


if __name__ == "__main__":
    print("Arguments received:", sys.argv)

    if len(sys.argv) != 9:
        print(
            "Usage: assembly_mlst_kraken2_ratio.py <sample_id> <total_reads> <kraken_file> <output_tsv> "
            "<target_species> <target_threshold> <legio_threshold> <legia_threshold>",
            file=sys.stderr
        )
        sys.exit(1)

    sample_id = sys.argv[1]
    total_reads = sys.argv[2]
    kraken_file = sys.argv[3]
    output_tsv = sys.argv[4]
    target_species = sys.argv[5]
    target_threshold = sys.argv[6]
    legio_threshold = sys.argv[7]
    legia_threshold = sys.argv[8]

    process_kraken(
        sample_id,
        total_reads,
        kraken_file,
        output_tsv,
        target_species,
        target_threshold,
        legio_threshold,
        legia_threshold
    )
