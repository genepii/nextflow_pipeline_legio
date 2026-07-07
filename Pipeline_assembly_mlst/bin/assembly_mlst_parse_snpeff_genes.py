#!/usr/bin/env python3

"""
assembly_mlst_parse_snpeff_genes.py

Parse SnpEff-annotated VCF files and extract gene-level annotations
for downstream AMR analysis.

Adds AF (allele frequency) extracted from VCF INFO field.
"""

import argparse
import gzip
import re
import pandas as pd


# Output schema (final TSV columns)
OUTPUT_COLUMNS = [
    "sample_id",
    "gene_id",
    "transcript_id",
    "bio_type",
    "variants_impact_MODIFIER",
    "variants_effect_non_coding_transcript_variant",
    "AF"
]


def open_vcf(vcf_path):
    """
    Open VCF file, supporting both plain text and gzip compressed formats.

    Parameters
    ----------
    vcf_path : str
        Path to VCF file (.vcf or .vcf.gz)

    Returns
    -------
    file handle
        Iterable file object in text mode
    """
    if vcf_path.endswith(".gz"):
        return gzip.open(vcf_path, "rt")

    return open(vcf_path, "r")


def extract_af(info_field):
    """
    Extract allele frequency (AF) from VCF INFO field.

    Parameters
    ----------
    info_field : str
        INFO column of a VCF line

    Returns
    -------
    str
        AF value if present, otherwise "NA"
    """
    # Match AF= value anywhere in INFO, including after other tags
    match = re.search(r"(?:^|;)AF=([^;]+)", info_field)
    return match.group(1) if match else "NA"


def parse_ann_field(ann_field, sample_id, af_value, rows):
    """
    Parse SnpEff ANN field and append structured annotations.

    Each ANN entry corresponds to a transcript-level annotation.
    Multiple ANN entries may exist per variant.
    """

    ann_entries = ann_field.split(",")

    for ann in ann_entries:

        fields = ann.split("|")

        # SnpEff ANN has multiple subfields; skip malformed entries
        if len(fields) < 8:
            continue

        # Extract key biological annotation fields
        effect = fields[1]         # Variant effect (e.g. missense_variant)
        impact = fields[2]         # Impact level (LOW, MODERATE, MODIFIER, HIGH)
        gene_id = fields[4]        # Gene identifier
        transcript_id = fields[6]  # Transcript identifier
        bio_type = fields[7]       # Gene/transcript biotype

        # Build output row (one per ANN entry)
        rows.append({
            "sample_id": sample_id,
            "gene_id": gene_id if gene_id else "NA",
            "transcript_id": transcript_id if transcript_id else "NA",
            "bio_type": bio_type if bio_type else "NA",

            # Binary flags derived from annotation logic
            "variants_impact_MODIFIER": int(impact == "MODIFIER"),
            "variants_effect_non_coding_transcript_variant": int(
                effect == "non_coding_transcript_variant"
            ),

            # Allele frequency propagated from VCF INFO
            "AF": af_value
        })


def parse_args():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        prog="assembly_mlst_parse_snpeff_genes.py",
        description="Extract gene-level annotations from SnpEff VCF"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input SnpEff annotated VCF (.vcf or .vcf.gz)"
    )

    parser.add_argument(
        "-s", "--sample",
        required=True,
        help="Sample identifier"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )

    return parser.parse_args()


def main():
    """
    Main pipeline:
    - Parse VCF line by line
    - Extract AF from INFO field
    - Extract ANN annotations
    - Flatten into tabular gene-level structure
    """

    args = parse_args()

    rows = []

    # Iterate over VCF records
    with open_vcf(args.input) as handle:

        for line in handle:

            # Skip metadata/header lines
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            # Ensure valid VCF structure (at least 8 columns)
            if len(cols) < 8:
                continue

            info = cols[7]

            # Extract allele frequency from INFO field
            af_value = extract_af(info)

            # Extract SnpEff annotation block (ANN field)
            ann_match = re.search(r"ANN=([^;]+)", info)

            if not ann_match:
                continue

            # Parse ANN entries and append results
            parse_ann_field(
                ann_match.group(1),
                args.sample,
                af_value,
                rows
            )

    # Ensure downstream compatibility even if no variants found
    if not rows:
        rows.append({
            "sample_id": args.sample,
            "gene_id": "NA",
            "transcript_id": "NA",
            "bio_type": "NA",
            "variants_impact_MODIFIER": "NA",
            "variants_effect_non_coding_transcript_variant": "NA",
            "AF": "NA"
        })

    # Convert to DataFrame and remove duplicates
    df = pd.DataFrame(rows, columns=OUTPUT_COLUMNS).drop_duplicates()

    # Write final TSV
    df.to_csv(args.output, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()