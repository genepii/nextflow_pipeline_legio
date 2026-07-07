#!/usr/bin/env python3
# ------------------------------------------------------------
# merge_snpeff_genes.py
#
# Purpose:
# Aggregate SnpEff gene annotation files per sample.
#
# Features:
# - Ensures unique (GeneId, TranscriptId) per sample
# - Handles empty or NA/NA-only files
# - Optional column prefix via -t/--tag
# - Sample_ID is never prefixed (kept stable for joins)
# ------------------------------------------------------------

import argparse
import pandas as pd
import os


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )

    parser.add_argument(
        "-t", "--tag",
        type=str,
        default=None,
        help="Optional prefix added to output column names (except Sample_ID)"
    )

    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input SnpEff gene annotation files"
    )

    return parser.parse_args()


def extract_sample_id(path):
    """
    Extract sample ID from file name.
    Example:
        ABC.SnpEff.genes.txt -> ABC
    """
    base = os.path.basename(path)
    return base.replace(".SnpEff.genes.txt", "")


def parse_file(path):
    """
    Parse a single SnpEff genes file.

    Returns:
        list of rows [sample_id, gene_id, transcript_id, impact, noncoding]
    """
    sample_id = extract_sample_id(path)
    rows = []
    has_valid_row = False

    with open(path) as f:
        for line in f:

            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            # Ensure minimum column length
            parts = parts + ["NA"] * (6 - len(parts))

            gene_id = parts[1]
            transcript_id = parts[2]
            impact = parts[4]
            noncoding = parts[5]

            rows.append([
                sample_id,
                gene_id,
                transcript_id,
                impact,
                noncoding
            ])

            # Track if at least one real annotation exists
            if gene_id != "NA" or transcript_id != "NA":
                has_valid_row = True

    # If file is empty or only NA/NA pairs exist, return neutral row
    if not has_valid_row:
        return [[sample_id, None, None, 0, 0]]

    return rows


def main():
    args = parse_args()

    # Base column names (always used internally)
    base_cols = ["Sample_ID", "GeneId", "TranscriptId", "Impact", "NonCoding"]

    rows = []

    # Parse all input files
    for file_path in args.inputs:
        rows.extend(parse_file(file_path))

    # Build dataframe with stable schema
    df = pd.DataFrame(rows, columns=base_cols)

    # Convert numeric fields safely
    df["Impact"] = pd.to_numeric(df["Impact"], errors="coerce").fillna(0)
    df["NonCoding"] = pd.to_numeric(df["NonCoding"], errors="coerce").fillna(0)

    # Remove fully NA gene/transcript pairs
    df_valid = df[~((df["GeneId"].isna()) & (df["TranscriptId"].isna()))]

    # ------------------------------------------------------------
    # Aggregate at gene/transcript pair level
    # ------------------------------------------------------------
    df_pair = (
        df_valid
        .groupby(["Sample_ID", "GeneId", "TranscriptId"], as_index=False)
        .agg(
            Impact=("Impact", "sum"),
            NonCoding=("NonCoding", "sum")
        )
    )

    # ------------------------------------------------------------
    # Final aggregation per sample
    # ------------------------------------------------------------
    summary = (
        df_pair.groupby("Sample_ID")
        .agg(
            Nb_Mutated_Genes=("GeneId", "size"),
            Nb_with_Impact=("Impact", "sum"),
            Nb_Non_Coding=("NonCoding", "sum")
        )
        .reset_index()
    )

    # Ensure all samples are present even if empty
    all_samples = pd.DataFrame({
        "Sample_ID": [extract_sample_id(p) for p in args.inputs]
    })

    summary = all_samples.merge(summary, on="Sample_ID", how="left")
    summary = summary.fillna(0)

    # ------------------------------------------------------------
    # Apply optional prefix to ALL columns except Sample_ID
    # ------------------------------------------------------------
    if args.tag:
        summary.columns = [
            c if c == "Sample_ID" else f"{args.tag}_{c}"
            for c in summary.columns
        ]

    # Write final TSV
    summary.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()