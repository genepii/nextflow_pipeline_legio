#!/usr/bin/env python3

"""
Merge FASTA and FASTQ metadata tables.
If no Comparison value => default.
"""


import argparse
import pandas as pd


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Merge metadata tables"
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Metadata table containing FASTA paths"
    )

    parser.add_argument(
        "--fastq",
        required=True,
        help="Metadata table containing FASTQ paths"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output metadata table"
    )

    return parser.parse_args()


def main():

    args = parse_arguments()

    # Load metadata tables
    fasta_df = pd.read_csv(
        args.fasta,
        sep="\t",
        dtype=str
    )

    fastq_df = pd.read_csv(
        args.fastq,
        sep="\t",
        dtype=str
    )

    # Fill missing values in Comparison column
    if "Comparison" in fasta_df.columns:
        fasta_df["Comparison"] = (
            fasta_df["Comparison"]
            .replace(r"^\s*$", "default", regex=True)
            .fillna("default")
        )

    # Fill remaining missing values
    fasta_df = fasta_df.fillna("NA")
    fastq_df = fastq_df.fillna("NA")

    # Keep only FASTQ columns to avoid duplicating metadata
    fastq_df = fastq_df[
        [
            "ID",
            "READS_1",
            "READS_2"
        ]
    ]

    # Merge FASTQ paths into FASTA metadata
    metadata = fasta_df.merge(
        fastq_df,
        on="ID",
        how="left"
    )

    # Write complete metadata table
    metadata.to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
