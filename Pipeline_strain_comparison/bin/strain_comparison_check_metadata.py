#!/usr/bin/env python3

"""
Validate metadata completeness and split samples according to
available FASTA and FASTQ files.

Outputs:
- Samples with missing files and missing file description
- Samples with an available FASTA file
- Samples with available paired FASTQ files
"""


import argparse
import pandas as pd


def parse_arguments():

    """
    Parse command line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Check metadata completeness"
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Input metadata TSV file"
    )

    parser.add_argument(
        "--missing",
        required=True,
        help="Output TSV containing samples with missing files"
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Output TSV containing samples with FASTA files"
    )

    parser.add_argument(
        "--fastq",
        required=True,
        help="Output TSV containing samples with paired FASTQ files"
    )

    return parser.parse_args()


def get_missing(row):

    """
    Identify missing sequencing files for one sample.

    Returns:
        Semicolon-separated list of missing files.
    """

    missing = []

    if row["FASTA"] == "NA":
        missing.append("FASTA")

    if row["READS_1"] == "NA":
        missing.append("READS_1")

    if row["READS_2"] == "NA":
        missing.append("READS_2")

    return ";".join(missing)


def main():

    args = parse_arguments()


    # Load metadata containing FASTA and FASTQ paths
    metadata = pd.read_csv(
        args.metadata,
        sep="\t",
        dtype=str
    ).fillna("NA")


    #
    # Identify samples with at least one missing file
    #
    missing = metadata[
        (metadata["FASTA"] == "NA") |
        (metadata["READS_1"] == "NA") |
        (metadata["READS_2"] == "NA")
    ].copy()


    # Add a description of missing files
    missing["MISSING"] = missing.apply(
        get_missing,
        axis=1
    )


    # Keep only metadata fields useful for reporting
    missing = missing[
        [
            "ID",
            "ST",
            "Year",
            "Origin",
            "MISSING"
        ]
    ]


    #
    # Select samples with an available FASTA assembly
    #
    fasta = metadata[
        metadata["FASTA"] != "NA"
    ]


    #
    # Select samples with both paired FASTQ files available
    #
    fastq = metadata[
        (metadata["READS_1"] != "NA") &
        (metadata["READS_2"] != "NA")
    ]


    # Write samples with missing files
    missing.to_csv(
        args.missing,
        sep="\t",
        index=False
    )


    # Write samples with FASTA files
    fasta.to_csv(
        args.fasta,
        sep="\t",
        index=False
    )


    # Write samples with paired FASTQ files
    fastq.to_csv(
        args.fastq,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
