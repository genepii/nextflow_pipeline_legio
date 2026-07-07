#!/usr/bin/env python3
# ------------------------------------------------------------
# merge_alleles.py
#
# Merge chewBBACA allele TSV files into a single matrix.
# Missing values are filled with NA.
#
# Usage:
#   python merge_alleles.py -o output.tsv file1.tsv file2.tsv ...
# ------------------------------------------------------------

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge chewBBACA allele TSV files into a matrix"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )
    parser.add_argument(
        "files",
        nargs="+",
        help="Input TSV files"
    )
    return parser.parse_args()


def load_file(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)

    # Ensure Sample_ID exists
    if "Sample_ID" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Sample_ID"})

    return df


def main():
    args = parse_args()

    dfs = []
    for f in args.files:
        df = load_file(f)
        dfs.append(df)

    # Full outer join on Sample_ID (union of all loci)
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on="Sample_ID", how="outer")

    # Fill missing values
    merged = merged.fillna("NA")

    # Optional: stable sorting
    merged = merged.sort_values("Sample_ID")

    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()