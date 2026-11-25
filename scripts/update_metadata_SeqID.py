#!/usr/bin/env python3

"""
Update Seq_ID and Analyse_ID using the most recent matching
Summary_table.tsv files.

Directory structure:
    /genomes/XXXX/YYYY_Assembly-MLST/Summary_table.tsv

Selection priority:
    1. Most recent YYYY (analysis date).
    2. Most recent XXXX (sequencing date).
"""

import argparse
import re
from pathlib import Path

import pandas as pd


DATE_PATTERN = re.compile(r"^\d{8}$")


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Update Seq_ID and Analyse_ID from Summary_table.tsv files."
    )
    parser.add_argument("--input", required=True, help="Input TSV file.")
    parser.add_argument("--genomes-dir", required=True, help="Genomes root directory.")
    parser.add_argument("--output", required=True, help="Output TSV file.")
    return parser.parse_args()


def find_summary_files(genomes_dir):
    """Find valid Summary_table.tsv files and extract their dates."""
    results = []

    for summary_file in Path(genomes_dir).glob(
        "????????/????????_Assembly-MLST/Summary_table.tsv"
    ):
        seq_id = summary_file.parent.parent.name
        analyse_id = summary_file.parent.name.replace(
            "_Assembly-MLST", ""
        )

        if DATE_PATTERN.match(seq_id) and DATE_PATTERN.match(analyse_id):
            results.append({
                "seq_id": seq_id,
                "analyse_id": analyse_id,
                "path": summary_file,
            })

    return results


def build_id_index(summary_files):
    """Build an index of sample IDs and their associated dates."""
    id_index = {}

    for summary in summary_files:
        try:
            df = pd.read_csv(summary["path"], sep="\t", dtype=str)
        except (OSError, pd.errors.ParserError) as error:
            print(f"WARNING: Cannot read {summary['path']}: {error}")
            continue

        if "Sample_ID" not in df.columns:
            print(f"WARNING: Missing ID column in {summary['path']}")
            continue

        for sample_id in df["Sample_ID"].dropna().unique():
            sample_id = sample_id.strip()

            if not sample_id:
                continue

            candidate = (
                summary["analyse_id"],
                summary["seq_id"],
            )

            if sample_id not in id_index or candidate > id_index[sample_id]["dates"]:
                id_index[sample_id] = {
                    "dates": candidate,
                    "seq_id": summary["seq_id"],
                    "analyse_id": summary["analyse_id"],
                }

    return id_index


def main():
    """Run the metadata update."""
    args = parse_arguments()

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    required_columns = {"ID", "Seq_ID", "Analyse_ID"}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {', '.join(sorted(missing_columns))}"
        )

    summary_files = find_summary_files(args.genomes_dir)
    id_index = build_id_index(summary_files)

    updated_count = 0

    for index, row in df.iterrows():
        sample_id = str(row["ID"]).strip()

        if sample_id in id_index:
            df.at[index, "Seq_ID"] = id_index[sample_id]["seq_id"]
            df.at[index, "Analyse_ID"] = id_index[sample_id]["analyse_id"]
            updated_count += 1

    df.to_csv(args.output, sep="\t", index=False, na_rep="NA")

    print(f"Summary files found: {len(summary_files)}")
    print(f"IDs indexed: {len(id_index)}")
    print(f"Rows updated: {updated_count}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
