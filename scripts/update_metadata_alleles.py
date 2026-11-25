#!/usr/bin/env python3

"""
Update lag, lpeA and lpeB metadata using an allele matrix.

The script:
    - Reads a metadata TSV file.
    - Reads an allele matrix TSV file.
    - Matches samples using the ID column.
    - Updates lag, lpeA and lpeB from the corresponding lpLIT columns.
    - Keeps existing values when source values are missing.
    - Replaces allele value 0 with LNF.
    - Preserves all other metadata columns.
"""

import argparse
import pandas as pd


COLUMN_MAPPING = {
    "lag": "lpLIT0656",
    "lpeA": "lpLIT2833",
    "lpeB": "lpLIT2834",
}


MISSING_VALUES = {"", "NA", "NaN", "nan", "None", "NULL"}


def is_missing(value):
    """Check whether a value represents missing data."""
    if pd.isna(value):
        return True

    return str(value).strip() in MISSING_VALUES


def normalize_value(value):
    """Convert zero values to LNF and preserve other values."""
    value = str(value).strip()

    if is_missing(value):
        return None

    try:
        if float(value) == 0:
            return "LNF"
    except ValueError:
        pass

    return value


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Update metadata using an allele matrix."
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Input metadata TSV file.",
    )

    parser.add_argument(
        "--alleles",
        required=True,
        help="Input allele matrix TSV file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output metadata TSV file.",
    )

    return parser.parse_args()


def main():
    """Run the metadata update workflow."""
    args = parse_arguments()

    metadata = pd.read_csv(
        args.metadata,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    alleles = pd.read_csv(
        args.alleles,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    # Validate required columns.
    required_metadata = {"ID", *COLUMN_MAPPING.keys()}
    required_alleles = {"ID", *COLUMN_MAPPING.values()}

    missing_metadata = required_metadata - set(metadata.columns)
    missing_alleles = required_alleles - set(alleles.columns)

    if missing_metadata:
        raise ValueError(
            f"Missing columns in metadata file: {sorted(missing_metadata)}"
        )

    if missing_alleles:
        raise ValueError(
            f"Missing columns in allele file: {sorted(missing_alleles)}"
        )

    # Prevent ambiguous sample matching.
    if alleles["ID"].duplicated().any():
        duplicates = alleles.loc[
            alleles["ID"].duplicated(keep=False), "ID"
        ].unique()

        raise ValueError(
            f"Duplicate IDs in allele file: {list(duplicates)}"
        )

    # Index allele data by sample ID.
    alleles = alleles.set_index("ID")

    # Update metadata values.
    updated_count = 0

    for metadata_column, allele_column in COLUMN_MAPPING.items():

        source_values = metadata["ID"].map(alleles[allele_column])

        for index, source_value in source_values.items():

            if is_missing(source_value):
                continue

            new_value = normalize_value(source_value)

            if new_value is not None:
                metadata.at[index, metadata_column] = new_value
                updated_count += 1

    # Write output.
    metadata.to_csv(
        args.output,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(f"Metadata successfully updated: {args.output}")
    print(f"Number of updated cells: {updated_count}")


if __name__ == "__main__":
    main()
