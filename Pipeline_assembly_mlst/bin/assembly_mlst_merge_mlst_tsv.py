#!/usr/bin/env python3

"""
Merge two TSV files with full column union.

Rules implemented:
- First column of each input is ALWAYS treated as ID (header ignored)
- Output first column is always named "ID"
- All columns from both files are kept (union of columns)
- Missing values are filled with empty strings
- Duplicate rows are removed
- Warning file is always created:
    - If no shared non-ID column exists: only message is written
    - Otherwise: list of columns missing per file
"""

import argparse
import pandas as pd


def main():
    """Main entry point: load TSVs, harmonize columns, merge and export."""

    # -------------------------
    # Argument parsing
    # -------------------------
    parser = argparse.ArgumentParser(
        description="Merge two TSV files keeping all columns."
    )

    parser.add_argument("--tsv-user", required=True)
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--output-name", required=True)
    parser.add_argument("--warning", required=True)

    args = parser.parse_args()

    # -------------------------
    # Load TSV files
    # -------------------------
    df_user = pd.read_csv(
        args.tsv_user,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    df_new = pd.read_csv(
        args.tsv,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    # -------------------------
    # Normalize ID column
    # First column is always ID regardless of header name
    # -------------------------
    df_user = df_user.rename(columns={df_user.columns[0]: "ID"})
    df_new = df_new.rename(columns={df_new.columns[0]: "ID"})

    # -------------------------
    # Column sets
    # -------------------------
    cols_user = list(df_user.columns)
    cols_new = list(df_new.columns)

    set_user = set(cols_user)
    set_new = set(cols_new)

    # -------------------------
    # Identify shared non-ID columns
    # (used only for warning logic, NOT filtering anymore)
    # -------------------------
    shared_cols = [c for c in cols_user if c in set_new and c != "ID"]

    # -------------------------
    # Warning collection
    # -------------------------
    warnings = []

    # Missing columns in new file
    for col in cols_user:
        if col not in set_new:
            warnings.append({"Header": col, "Absent_in": args.tsv})

    # Missing columns in user file
    for col in cols_new:
        if col not in set_user:
            warnings.append({"Header": col, "Absent_in": args.tsv_user})

    # Special case: no shared non-ID columns
    if len(shared_cols) == 0:
        with open(args.warning, "w") as f:
            f.write("No shared column was found between the input files.\n")

    else:
        pd.DataFrame(warnings, columns=["Header", "Absent_in"]).to_csv(
            args.warning,
            sep="\t",
            index=False,
        )

    # -------------------------
    # COLUMN UNION (NEW BEHAVIOR)
    # Keep all columns from both files
    # -------------------------
    all_cols = list(dict.fromkeys(cols_user + cols_new))

    # Ensure ID is first column
    all_cols = ["ID"] + [c for c in all_cols if c != "ID"]

    # -------------------------
    # Reindex both datasets to full schema
    # Missing columns -> empty string
    # -------------------------
    df_user = df_user.reindex(columns=all_cols, fill_value="")
    df_new = df_new.reindex(columns=all_cols, fill_value="")

    # -------------------------
    # Merge datasets
    # -------------------------
    merged = pd.concat(
        [df_user, df_new],
        ignore_index=True,
    )

    # Remove duplicate rows
    merged = merged.drop_duplicates()

    # -------------------------
    # Write output TSV
    # -------------------------
    merged.to_csv(
        args.output_name,
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
