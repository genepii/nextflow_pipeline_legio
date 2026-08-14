#!/usr/bin/env python3

"""
Merge two TSV files with full column union.

Rules implemented:
- First column of each input is ALWAYS treated as ID (header ignored).
- Output first column is always named "ID".
- All columns from both files are kept (union of columns).
- Missing values are filled with 0.
- If the same ID exists in both files:
    - If the profiles are identical, only one row is kept.
    - If the profiles are different, only the row from --tsv is kept.
- Rows with IDs present in only one file are kept.
- Duplicate rows are removed.
- Warning file is always created:
    - If no shared non-ID column exists: only a message is written.
    - Otherwise: list of columns missing per file.
"""

import argparse

import pandas as pd


def main():
    """Load TSV files, harmonize columns, resolve conflicting IDs, and export."""

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
    ).replace("", "0")

    df_new = pd.read_csv(
        args.tsv,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    ).replace("", "0")

    # -------------------------
    # Normalize ID column
    # First column is always treated as ID,
    # regardless of its original header.
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
    # -------------------------
    shared_cols = [
        col
        for col in cols_user
        if col in set_new and col != "ID"
    ]

    # -------------------------
    # Warning collection
    # -------------------------
    warnings = []

    # Columns present in the user file but missing from the new file
    for col in cols_user:
        if col not in set_new:
            warnings.append(
                {
                    "Header": col,
                    "Absent_in": args.tsv,
                }
            )

    # Columns present in the new file but missing from the user file
    for col in cols_new:
        if col not in set_user:
            warnings.append(
                {
                    "Header": col,
                    "Absent_in": args.tsv_user,
                }
            )

    # -------------------------
    # Write warning file
    # -------------------------
    if len(shared_cols) == 0:
        with open(args.warning, "w") as warning_file:
            warning_file.write(
                "No shared column was found between the input files.\n"
            )
    else:
        pd.DataFrame(
            warnings,
            columns=["Header", "Absent_in"],
        ).to_csv(
            args.warning,
            sep="\t",
            index=False,
        )

    # -------------------------
    # Build the complete column union
    # Keep the order from --tsv-user first,
    # then append columns only present in --tsv.
    # -------------------------
    all_cols = list(dict.fromkeys(cols_user + cols_new))

    # Ensure ID is always the first column
    all_cols = ["ID"] + [
        col for col in all_cols
        if col != "ID"
    ]

    # -------------------------
    # Reindex both datasets to the full schema
    # Missing columns are filled with empty strings.
    # -------------------------
    df_user = df_user.reindex(
        columns=all_cols,
        fill_value="0",
    )

    df_new = df_new.reindex(
        columns=all_cols,
        fill_value="0",
    )

    # -------------------------
    # Resolve IDs present in both files
    #
    # For a shared ID:
    # - If profiles are identical, keep the row once.
    # - If profiles differ, discard the --tsv-user row
    #   and keep only the --tsv row.
    #
    # The comparison is performed on all non-ID columns
    # after harmonizing both files to the same schema.
    # -------------------------
    shared_ids = set(df_user["ID"]) & set(df_new["ID"])

    if shared_ids:
        # Keep IDs that are present in both files out of the
        # user dataset. The corresponding row from --tsv will
        # always be retained below.
        df_user = df_user[
            ~df_user["ID"].isin(shared_ids)
        ]

    # -------------------------
    # Merge datasets
    #
    # IDs shared by both files are now represented only by
    # their row from --tsv.
    # IDs present in only one file are retained.
    # -------------------------
    merged = pd.concat(
        [df_user, df_new],
        ignore_index=True,
    )

    # -------------------------
    # Remove exact duplicate rows
    # -------------------------
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