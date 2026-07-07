#!/usr/bin/env python3
# ------------------------------------------------------------
# assembly_mlst_summary_table.py
#
# Merge multiple pipeline summary tables into a single report.
#
# Features:
# - Outer join on Sample_ID
# - Depth computation (Total_bases / genome size)
# - Controlled column ordering
# - Guaranteed presence of all protected columns
# - Conditional removal of fully empty columns (non-protected only)
# - Insertion of MompS_* block between neuA and lag
# ------------------------------------------------------------

import argparse
from functools import reduce
import pandas as pd


# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------
def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge summary tables into a single report"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file"
    )

    parser.add_argument(
        "files",
        nargs="+",
        help="Input summary TSV files"
    )

    return parser.parse_args()


# ------------------------------------------------------------
# Load TSV safely
# ------------------------------------------------------------
def load_table(path):
    """
    Load a TSV file as a pandas DataFrame and ensure Sample_ID exists.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)

    if "Sample_ID" not in df.columns:
        raise ValueError(f"Missing Sample_ID in {path}")

    return df


# ------------------------------------------------------------
# Ensure mandatory columns exist
# ------------------------------------------------------------
def ensure_base_columns(df, base_columns):
    """
    Ensure that all columns in base_columns exist in the DataFrame.
    Missing columns are created and filled with 'NA'.
    """
    for col in base_columns:
        if col not in df.columns:
            df[col] = "NA"
    return df


# ------------------------------------------------------------
# Drop columns fully NA unless protected
# ------------------------------------------------------------
def drop_all_na_columns(df, protected_columns):
    """
    Remove columns fully empty (NaN or 'NA'),
    except columns explicitly protected.
    """

    keep = []

    for col in df.columns:

        # Always keep protected columns
        if col in protected_columns:
            keep.append(col)
            continue

        series = df[col]

        is_all_na = series.isna().all() or series.astype(str).eq("NA").all()

        if not is_all_na:
            keep.append(col)

    return df[keep]


# ------------------------------------------------------------
# Column ordering logic
# ------------------------------------------------------------
def reorder_columns(df):
    """
    Reorder columns according to a fixed schema and insert MompS block.
    """

    base_order = [
        "Sample_ID",
        "Total_bases",
        "Total_length",
        "Depth",
        "Number_of_contigs",
        "GC",
        "N50",
        "auN",
        "L90",
        "Legionella_pneumophila_percent",
        "Legionella_spp_percent",
        "Kraken2_results",
        "FastANI_strain",
        "FastANI_value",
        "ST",
        "flaA",
        "pilE",
        "asd",
        "mip",
        "mompS",
        "proA",
        "neuA",
        "lag",
        "lpeA",
        "lpeB",
        "AMR_Nb_Mutated_Genes",
        "AMR_Nb_with_Impact",
        "AMR_Nb_Non_Coding"
    ]

    momps_block = [
        "MompS_ST",
        "MompS_flaA",
        "MompS_pilE",
        "MompS_asd",
        "MompS_mip",
        "MompS_mompS",
        "MompS_proA",
        "MompS_neuA"
    ]

    # Ensure all base columns exist before ordering
    df = ensure_base_columns(df, base_order)

    # Keep only columns that exist
    base_existing = [c for c in base_order if c in df.columns]
    momps_existing = [c for c in momps_block if c in df.columns]

    # Locate insertion points
    idx_neuA = base_existing.index("neuA") + 1
    idx_lag = base_existing.index("lag")

    before_momps = base_existing[:idx_neuA]
    after_momps = base_existing[idx_lag:]

    used = set(before_momps + after_momps + momps_existing)
    remaining = [c for c in df.columns if c not in used]

    return df[before_momps + momps_existing + after_momps + remaining]


# ------------------------------------------------------------
# Main pipeline
# ------------------------------------------------------------
def main():
    """Main execution pipeline."""
    args = parse_args()

    # Load all input tables
    dfs = [load_table(f) for f in args.files]

    # Merge on Sample_ID (outer join)
    merged = reduce(
        lambda l, r: pd.merge(l, r, on="Sample_ID", how="outer"),
        dfs
    )

    # --------------------------------------------------------
    # Genome size detection
    # --------------------------------------------------------
    genome_col = next(
        (c for c in ["Total_length", "Genome_Size", "Total_length_bp"]
         if c in merged.columns),
        None
    )

    merged["Genome_Size"] = pd.to_numeric(
        merged[genome_col], errors="coerce"
    ) if genome_col else pd.NA

    merged["Total_bases"] = pd.to_numeric(
        merged.get("Total_bases", pd.Series()), errors="coerce"
    )

    # Compute sequencing depth
    merged["Depth"] = merged["Total_bases"] / merged["Genome_Size"]

    merged.drop(columns=["Genome_Size"], inplace=True, errors="ignore")

    # --------------------------------------------------------
    # Define protected columns (never removed)
    # Guarantee all base columns exist BEFORE cleaning
    # --------------------------------------------------------
    protected_columns = [
        "Sample_ID",
        "Total_bases",
        "Total_length",
        "Depth",
        "Number_of_contigs",
        "GC",
        "N50",
        "auN",
        "L90",
        "Legionella_pneumophila_percent",
        "Legionella_spp_percent",
        "Kraken2_results",
        "FastANI_strain",
        "FastANI_value",
        "ST",
        "flaA",
        "pilE",
        "asd",
        "mip",
        "mompS",
        "proA",
        "neuA",
        "lag",
        "lpeA",
        "lpeB",
        "AMR_Nb_Mutated_Genes",
        "AMR_Nb_with_Impact",
        "AMR_Nb_Non_Coding"
    ]

    merged = ensure_base_columns(merged, protected_columns)

    # Remove fully empty columns (except protected ones)
    merged = drop_all_na_columns(merged, protected_columns)

    # Fill remaining missing values for consistency
    merged = merged.fillna("NA")

    # Sort rows
    merged = merged.sort_values("Sample_ID")

    # Final column ordering
    merged = reorder_columns(merged)

    # Write output
    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()