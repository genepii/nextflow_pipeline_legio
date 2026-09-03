#!/usr/bin/env python3

"""
Recursively search FASTA files using configurable naming rules.

Search rules:
    - FASTA files are searched recursively from all matching directories.
    - Input directories may contain shell-style wildcards such as '*'.
    - Multiple wildcards are supported in a single path.
    - Paths without wildcards are also supported.
    - Files located inside directories named "scaffolds" are ignored.
    - The newest matching FASTA file is selected for each sample.

Examples:
    /data/genomes
    /data/genomes/*/contigs
    /data/genomes/*/contigs/*/500b

Outputs:
    - output:
        Metadata TSV containing paths to copied FASTA files.
    - RealPath_<output>:
        Metadata TSV containing original FASTA file paths.
"""

from pathlib import Path
import argparse
import glob
import shutil

import pandas as pd


def parse_arguments():
    """Parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Find and copy FASTA files"
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Input metadata TSV file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output metadata TSV file"
    )

    parser.add_argument(
        "--copy",
        required=True,
        help="Directory where selected FASTA files are copied"
    )

    parser.add_argument(
        "--directories",
        nargs="+",
        required=True,
        help=(
            "FASTA search directories or glob patterns. "
            "For example: /data/genomes or /data/genomes/*/contigs"
        )
    )

    parser.add_argument(
        "--prefix",
        nargs="+",
        default=["false"],
        help="FASTA filename prefixes"
    )

    parser.add_argument(
        "--suffix",
        nargs="+",
        default=["false"],
        help="FASTA filename suffixes"
    )

    parser.add_argument(
        "--extensions",
        nargs="+",
        default=["false"],
        help="FASTA filename extensions"
    )

    return parser.parse_args()


def build_sample_patterns(
    sample_id,
    prefixes,
    suffixes,
    extensions
):
    """
    Generate possible FASTA filenames.

    Filename format:
        {prefix}{ID}{suffix}{extension}
    """

    if prefixes == ["false"]:
        prefixes = [""]

    if suffixes == ["false"]:
        suffixes = [""]

    if extensions == ["false"]:
        extensions = [""]

    patterns = []

    for prefix in prefixes:
        for suffix in suffixes:
            for extension in extensions:
                patterns.append(
                    f"{prefix}{sample_id}{suffix}{extension}"
                )

    return patterns


def expand_directories(directory_patterns):
    """
    Expand directory paths containing glob patterns.

    Both regular paths and paths containing wildcards are supported.

    Examples:
        /data/genomes
        /data/genomes/*/contigs
        /data/genomes/*/contigs/*/500b

    Returns:
        A list of existing directories matching the input patterns.
    """

    directories = []

    for pattern in directory_patterns:

        # glob() returns the path itself when there is no wildcard,
        # provided that the path exists.
        matches = glob.glob(
            pattern,
            recursive=True
        )

        for match in matches:

            path = Path(match)

            if path.is_dir():
                directories.append(path)

    # Remove duplicates while preserving the original order.
    unique_directories = list(
        dict.fromkeys(directories)
    )

    return unique_directories


def find_fasta_files(
    directories,
    sample_patterns
):
    """
    Search FASTA files recursively.

    Input directories may contain glob patterns. All directories
    matching these patterns are searched recursively.

    Files located inside any directory named "scaffolds"
    are ignored.

    If several FASTA files match the same sample, the newest
    file based on modification time is selected.

    Returns:
        Dictionary mapping sample IDs to selected FASTA paths.
    """

    index = {}

    pattern_lookup = {}

    for sample, patterns in sample_patterns.items():
        for pattern in patterns:
            pattern_lookup[pattern] = sample

    # Expand all directory patterns before starting the recursive search.
    roots = expand_directories(directories)

    for root in roots:

        for file in root.rglob("*"):

            if not file.is_file():
                continue

            # Ignore files located inside any "scaffolds" directory.
            if "scaffolds" in file.parts:
                continue

            sample = pattern_lookup.get(file.name)

            if sample is None:
                continue

            mtime = file.stat().st_mtime

            # Keep the newest matching FASTA file for each sample.
            if (
                sample not in index
                or mtime > index[sample][0]
            ):
                index[sample] = (
                    mtime,
                    file
                )

    return {
        sample: path
        for sample, (_, path) in index.items()
    }


def copy_fasta_files(
    fasta_files,
    output_directory
):
    """
    Copy selected FASTA files to the output directory.

    Returns:
        copied_files:
            Paths of copied FASTA files.
        original_files:
            Original FASTA file paths.
    """

    output_directory = Path(output_directory)

    output_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    copied_files = {}
    original_files = {}

    for sample, fasta in fasta_files.items():

        destination = output_directory / fasta.name

        shutil.copy2(
            fasta,
            destination
        )

        copied_files[sample] = str(
            destination.resolve()
        )

        original_files[sample] = str(
            fasta.resolve()
        )

    return copied_files, original_files


def main():

    args = parse_arguments()

    metadata = pd.read_csv(
        args.metadata,
        sep="\t",
        dtype=str
    ).fillna("NA")

    sample_patterns = {}

    for sample_id in metadata["ID"]:

        sample_patterns[sample_id] = build_sample_patterns(
            sample_id,
            args.prefix,
            args.suffix,
            args.extensions
        )

    fasta_files = find_fasta_files(
        args.directories,
        sample_patterns
    )

    fasta_index, realpath_index = copy_fasta_files(
        fasta_files,
        args.copy
    )

    metadata["FASTA"] = metadata["ID"].map(
        lambda sample: fasta_index.get(
            sample,
            "NA"
        )
    )

    metadata.to_csv(
        args.output,
        sep="\t",
        index=False
    )

    realpath_metadata = metadata.copy()

    realpath_metadata["FASTA"] = realpath_metadata["ID"].map(
        lambda sample: realpath_index.get(
            sample,
            "NA"
        )
    )

    output_path = Path(args.output)

    realpath_output = output_path.with_name(
        f"RealPath_{output_path.name}"
    )

    realpath_metadata.to_csv(
        realpath_output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
