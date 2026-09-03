#!/usr/bin/env python3

"""
Recursively search FASTQ files using configurable naming rules.

Select the newest R1 and R2 files for each sample and copy them
to the output directory.

Input directories can be:
    - regular paths without wildcards
    - paths containing one or more '*' wildcards

Examples:
    /data/fastq
    /data/fastq/*/contigs
    /data/fastq/*/contigs/*/500b

Generate two metadata files:
    - output:
        paths to copied FASTQ files
    - RealPath_<output>:
        original FASTQ file paths
"""

from pathlib import Path
import argparse
import fnmatch
import glob
import re
import shutil

import pandas as pd


def parse_arguments():
    """Parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Find and copy FASTQ R1/R2 files"
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
        help="Directory where selected FASTQ files are copied"
    )

    parser.add_argument(
        "--directories",
        nargs="+",
        required=True,
        help="FASTQ search directories or paths containing '*'"
    )

    parser.add_argument(
        "--prefix",
        nargs="+",
        default=["false"],
        help="FASTQ filename prefixes"
    )

    parser.add_argument(
        "--suffix",
        nargs="+",
        default=["false"],
        help="FASTQ filename suffix patterns"
    )

    parser.add_argument(
        "--extensions",
        nargs="+",
        default=["false"],
        help="FASTQ filename extension patterns"
    )

    return parser.parse_args()


def expand_R_pattern(pattern):
    """
    Expand R{1,2} notation into R1 and R2 patterns.
    """

    if "R{1,2}" not in pattern:
        return [pattern]

    return [
        pattern.replace("R{1,2}", "R1"),
        pattern.replace("R{1,2}", "R2")
    ]


def build_sample_patterns(
    sample_id,
    prefixes,
    suffixes,
    extensions
):
    """
    Generate possible FASTQ filenames for one sample.

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
            for expanded_suffix in expand_R_pattern(suffix):
                for extension in extensions:
                    patterns.append(
                        f"{prefix}{sample_id}{expanded_suffix}{extension}"
                    )

    return patterns


def detect_read(filename):
    """Identify FASTQ read direction."""

    if re.search(r"R1", filename):
        return "R1"

    if re.search(r"R2", filename):
        return "R2"

    if re.search(r"_1(\D|$)", filename):
        return "R1"

    if re.search(r"_2(\D|$)", filename):
        return "R2"

    return None


def expand_directories(directories):
    """
    Expand directory paths containing '*' wildcards.

    Paths without wildcards are also accepted.

    Examples:
        /data/fastq
        /data/fastq/*/contigs
        /data/fastq/*/contigs/*/500b

    Returns:
        List of existing directories matching the input paths.
    """

    expanded_directories = []

    for directory in directories:

        # Expand '*' if present.
        # If there is no wildcard, glob() simply returns the
        # existing path.
        matches = glob.glob(directory)

        for match in matches:

            path = Path(match)

            if path.is_dir():
                expanded_directories.append(path)

    # Remove duplicate directories while preserving their order.
    return list(
        dict.fromkeys(expanded_directories)
    )


def find_fastq_files(
    directories,
    sample_patterns
):
    """
    Search FASTQ files recursively.

    Input directories can contain '*' wildcards. Every directory
    matching the provided patterns is searched recursively.

    Keep the newest matching R1/R2 file for each sample.
    """

    index = {}

    pattern_lookup = {}

    for sample, patterns in sample_patterns.items():
        for pattern in patterns:
            pattern_lookup[pattern] = sample

    # Expand directory patterns before starting the recursive search.
    search_directories = expand_directories(
        directories
    )

    for root in search_directories:

        for file in root.rglob("*"):

            if not file.is_file():
                continue

            sample = None

            # Match filenames using shell-style wildcards.
            for pattern, sample_id in pattern_lookup.items():

                if fnmatch.fnmatch(
                    file.name,
                    pattern
                ):
                    sample = sample_id
                    break

            if sample is None:
                continue

            read = detect_read(file.name)

            if read is None:
                continue

            key = (
                sample,
                read
            )

            mtime = file.stat().st_mtime

            # Keep the newest file for each sample/read pair.
            if (
                key not in index
                or mtime > index[key][0]
            ):
                index[key] = (
                    mtime,
                    file
                )

    return {
        key: path
        for key, (_, path) in index.items()
    }


def copy_fastq_files(
    fastq_files,
    output_directory
):
    """
    Copy selected FASTQ files to output directory.

    Return:
        - copied file paths indexed by sample/read
        - original file paths indexed by sample/read
    """

    output_directory = Path(output_directory)

    output_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    copied_files = {}
    original_files = {}

    for key, fastq in fastq_files.items():

        destination = output_directory / fastq.name

        shutil.copy2(
            fastq,
            destination
        )

        copied_files[key] = str(
            destination.resolve()
        )

        original_files[key] = str(
            fastq.resolve()
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

    fastq_files = find_fastq_files(
        args.directories,
        sample_patterns
    )

    fastq_index, realpath_index = copy_fastq_files(
        fastq_files,
        args.copy
    )

    reads_1 = []
    reads_2 = []

    real_reads_1 = []
    real_reads_2 = []

    for sample_id in metadata["ID"]:

        reads_1.append(
            fastq_index.get(
                (sample_id, "R1"),
                "NA"
            )
        )

        reads_2.append(
            fastq_index.get(
                (sample_id, "R2"),
                "NA"
            )
        )

        real_reads_1.append(
            realpath_index.get(
                (sample_id, "R1"),
                "NA"
            )
        )

        real_reads_2.append(
            realpath_index.get(
                (sample_id, "R2"),
                "NA"
            )
        )

    metadata["READS_1"] = reads_1
    metadata["READS_2"] = reads_2

    metadata.to_csv(
        args.output,
        sep="\t",
        index=False
    )

    realpath_metadata = metadata.copy()

    realpath_metadata["READS_1"] = real_reads_1
    realpath_metadata["READS_2"] = real_reads_2

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
