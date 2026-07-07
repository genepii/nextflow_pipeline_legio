"""
VCF processing script (HPC / Nextflow compatible)

Purpose:
- Recalculate allele frequency (AO/DP)
- Filter by variant type and frequency bounds
- Rebuild INFO field in VCF-like output

Input: VCF file
Output: processed VCF-like file

Constraints:
- Python 2.7
- Stateless, reproducible
- Memory-efficient streaming
"""

from __future__ import print_function
import sys
import getopt


def usage():
    # Display command-line usage instructions
    print(
        "usage: assembly_mlst_recalc_fbvcf.py "
        "-i input.vcf -o output.vcf "
        "-t snp,del,ins,mnp,complex "
        "-n 0.0 -x 1.0"
    )


def parse_args(argv):
    # Parse and validate CLI arguments
    cfg = {
        "input_file": None,
        "output_file": None,
        "call_types": ["snp", "del", "ins", "mnp", "complex"],
        "min_freq": 0.0,
        "max_freq": 1.0
    }

    try:
        opts, _ = getopt.getopt(
            argv,
            "hi:o:t:n:x:",
            ["help", "input", "output", "calltype", "minfreq", "maxfreq"]
        )

        for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit(0)

            elif opt in ("-i", "--input"):
                cfg["input_file"] = arg

            elif opt in ("-o", "--output"):
                cfg["output_file"] = arg

            elif opt in ("-t", "--calltype"):
                # Convert comma-separated list into Python list
                cfg["call_types"] = arg.split(",")

            elif opt in ("-n", "--minfreq"):
                cfg["min_freq"] = float(arg)

            elif opt in ("-x", "--maxfreq"):
                cfg["max_freq"] = float(arg)

        # Ensure required arguments are provided
        if not cfg["input_file"] or not cfg["output_file"]:
            usage()
            sys.exit(1)

        return cfg

    except getopt.GetoptError:
        # Handle invalid CLI arguments
        usage()
        sys.exit(2)


def truncate_float(value, decimals):
    # Truncate float without rounding to fixed decimals
    s = str(value)

    # Handle scientific notation cases
    if "e" in s or "E" in s:
        return float(("{0:." + str(decimals) + "f}").format(value))

    # Handle integer values
    if "." not in s:
        return float(s)

    # Manual truncation of decimal part
    i, _, d = s.partition(".")
    return float(i + "." + (d + "0" * decimals)[:decimals])


def parse_info_field(info_items):
    # Parse INFO column into key/value arrays
    keys = []
    values = []

    for item in info_items:
        k, v = item.split("=", 1)
        keys.append(k)
        values.append(v)

    return keys, values


def process_vcf(cfg):
    # Precompute frequency thresholds (HPC optimization)
    min_f = truncate_float(cfg["min_freq"], 3)
    max_f = truncate_float(cfg["max_freq"], 3)

    header_lines = []
    output_lines = []

    # Stream input file (memory efficient for large VCF)
    with open(cfg["input_file"], "r") as f:
        for line in f:
            line = line.rstrip("\n")

            # Separate header from data
            if line.startswith("#"):
                header_lines.append(line)
                continue

            fields = line.split("\t")

            # Parse INFO field
            info_items = fields[7].split(";")
            keys, values = parse_info_field(info_items)

            # Strict validation of expected INFO structure
            if (keys[3] != "AF" or keys[5] != "AO" or keys[7] != "DP" or
                keys[22] != "PAO" or keys[25] != "PRO" or keys[28] != "RO" or
                keys[40] != "TYPE"):
                sys.exit("Fields wrongly formatted")

            # Allele frequency computation
            af = float(int(values[5])) / float(values[7])
            af_trunc = truncate_float(af, 3)

            # Apply variant type logic
            if values[40] not in cfg["call_types"]:
                values[3] = str(af_trunc)
            else:
                # Clamp frequency within bounds
                if af_trunc < min_f:
                    values[3] = str(min_f)
                elif af_trunc > max_f:
                    values[3] = str(max_f)
                else:
                    values[3] = str(af_trunc)

            # Rebuild INFO field
            new_info = []
            for i in range(len(keys)):
                new_info.append(keys[i] + "=" + values[i])

            fields[7] = ";".join(new_info)

            # Store processed record
            output_lines.append("\t".join(fields))

    return header_lines, output_lines


def write_output(path, header, records):
    # Write final output file
    with open(path, "w") as out:
        for h in header:
            out.write(h + "\n")
        for r in records:
            out.write(r + "\n")


def main(argv):
    # Main execution pipeline
    cfg = parse_args(argv)
    header, records = process_vcf(cfg)
    write_output(cfg["output_file"], header, records)


if __name__ == "__main__":
    main(sys.argv[1:])