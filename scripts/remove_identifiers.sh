#!/usr/bin/env bash

set -Eeuo pipefail

# Display script usage.
usage() {
    cat <<'EOF'
Usage:
    remove_identifiers.sh -i|--list FILE [-d|--db DIRECTORY]
    remove_identifiers.sh -h|--help

Arguments:
    -i, --list FILE       File containing identifiers (one per line).
    -d, --db DIRECTORY    Directory containing files to process.
                          Defaults to the current working directory.
    -h, --help            Display this help message.

Description:
    Removes every line containing an exact identifier followed by
    whitespace or the end of the line.

    Only files directly inside the target directory are processed.
    Subdirectories are ignored.

    Identifiers may only contain letters, digits and hyphens.
EOF
}

# Print an error message and exit.
die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

# Parse command-line arguments.
list_file=""
db_dir="."

while (($# > 0)); do
    case "$1" in
        -i|--list)
            (($# >= 2)) || die "Missing value for $1"
            list_file="$2"
            shift 2
            ;;

        -d|--db)
            (($# >= 2)) || die "Missing value for $1"
            db_dir="$2"
            shift 2
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            die "Unknown argument: $1. Use --help for usage."
            ;;
    esac
done

# Validate input arguments.
[[ -n "$list_file" ]] || die "Missing required argument: -i|--list"
[[ -f "$list_file" && -r "$list_file" ]] || die "List file does not exist or is not readable: $list_file"
[[ -d "$db_dir" && -r "$db_dir" && -w "$db_dir" ]] || die "Database directory does not exist or is not accessible: $db_dir"

# Resolve absolute paths.
list_file=$(realpath "$list_file")
db_dir=$(realpath "$db_dir")

# Create temporary directory.
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

tmp_list="$tmp_dir/identifiers.txt"
tmp_found="$tmp_dir/found.txt"

# Load identifiers, validate characters and remove duplicates.
declare -A identifiers=()

while IFS= read -r identifier || [[ -n "$identifier" ]]; do

    # Ignore empty lines.
    [[ -n "$identifier" ]] || continue

    # Validate identifier characters.
    [[ "$identifier" =~ ^[a-zA-Z0-9-]+$ ]] ||
        die "Invalid identifier in list file: '$identifier'. Only letters, digits and hyphens are allowed."

    identifiers["$identifier"]=1

done < "$list_file"

nb_asked=${#identifiers[@]}

((nb_asked > 0)) || die "No valid identifiers found in the list."

# Write unique identifiers to temporary file.
printf '%s\n' "${!identifiers[@]}" > "$tmp_list"

# Process files directly inside the database directory.
nb_files=0
nb_modified_files=0
nb_deleted=0

while IFS= read -r -d '' file; do

    [[ -f "$file" ]] || continue
    [[ -r "$file" && -w "$file" ]] || die "File is not readable or writable: $file"

    ((nb_files += 1))

    tmp_output="$tmp_dir/output"
    tmp_count="$tmp_dir/count"
    tmp_file_found="$tmp_dir/found_current"

    # Reset temporary files for the current input file.
    : > "$tmp_file_found"

    # Count matching identifiers and remove matching lines.
    if ! awk \
        -v found_file="$tmp_file_found" \
        -v count_file="$tmp_count" '
        NR == FNR {
            identifiers[$0] = 1
            next
        }
        {
            line = $0
            delete_line = 0

            for (id in identifiers) {
                if (line ~ ("(^|[[:space:]])" id "([[:space:]]|$)")) {
                    delete_line = 1
                    found[id] = 1
                }
            }

            if (delete_line) {
                deleted++
            } else {
                print
            }
        }
        END {
            for (id in found) {
                print id >> found_file
            }

            print deleted + 0 > count_file
        }
    ' "$tmp_list" "$file" > "$tmp_output"; then
        die "Failed to process file: $file"
    fi

    # Read the number of deleted lines.
    file_deleted=$(<"$tmp_count")

    # Collect identifiers found in the current file.
    cat "$tmp_file_found" >> "$tmp_found"

    # Replace original file only if matching lines were found.
    if ((file_deleted > 0)); then
        cat "$tmp_output" > "$file"
        ((nb_modified_files += 1))
        ((nb_deleted += file_deleted))
    fi

done < <(find "$db_dir" -maxdepth 1 -type f -print0)

# Remove duplicate found identifiers.
if [[ -s "$tmp_found" ]]; then
    sort -u "$tmp_found" -o "$tmp_found"
fi

# Count unique found identifiers.
nb_found=0

if [[ -s "$tmp_found" ]]; then
    nb_found=$(wc -l < "$tmp_found")
    nb_found=${nb_found//[[:space:]]/}
fi

# Display identifier status in alphabetical order.
printf '\n'

while IFS= read -r identifier; do

    if grep -Fxq -- "$identifier" "$tmp_found"; then
        printf '%-20s ......... FOUND\n' "$identifier"
    else
        printf '%-20s ......... NOT FOUND\n' "$identifier"
    fi

done < <(sort "$tmp_list")

# Display summary.
printf '\n========== SUMMARY ==========\n'
printf 'Identifiers requested : %d\n' "$nb_asked"
printf 'Identifiers found     : %d\n' "$nb_found"
printf 'Lines deleted         : %d\n' "$nb_deleted"
printf 'Files processed       : %d\n' "$nb_files"
printf 'Files modified        : %d\n' "$nb_modified_files"
printf '=============================\n'
