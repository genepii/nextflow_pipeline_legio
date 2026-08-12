#!/bin/bash

################################################################################
#                                                                              #
# start_strain_comparison.sh version 1                                         #
#                                                                              #
# Aurelie PETICCA, last update: 2026-08                                        #
# Christophe GINEVRA                                                           #
#                                                                              #
# Aim: Launch for Comparison of Lp strains nextflow pipeline                   #
#                                                                              #
# Usage:  start_strain_comparison.sh comparison_ID [options]                   #
#                                                                              #
################################################################################

# Stop the script as soon as a command returns a non-zero exit code
set -u

################################################################################
# Help screen
display_help() { 
 	echo "Usage: $0 -d comparison_ID [options]" >&2
	echo >&2
 	echo >&2
 	echo "   -d,--compare   [str]           COMPARISON_ID, required" >&2
 	echo "   -c,--config    [path]          nextflow config file to use, by default :
                                                <nextflow_folder>/config/strain_comparison.config" >&2
 	echo "   -w,--work      [path]          folder where all the temporary files will be written, by default :
                                                /srv/scratch/iai/bachcl/result/Legionella/Strain_comparison/{comparison_ID}/{analyse_id}_Strain-Comparison" >&2
 	echo "   -m,--tmp       [path]          temporary folder where the input files will be stored, by default : 
                                                /srv/scratch/iai/bachcl/Raw_fastq/Legionella/Strain_comparison/{comparison_ID}" >&2
 	echo "   -o,--output    [path]          folder where the final output files will be written, by default :
                                                /srv/net/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_results/Strain_comparison/{comparison_ID}/{analyse_id}_Strain-Comparison" >&2
 	echo "   -t,--adapters  [True/False]    remove Illumina adaptaters, by default : True" >&2
 	echo "   -z,--zoom      [str]           enable ReporTree zoom functionality on the samples analysed (analyse, by default), none, or sample(s) ID to zoom" >&2
 	echo "   -a,--meta      [path]          TSV file containing metadata for all the samples; must contain columns ID, ST, Year, Origin and Comparison; by default :
                                                /srv/net/cluqumngs/TMP_IAI/04_CNR_Legionella/Input_analysis_nextflow/Strains_{comparison_ID}.txt" >&2
    echo "   -s,--dbsnippy  [path]          path containing previous Snippy results and the reference used; see the manual for the hierarchy; by default :
                                                /srv/scratch/iai/bachcl/db/legio/Snippy" >&2
	echo >&2
 	echo "   -h, --help                     write this report and exit" >&2
    echo >&2
}

usage() {
    echo "Usage: $0 -d <comparison_id> [options]"
    echo ""
    echo "Options:"
    echo "  -d, --compare    Comparison ID (required)"
    echo "  -c, --config     Config file"
    echo "  -o, --output     Output folder"
    echo "  -m, --tmp        Temporary folder"
    echo "  -w, --work       Work folder"
    echo "  -t, --adapters   Remove Illumina adapters (True/False)"
    echo "  -z, --zoom       ReporTree zoom"
    echo "  -a, --meta       Metadata TSV file"
    echo "  -s, --dbsnippy   Snippy results/reference"
    echo "  -h, --help       Help"
}

################################################################################
# Configuration

## Variables init
analyse_id=""
comparison_id=""
config_file=""
output_folder_prefix=""
tmp_folder_prefix=""
work_folder_prefix=""
adapters=""
zoom=""
metadata_user=""
db_snippy=""

## Default values
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="${script_dir}/config/strain_comparison.config"

output_folder_prefix="/srv/net/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_results/Strain_comparison"
work_folder_prefix="/srv/scratch/iai/bachcl/result/Legionella/Strain_comparison"
tmp_folder_prefix="/srv/scratch/iai/bachcl/Raw_fastq/Legionella/Strain_comparison"
adapters="true"
zoom="analyse"
metadata_user_prefix="/srv/net/cluqumngs/TMP_IAI/04_CNR_Legionella/Input_analysis_nextflow/Strains_"
db_snippy="/srv/scratch/iai/bachcl/db/legio/Snippy"
analyse_id=$(date +%Y%m%d)

################################################################################
# User values

## Check args presence
if [ $# -eq 0 ]; then
    echo "ERROR: --compare is required"
    usage
    echo ""
    exit 1
fi

## Parsing loop
while [ $# -gt 0 ]; do
    case "$1" in

        -d|--compare)
            comparison_id="${2:?ERROR: missing value for --compare}"
            shift 2
            ;;

        -c|--config)
            config_file="${2:?ERROR: missing value for --config}"
            shift 2
            ;;

        -o|--output)
            output_folder_prefix="${2:?ERROR: missing value for --output}"
            shift 2
            ;;

        -m|--tmp)
            tmp_folder_prefix="${2:?ERROR: missing value for --tmp}"
            shift 2
            ;;

        -w|--work)
            work_folder_prefix="${2:?ERROR: missing value for --work}"
            shift 2
            ;;

        -t|--adapters)
            adapters="${2:?ERROR: missing value for --adapters}"
            shift 2
            ;;

        -a|--meta)
            metadata_user="${2:?ERROR: missing value for --meta}"
            shift 2
            ;;

        -z|--zoom)
            zoom="${2:?ERROR: missing value for --zoom}"

            if [[ "$zoom" =~ ^(none|analyse|[A-Za-z0-9-]{6,100})$ ]]; then
                :
            else
                echo "ERROR: invalid value for --zoom (allowed: analyse, none or sample ID list)" >&2
                exit 1
            fi
            shift 2
            ;;

        -s|--dbsnippy)
            db_snippy="${2:?ERROR: missing value for --dbsnippy}"
            shift 2
            ;;

        -h|--help)
            display_help
            exit 0
            ;;

        --)
            shift
            break
            ;;

        -*)
            echo "ERROR: Unknown option: $1" >&2
            usage
            exit 1
            ;;

        *)
            echo "ERROR: Unexpected argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

## Required argument check
if [[ -z "${comparison_id}" ]]; then
    echo "ERROR: --compare is required"
    usage
    echo ""
    exit 1
fi

## Generate folder names
output_folder="${output_folder_prefix}/${comparison_id}/${analyse_id}_Strain-Comparison"
tmp_folder="${tmp_folder_prefix}/${comparison_id}"
work_folder="${work_folder_prefix}/${comparison_id}/${analyse_id}_Strain-Comparison/work"
result_folder="${work_folder_prefix}/${comparison_id}/${analyse_id}_Strain-Comparison"


## Metadata
if [[ -z "${metadata_user}" ]]; then
    metadata_user="${metadata_user_prefix}${comparison_id}.txt"
fi

if [[ ! -f "${metadata_user}" ]]; then
    echo "ERROR: file not found: ${metadata_user}" >&2
    exit 1
fi

### Check if TSV
header=$(head -n1 "${metadata_user}")
[[ "${header}" == *$'\t'* ]] || {
    echo "ERROR: ${metadata_user} is not a TSV file." >&2
    exit 1
}

### Check mandatory columns
for col in ID ST Year Origin Comparison; do
    grep -qw "${col}" <<< "$(tr '\t' '\n' <<< "${header}")" || {
        echo "ERROR: Missing column '${col}' in ${metadata_user}" >&2
        exit 1
    }
done


################################################################################
# Start script
## Variables for launching Nextflow
pipeline_file="${script_dir}/workflow_strain_comparison.nf"
nf_exec="${script_dir}/../nextflow_25.10.4"

echo "START -----------------------------------------------------------------------------------------------------------------"
echo ""

## Copy raw data from input server to calculation engine
echo "--- STRAIN COMPARISON ANALYSIS STARTING -------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "The strain comparison analysis ${comparison_id} has started." \
# | mail -s "Strain comparison ${comparison_id}" your.email@example.com

if ! k5start -U -f /home/chu-lyon.fr/ginevrach/login.kt \
    -- "${nf_exec}" \
    -C "${config_file}" \
    run "${pipeline_file}" \
    --suffix "${comparison_id}" \
    --input "${tmp_folder}" \
    -w "${work_folder}" \
    -profile local \
    --result "${result_folder}" \
    --adapters "${adapters}" \
    --rep_zoom "${zoom}" \
    --metadata "${metadata_user}" \
    --snippy_dir "${db_snippy}" \
    -with-trace "${result_folder}/trace_${comparison_id}_${analyse_id}.txt" \
    -with-report "${result_folder}/report_${comparison_id}_${analyse_id}.html"
then
    LOG="error"
fi

echo "--- FINISHED ----------------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

# Copy results from calculation engine to storage server
echo "--- SAVING OUTPUT DATA -----------------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

mkdir -p "${output_folder}"
rsync -avQ \
    --exclude='dev' \
    --exclude='work' \
    "${result_folder}/" "${output_folder}/"
echo ""
### Synchronising every files/folders, dev and work not needed

## Add new Snippy results vs ST to Database
shopt -s nullglob
for sample_dir in "${result_folder}"/dev/1_Snippy-ST/ST*/*; do
    [[ -d "${sample_dir}" ]] || continue

    # Extract ST${comparison} and sample_id from path
    st_dir=$(basename "$(dirname "${sample_dir}")")
    sample_id=$(basename "${sample_dir}")

    destination="${db_snippy}/${st_dir}/${sample_id}"

    rsync -avQ \
        --exclude='*.bam' \
        --exclude='*.bam.bai' \
        "${sample_dir}/" \
        "${destination}/"
done

echo "--- FINISHED ----------------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""


## Remove results from calculation engine
echo "Deleting... ${work_folder}"
rm -r "${work_folder}"
rm -r "${result_folder}/dev/0-1_Trimmed"            # Warning: Delete Trimmed Reads for space
rm -f "${result_folder}"/dev/1_Snippy-*/*/*/*.bam*  # Warning: Delete Bam files for space
echo "Deleting... ${tmp_folder}"
rm -r "${tmp_folder}"

echo "--- FINISHED - to DELETE ----------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "The strain comparison analysis ${comparison_id} is available in ${output_folder}" \
# | mail -s "Strain comparison ${comparison_id}" your.email@example.com

echo "END -------------------------------------------------------------------------------------------------------------------"