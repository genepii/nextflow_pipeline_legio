#!/bin/bash
#SBATCH --time 24:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --gres submitter:1
#SBATCH --output /srv/scratch/peticcaau/SLURM/%j.out
#SBATCH --error /srv/scratch/peticcaau/SLURM/%j.err
#SBATCH --partition diag_iai 
#SBATCH --nodelist di6182su

set -e

################################################################################
#                                                                              #
# start_elgato_sbt.sh version 1                                                #
#                                                                              #
# Aurelie PETICCA, last update: 2026-05                                        #
# Christophe GINEVRA                                                           #
#                                                                              #
# Aim: SLURM launch for El Gato Nested SBT nextflow pipeline                   #
#                                                                              #
# Usage:  sbatch start_elgato_sbt.sh -s NGS-WEB_ID [options]                   #
#                                                                              #
################################################################################

# Help screen
display_help() { 
    echo "Usage: sbatch start_elgato_sbt.sh -s NGS-WEB_ID [options]" >&2
    echo >&2
    echo >&2
    echo "   -s             [str]           NGS-WEB series, locally a date in format YYYYMMDD (required)" >&2
    echo "   -r             [str]           Slurm job id to rerun (optional)" >&2
    echo >&2
    echo "   -h, --help                     print this help message and exit" >&2
    echo >&2
}

usage() {
    echo "Usage: sbatch start_elgato_sbt.sh -s <YYYYMMDD> [options]"
    echo ""
    echo "Options:"
    echo "  -s             NGS-WEB series (required)"
    echo "  -r             Slurm job id"
    echo "  -h, --help     Show help"
}

################################################################################
# Configuration
## Variables init
ngsweb_id=""
script_dir=""
config_file=""
pipeline_nf=""
nf_exec=""

## Nextflow default values
script_path="$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')"
script_dir="$(dirname "$script_path")"
config_file="${script_dir}/../config/elgato_sbt.config"
pipeline_nf="${script_dir}/../../workflow_elgato_sbt.nf"
dump_nf="${script_dir}/../dump_params.nf"
nf_exec="${script_dir}/../../nextflow_25.10.4"

## User values
### Check args presence
if [ $# -eq 0 ]; then
    echo "ERROR: NGS-WEB sequencing ID is required"
    usage
	echo ""
    exit 1
fi

### Parsing loop
while getopts ":r:s:h" opt; do
    case ${opt} in
    r)
        # [optional] SLURM job number of a failed run to be resumed
        slurm_job="${OPTARG}"
        ;;
    s)
        # Name of the NGS-WEB series from which to retrieve the FASTQ files
        ngsweb_id="${OPTARG}"
        ;;
    h)
        # Help option
        display_help
        exit 0
        ;;
    :)
        echo "ERROR: Option -${OPTARG} requires an argument" 1>&2
        exit 1
        ;;
    \?)
        echo "ERROR: Unknown option -${OPTARG}" 1>&2
        usage
        exit 1
        ;;
    esac
done

### Required argument check
if [[ -z "${ngsweb_id}" ]]; then
    echo "ERROR: NGS-WEB sequencing ID is required"
    usage
	echo ""
    exit 1
fi

## Retrieving paths from nextflow config, need to start with "/"
mapfile -t paths < <(
    "${nf_exec}" run "${dump_nf}" \
        -c "${config_file}" \
        --suffix "${ngsweb_id}" \
        2>&1 | grep '^/'
)
input_folder="${paths[0]}"
output_folder="${paths[1]}"
save_folder="${paths[2]}"
tmp_folder="${paths[3]}"
work_folder="${paths[4]}"
result_folder="${paths[5]}"

# If path empty error
for d in "${output_folder}" "${save_folder}" "${tmp_folder}" "${work_folder}" "${result_folder}"; do
    if [[ -z "$d" ]]; then
        echo "ERROR: One of the required paths is empty" >&2
        exit 1
    fi
    mkdir -p "$d"
done

# Check input_folder contains FASTQ files
if [[ -z "${input_folder}" ]]; then
    echo "ERROR: No input_folder given" >&2
    exit 1
fi

shopt -s nullglob

fastq_files=( "${input_folder}"/*.fastq "${input_folder}"/*.fastq.gz )

if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "ERROR: No FASTQ or FASTQ.GZ files found in input_folder: ${input_folder}" >&2
    exit 1
fi


################################################################################
# start script
echo "START -----------------------------------------------------------------------------------------------------------------"
echo ""

## Copy raw data from input server to calculation engine
echo "--- SAVING INPUT DATA -----------------------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

mkdir -p "${tmp_folder}"
chmod -R 777 "${tmp_folder}"
rsync -avQ --ignore-existing \
    --include='*.fastq' \
    --include='*.fastq.gz' \
    --exclude='*' \
    "${input_folder}/" "${tmp_folder}/"
echo ""

echo "--- FINISHED - to TMP FOLDER ----------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

## Copy raw data from input server to storage server
mkdir -p "${save_folder}"
rsync -avQ --ignore-existing \
    --include='*.fastq' \
    --include='*.fastq.gz' \
    --exclude='*' \
    "${input_folder}/" "${save_folder}/"
echo ""

echo "--- FINISHED - to SAVE FOLDER ---------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

## Start Nested analysis
echo "--- EL GATO NESTED SBT ANALYSIS STARTING ------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "L'analyse EL GATO NESTED SBT du run Legionella-Nested-SBT-${ngsweb_id} est en cours" \
# | mail -s "Analyse EL GATO NESTED SBT Legionella-Nested-SBT-${ngsweb_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

"${nf_exec}" \
    -C "${config_file}" \
    run "${pipeline_nf}" \
    --suffix "${ngsweb_id}" \
    -profile slurm \
    -with-trace "${result_folder}/trace_${ngsweb_id}.txt" \
    -with-report "${result_folder}/report_${ngsweb_id}.html"

echo "--- FINISHED --------------------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

# Copy results from calculation engine to storage server
echo "--- SAVING OUTPUT DATA and REMOVING TMP DATA ------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

mkdir -p "${output_folder}"
rsync -avQ \
    --exclude='dev' \
    --exclude='work' \
    "$result_folder/" "$output_folder/"
echo ""
### NB : synchronising every files/folders, dev and work not needed

echo "--- FINISHED - to SAVE FOLDER ---------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

## Remove results from calculation engine
echo "Deleting... ${work_folder}"
rm -r "${work_folder}"
echo "Deleting... ${tmp_folder}"
rm -r "${tmp_folder}"

echo "--- FINISHED - to DELETE --------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "L'analyse EL GATO NESTED SBT du run Legionella-Nested-SBT-${ngsweb_id} est disponible ici : ${output_folder}" \
# | mail -s "Analyse EL GATO NESTED SBT Legionella-Nested-SBT-${ngsweb_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

echo "END -------------------------------------------------------------------------------------------------------------------"
