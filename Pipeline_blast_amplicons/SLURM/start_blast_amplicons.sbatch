#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --gres=submitter:1
#SBATCH --output=/srv/net/cluqumngs/TMP_IAI/redmine/65198/%j.out
#SBATCH --error=/srv/net/cluqumngs/TMP_IAI/redmine/65198/%j.err
#SBATCH --comment="Blastn Amplicons pipeline, SLURM out in redmine/65198"
#SBATCH --job-name=Blast-amplicons
#SBATCH --partition=diag_iai
#SBATCH --nodelist=di6182su

set -e

# In local storage allocated by SLURM 
cd $SCRATCH

################################################################################
#                                                                              #
# start_blast_amplicons.sbatch version 1                                       #
#                                                                              #
# Aurelie PETICCA, last update: 2026-05                                        #
# Christophe GINEVRA                                                           #
#                                                                              #
# Aim: SLURM launch for Blastn Amplicons nextflow pipeline                     #
#                                                                              #
# Usage:  sbatch start_blast_amplicons.sbatch -s NGS-WEB_serie [options]       #
# WARNING : if launched via NGS-WEB, only the copy in /NGS-WEB/ is used        #
#                                                                              #
################################################################################

# Help screen
display_help() { 
    echo "Usage: sbatch start_blast_amplicons.sbatch -s NGS-WEB_serie [options]" >&2
    echo >&2
    echo >&2
    echo "   -s             [str]           NGS-WEB series, ends with date in format YYYYMMDD (required)" >&2
 	echo "   -c             [path]          nextflow config file to use, by default : 
                                                <nextflow_folder>/config/blast_amplicons.config" >&2
    echo >&2
    echo "   -h                             print this help message and exit" >&2
    echo >&2
}

usage() {
    echo "Usage: sbatch start_blast_amplicons.sbatch -s <YYYYMMDD> [options]" >&2
    echo >&2
    echo "Options:" >&2
    echo "  -s             NGS-WEB series (required)" >&2
    echo "  -c             Config file" >&2
    echo "  -h             Show help" >&2
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
# script_path="$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')"
# script_dir="$(dirname "$script_path")"
script_dir="/srv/scratch/DEPOT/pipelines/nextflow_pipeline_legio/Pipeline_blast_amplicons/SLURM"
config_file="${script_dir}/config/blast_amplicons.config"
pipeline_nf="${script_dir}/../workflow_blast_amplicons.nf"
dump_nf="${script_dir}/dump_params.nf"
nf_exec="${script_dir}/../../nextflow_25.10.4"

### Launch-specific configuration  
clusterOptions="--account=\"${SLURM_JOB_ACCOUNT}\" \
--partition=\"${SLURM_JOB_PARTITION}\" \
--nodelist=\"${SLURM_JOB_NODELIST}\""

## User values
### Check args presence
if [ $# -eq 0 ]; then
    echo "ERROR: NGS-WEB sequencing ID is required" >&2
    usage >&2
	echo >&2
    exit 1
fi

### Parsing loop
while getopts ":s:hc:" opt; do
    case ${opt} in
    s)
        # Name of the NGS-WEB series from which to retrieve the FASTQ files
        ngsweb_serie="${OPTARG}"
        # Extract last 8 characters
        ngsweb_id="${ngsweb_serie: -8}"
        ;;
    h)
        # Help option
        display_help
        exit 0
        ;;
    c)
        # Config file for dev
        config_file="${OPTARG}"
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
if [[ -z "${ngsweb_serie}" ]]; then
    echo "ERROR: NGS-WEB serie ID is required" >&2
    usage
	echo >&2
    exit 1
elif [[ ! "$ngsweb_id" =~ ^[0-9]{8}$ ]]; then
    echo "Error: NGS-WEB serie ID must end with 8 digits (got: $ngsweb_serie)" >&2
    usage
	echo >&2
    exit 1
fi

## Retrieving paths from nextflow config, need to start with "/"
mapfile -t paths < <(
    "${nf_exec}" run "${dump_nf}" \
        -c "${config_file}" \
        --suffix "${ngsweb_id}" \
        --clusterOptions "${clusterOptions}" \
        2>&1 | grep -E '^.*=' | cut -d= -f2-
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
        echo "output_folder : ${output_folder}" >&2
        echo "save_folder   : ${save_folder}" >&2
        echo "tmp_folder    : ${tmp_folder}" >&2
        echo "work_folder   : ${work_folder}" >&2
        echo "result_folder : ${result_folder}" >&2
        echo >&2
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
echo

## Copy raw data from input server to calculation engine
echo "--- SAVING INPUT DATA -----------------------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo

mkdir -p "${tmp_folder}"
chmod -R 777 "${tmp_folder}"
rsync -avQ --ignore-existing \
    --include='*.fastq' \
    --include='*.fastq.gz' \
    --exclude='*' \
    "${input_folder}/" "${tmp_folder}/"
echo

echo "--- FINISHED - to TMP FOLDER ----------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo

## Copy raw data from input server to storage server
mkdir -p "${save_folder}"
rsync -avQ --ignore-existing \
    --include='*.fastq' \
    --include='*.fastq.gz' \
    --exclude='*' \
    "${input_folder}/" "${save_folder}/"
echo

echo "--- FINISHED - to SAVE FOLDER ---------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo

## Start Blast analysis
echo "--- BLAST AMPLICONS ANALYSIS STARTING ------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo

### Notification
# echo "L'analyse BLAST AMPLICONS du run Legionella-Amplicons-${ngsweb_id} est en cours" \
# | mail -s "Analyse BLAST Legionella-Amplicons-${ngsweb_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

"${nf_exec}" \
    -C "${config_file}" \
    run "${pipeline_nf}" \
    --suffix "${ngsweb_id}" \
    -profile slurm \
    --clusterOptions "${clusterOptions}" \
    -with-trace "${SCRATCH}/nextflow_${ngsweb_id}.txt" \
    -with-report "${SCRATCH}/nextflow_${ngsweb_id}.html"

echo "--- FINISHED --------------------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo

# Copy results from calculation engine to storage server
echo "--- SAVING OUTPUT DATA and REMOVING TMP DATA ------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo

mkdir -p "${output_folder}"
rsync -avQ \
    --exclude='dev' \
    --exclude='work' \
    "$result_folder/" "$output_folder/"
echo
### NB : synchronising every files/folders, dev and work not needed

echo "--- FINISHED - to SAVE FOLDER ---------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo

## Remove results from calculation engine
echo "Deleting... ${work_folder}"
rm -r "${work_folder}"
echo "Deleting... ${tmp_folder}"
rm -r "${tmp_folder}"

echo "--- FINISHED - to DELETE --------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo

### Notification
# echo "L'analyse BLAST AMPLICONS du run Legionella-Amplicons-${ngsweb_id} est disponible ici : ${output_folder}" \
# | mail -s "Analyse BLAST Legionella-Amplicons-${ngsweb_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

echo "END -------------------------------------------------------------------------------------------------------------------"
