#!/bin/bash

################################################################################
#                                                                              #
# start_qiime2.sh version 1                                                    #
#                                                                              #
# Aurelie PETICCA, last update: 2026-04                                        #
#                                                                              #
# Aim: Launch for Qiime2 nextflow pipeline                                     #
#                                                                              #
# Usage:  start_qiime2.sh sequencing_ID [options]                              #
#                                                                              #
################################################################################

# Stop the script as soon as a command returns a non-zero exit code
set -u

################################################################################
# Help screen
display_help() { 
 	echo "Usage: $0 -d sequencing_ID [options]">&2
	echo >&2
 	echo >&2
 	echo "   -d,--seq_id    [str]           SEQUENCING_ID, locally a date in format YYYYMMDD, Mandatory" >&2
 	echo "   -i,--input     [path]          folder containing the sequencing data, by default : 
                                                /srv/net/cluqumngs/BDD_COMMUN/Illumina/FASTQ/Legionella-Amplicons-{sequencing_ID}" >&2
	echo "   -w,--work      [path]          folder where all the output files will be written, by default : 
                                                /srv/scratch/iai/bachcl/result/Legionella/23S-5S/{sequencing_ID}/{analyse_ID}" >&2
 	echo "   -m,--tmp       [path]          temporary folder where the input files will be stored, by default : 
                                                /srv/scratch/iai/bachcl/Raw_fastq/Legionella/23S-5S/{sequencing_ID}" >&2
 	echo "   -s,--save      [path]          folder where the input files will be saved, by default : 
                                                /srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/Raw_fastq/23S-5S/{sequencing_ID}" >&2
 	echo "   -o,--output    [path]          folder where the final output files will be written, by default : 
                                                /srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_results/23S-5S/{sequencing_ID}/{analyse_ID}" >&2
 	echo "   -pe,--paired   [True/False]    PE (True) or SE (False) Illumina sequencing, by default : True" >&2
 	echo "   -a,--all       [True/False]    analyse all the data in a single file (True) or separately (False), by default : True" >&2
 	echo "   -t,--trim      [True/False]    remove adaptaters and bad quality reads (QPhred < 30), by default : False" >&2
	echo >&2
 	echo "   -h, --help                     write this report and exit" >&2
    echo >&2
}

usage() {
    echo "Usage: $0 -d <seq_id> [options]"
    echo ""
    echo "Options:"
    echo "  -d, --seq_id   Sequencing ID (required)"
    echo "  -i, --input    Input folder"
    echo "  -o, --output   Output folder"
    echo "  -s, --save     Save folder"
    echo "  -m, --tmp      Temporary folder"
    echo "  -w, --work     Work folder"
    echo "  -pe, --paired  Paired-end (True/False)"
    echo "  -a, --all      One or separately (True/False)"
    echo "  -t, --trim     Trimming option (True/False)"
    echo "  -h, --help     Help"
}

################################################################################
# Configuration
## Variables init
sequencing_id=""
input_folder=""
output_folder_prefix=""
save_folder_prefix=""
tmp_folder_prefix=""
work_folder_prefix=""
paired_end=""
all_in_one=""
trim=""
analyse_id=""

## Default values
output_folder_prefix="/srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_results/23S-5S"
save_folder_prefix="/srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/Raw_fastq/23S-5S"
tmp_folder_prefix="/srv/scratch/iai/bachcl/Raw_fastq/Legionella/23S-5S"
work_folder_prefix="/srv/scratch/iai/bachcl/result/Legionella/23S-5S"
paired_end="true"
all_in_one="true"
trim="false"
analyse_id=$(date +%Y%m%d)

## User values
### Check args presence
if [ $# -eq 0 ]; then
    echo "ERROR: --seq_id is mandatory"
    usage
	echo ""
    exit 1
fi

### Parsing loop
while [ $# -gt 0 ]; do
    case "$1" in

        -d|--seq_id)
            sequencing_id="${2:?ERROR: missing value for --seq_id}"
            shift 2
            ;;

        -i|--input)
            input_folder="${2:?ERROR: missing value for --input}"
            shift 2
            ;;

        -o|--output)
            output_folder_prefix="${2:?ERROR: missing value for --output}"
            shift 2
            ;;

        -s|--save)
            save_folder_prefix="${2:?ERROR: missing value for --save}"
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

        -pe|--paired)
            paired_end="${2:?ERROR: missing value for --paired}"
            shift 2
            ;;

        -a|--all)
            all_in_one="${2:?ERROR: missing value for --all}"
            shift 2
            ;;

        -t|--trim)
            trim="${2:?ERROR: missing value for --trim}"
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

### Required argument check
if [[ -z "${sequencing_id}" ]]; then
    echo "ERROR: --seq_id is mandatory"
    usage
	echo ""
    exit 1
fi

## Generate folder names
if [[ -z "${input_folder}" ]]; then
    input_folder="/srv/net/cluqumngs/BDD_COMMUN/Illumina/FASTQ/Legionella-Amplicons-${sequencing_id}"
fi
output_folder="${output_folder_prefix}/${sequencing_id}/${analyse_id}"
save_folder="${save_folder_prefix}/${sequencing_id}"
tmp_folder="${tmp_folder_prefix}/${sequencing_id}"
work_folder="${work_folder_prefix}/${sequencing_id}/${analyse_id}/work"
result_folder="${work_folder_prefix}/${sequencing_id}/${analyse_id}"


################################################################################
# start script
## Variables for launching nextflow
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="${script_dir}/../config/qiime2.config"
pipeline_file="${script_dir}/../workflows/workflow_qiime2.nf"
nf_exec="${script_dir}/../nextflow_25.10.4"

echo "START -----------------------------------------------------------------------------------------------------------------"
echo ""

## Copy raw data from input server to calculation engine
echo "--- SAVING INPUT DATA -----------------------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

mkdir -p "${tmp_folder}"
chmod -R 777 "${tmp_folder}"
rsync -avQ --ignore-existing "${input_folder}/" "${tmp_folder}/"
echo ""

echo "--- FINISHED - to TMP FOLDER ----------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

## Copy raw data from input server to storage server
mkdir -p "${save_folder}"
rsync -avQ --ignore-existing "${input_folder}/" "${save_folder}/"
echo ""

echo "--- FINISHED - to SAVE FOLDER ---------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

## Start Qiime2 analysis
echo "--- QIIME2 ANALYSIS STARTING ----------------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "L'analyse du run Legionella-Amplicons-${sequencing_id} est en cours" \
# | mail -s "Analyse de Legionella-Amplicon-${sequencing_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

# TODO : modifier le chemin pour -f ?
k5start -U -f /home/chu-lyon.fr/ginevrach/login.kt \
    -- "${nf_exec}"  \
    -C "${config_file}"  \
    run "${pipeline_file}"  \
    --suffix "${sequencing_id}" \
    --input_dir "${tmp_folder}" \
    -w "${work_folder}" \
    --result "${result_folder}" \
    --paired_end "${paired_end}" \
    --all_in_one "${all_in_one}" \
    --trimming "${trim}" \
    -with-trace "${result_folder}/trace_${sequencing_id}_${analyse_id}.txt" \
    -with-report "${result_folder}/report_${sequencing_id}_${analyse_id}.html" \
    || LOG="error"

echo "--- FINISHED --------------------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

# Copy results from calculation engine to storage server
echo "--- SAVING OUTPUT DATA and REMOVING TMP DATA ------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

mkdir -p "${output_folder}"
rsync -avQ --exclude='*/' "$result_folder/" "$output_folder/"
echo ""
### NB : synchronising only the first level files, rest not needed

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
# echo "L'analyse du run Legionella-Amplicon-${sequencing_id} est disponible ici Z:\04_CNR_Legionella\NGS_results\23S-5S\\${sequencing_id}\\${analyse_id}" \
# | mail -s "Analyse de Legionella-Amplicon-${sequencing_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

echo "END -------------------------------------------------------------------------------------------------------------------"
