#!/bin/bash

################################################################################
#                                                                              #
# start_blast_amplicons.sh version 1                                           #
#                                                                              #
# Aurelie PETICCA, last update: 2026-04                                        #
# Christophe GINEVRA                                                           #
#                                                                              #
# Aim: Launch for Blastn Amplicons nextflow pipeline                           #
#                                                                              #
# Usage:  start_blast_amplicons.sh sequencing_ID [options]                     #
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
 	echo "   -c,--config    [path]          nextflow config file to use, by default : 
                                                <nextflow_folder>/config/blast_amplicons.config" >&2
 	echo "   -i,--input     [path]          folder containing the sequencing data, by default : 
                                                /srv/net/cluqumngs/BDD_COMMUN/Illumina/FASTQ/Legionella-Amplicons-{sequencing_ID}" >&2
	echo "   -w,--work      [path]          folder where all the output files will be written, by default : 
                                                /srv/scratch/iai/bachcl/result/Legionella/23S-5S/{sequencing_ID}/{analyse_id}_Blast-amplicons" >&2
 	echo "   -m,--tmp       [path]          temporary folder where the input files will be stored, by default : 
                                                /srv/scratch/iai/bachcl/Raw_fastq/Legionella/23S-5S/{sequencing_ID}" >&2
 	echo "   -s,--save      [path]          folder where the input files will be saved, by default : 
                                                /srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/Raw_fastq/23S-5S/{sequencing_ID}" >&2
 	echo "   -o,--output    [path]          folder where the final output files will be written, by default : 
                                                /srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_results/23S-5S/{sequencing_ID}/{analyse_id}_Blast-amplicons" >&2
 	echo "   -a,--adapter   [True/False]    remove Illumina adaptaters during trimming step, by default : True" >&2
 	echo "   -e,--deconta   [True/False]    decontamination of reads against a database, by default : False" >&2
 	echo "   -n,--down      [float]         percentage of reads retained for analysis, by default : 1 (=100%)" >&2
 	echo "   -k,--kraken    [True/False]    classifies sequence fragments with Kraken2, by default : True" >&2
	echo >&2
 	echo "   -h, --help                     write this report and exit" >&2
    echo >&2
}

usage() {
    echo "Usage: $0 -d <seq_id> [options]"
    echo ""
    echo "Options:"
    echo "  -d, --seq_id   Sequencing ID (required)"
    echo "  -c, --config   Config file"
    echo "  -i, --input    Input folder"
    echo "  -w, --work     Work folder"
    echo "  -m, --tmp      Temporary folder"
    echo "  -s, --save     Save folder"
    echo "  -o, --output   Output folder"
 	echo "  -a, --adapter  Remove adapters"
 	echo "  -e, --deconta  Decontamination"
 	echo "  -n, --down     Downsampling"
 	echo "  -k, --kraken   Kraken2 classification"
    echo "  -h, --help     Help"
}

################################################################################
# Configuration
## Variables init
sequencing_id=""
config_file=""
input_folder=""
output_folder_prefix=""
save_folder_prefix=""
tmp_folder_prefix=""
work_folder_prefix=""
paired_end=""
adapter=""
deconta=""
downsampling=""
down_to=""
kraken=""
analyse_id=""

## Default values
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="${script_dir}/../config/blast_amplicons.config"

output_folder_prefix="/srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_results/23S-5S"
save_folder_prefix="/srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/Raw_fastq/23S-5S"
tmp_folder_prefix="/srv/scratch/iai/bachcl/Raw_fastq/Legionella/23S-5S"
work_folder_prefix="/srv/scratch/iai/bachcl/result/Legionella/23S-5S"
paired_end="true"
adapter="true"
deconta="false"
downsampling="false"
down_to=1
kraken="true"
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

        -c|--config)
            config_file="${2:?ERROR: missing value for --config}"
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

        -a|--adapter)
            adapter="${2:?ERROR: missing value for --adapter}"
            shift 2
            ;;

        -e|--deconta)
            deconta="${2:?ERROR: missing value for --deconta}"
            shift 2
            ;;

        -n|--down)
            down_to="${2:?ERROR: missing value for --down}"

            #if down_to < 1
            if awk "BEGIN {exit !($down_to < 1)}"; then
                downsampling="true"
            else
                downsampling="false"
            fi

            shift 2
            ;;

        -k|--kraken)
            kraken="${2:?ERROR: missing value for --kraken}"
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
output_folder="${output_folder_prefix}/${sequencing_id}/${analyse_id}_Blast-amplicons"
save_folder="${save_folder_prefix}/${sequencing_id}"
tmp_folder="${tmp_folder_prefix}/${sequencing_id}"
work_folder="${work_folder_prefix}/${sequencing_id}/${analyse_id}_Blast-amplicons/work"
result_folder="${work_folder_prefix}/${sequencing_id}/${analyse_id}_Blast-amplicons"


################################################################################
# start script
## Variables for launching nextflow
pipeline_file="${script_dir}/../workflow_blast_amplicons.nf"
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

## Start Blast analysis
echo "--- BLAST AMPLICONS ANALYSIS STARTING ------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "L'analyse du run Legionella-Amplicons-${sequencing_id} est en cours" \
# | mail -s "Analyse de Legionella-Amplicon-${sequencing_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

k5start -U -f /home/chu-lyon.fr/ginevrach/login.kt \
    -- "${nf_exec}"  \
    -C "${config_file}"  \
    run "${pipeline_file}"  \
    --suffix "${sequencing_id}" \
    --input_dir "${tmp_folder}" \
    -w "${work_folder}" \
    --result "${result_folder}" \
    --paired_end "${paired_end}" \
    --adapter "${adapter}" \
    --decontamination "${deconta}" \
    --downsampling "${downsampling}" \
    --bbtools_downsampled "${down_to}" \
    --kraken2_assign "${kraken}" \
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
