#!/bin/bash

################################################################################
#                                                                              #
# start_assembly_mlst.sh version 1                                             #
#                                                                              #
# Aurelie PETICCA, last update: 2026-07                                        #
# Christophe GINEVRA                                                           #
#                                                                              #
# Aim: Launch for Assembly + MLST nextflow pipeline                            #
#                                                                              #
# Usage:  start_assembly_mlst.sh sequencing_ID [options]                       #
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
 	echo "   -d,--seq_id    [str]           SEQUENCING_ID, locally a date in format YYYYMMDD, Required" >&2
 	echo "   -c,--config    [path]          nextflow config file to use, by default : 
                                                <nextflow_folder>/config/assembly_mlst.config" >&2
 	echo "   -i,--input     [path]          folder containing the sequencing data, by default : 
                                                /srv/net/cluqumngs/BDD_COMMUN/Illumina/FASTQ/Legionella-{sequencing_ID}" >&2
	echo "   -w,--work      [path]          folder where all the output files will be written, by default : 
                                                /srv/scratch/iai/bachcl/result/Legionella/Genomes/{sequencing_ID}/{analyse_ID}_Assembly-MLST" >&2
 	echo "   -m,--tmp       [path]          temporary folder where the input files will be stored, by default : 
                                                /srv/scratch/iai/bachcl/Raw_fastq/Legionella/Genomes/{sequencing_ID}" >&2
 	echo "   -s,--save      [path]          folder where the input files will be saved, by default : 
                                                /srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/Raw_fastq/Genomes/{sequencing_ID}" >&2
 	echo "   -o,--output    [path]          folder where the final output files will be written, by default : 
                                                /srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_results/Genomes/{sequencing_ID}/{analyse_ID}_Assembly-MLST" >&2
 	echo "   -t,--adapters  [True/False]    remove Illumina adaptaters, by default : True" >&2
 	echo "   -e,--deconta   [True/False]    decontamination of reads against a database, by default : False" >&2
 	echo "   -n,--down      [float]         percentage of reads retained for analysis, by default : 1 (=100%)" >&2
 	echo "   -p,--momps     [True/False]    enable MLST research by mompS, by default : False" >&2
 	echo "   -f,--snpeff    [True/False]    enable SnpEff to search for changes other than AMR changes, by default : False" >&2
 	echo "   -z,--zoom      [str]           enable ReporTree zoom functionality on the samples analysed (analyse), all, none (by default), or sample(s) ID to zoom" >&2
 	echo "   -a,--meta      [path]          TSV file containing metadata for all the samples; must contain columns ID, Year, Origin and Linked_to; by default :
                                                /srv/net/cluqumngs/TMP_IAI/04_CNR_Legionella/Input_analysis_nextflow/Metadata_{sequencing_ID}.tsv" >&2
 	echo "   -r,--part      [path]          folder for cgMLST-based HC cluster file named Lp_{nb}genes_{partitions|MLSTchewbbaca}.tsv, by default :
                                                /srv/scratch/iai/bachcl/db/legio/ReporTree" >&2
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
    echo "  -o, --output   Output folder"
    echo "  -s, --save     Save folder"
    echo "  -m, --tmp      Temporary folder"
    echo "  -w, --work     Work folder"
 	echo "  -t, --adapters Remove Illumina adapters (True/False)"
 	echo "  -e, --deconta  Decontamination"
 	echo "  -n, --down     Downsampling"
 	echo "  -p, --momps    MompS MLST"
 	echo "  -f, --snpeff   SnpEff other than AMR"
 	echo "  -z, --zoom     ReporTree zoom functionality"
    echo "  -a, --meta     Metadata TSV file"
    echo "  -r, --part     HC cluster file"
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
adapters=""
deconta=""
downsampling=""
down_to=""
momps=""
snpeff_other=""
zoom=""
metadata_user=""
partition_folder=""
analyse_id=""

## Default values
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="${script_dir}/config/assembly_mlst.config"

output_folder_prefix="/srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_results/Genomes"
save_folder_prefix="/srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/Raw_fastq/Genomes"
tmp_folder_prefix="/srv/scratch/iai/bachcl/Raw_fastq/Legionella/Genomes"
work_folder_prefix="/srv/scratch/iai/bachcl/result/Legionella/Genomes"
adapters="true"
deconta="false"
downsampling="false"
down_to=1
momps="false"
snpeff_other="false"
zoom="none"
metadata_user_prefix="/srv/net/cluqumngs/TMP_IAI/04_CNR_Legionella/Input_analysis_nextflow/Metadata_"
partition_folder="/srv/net/cluqumngs/TMP_IAI/04_CNR_Legionella/NGS_Database/ReporTree"
analyse_id=$(date +%Y%m%d)

## User values
### Check args presence
if [ $# -eq 0 ]; then
    echo "ERROR: --seq_id is Required"
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

            # Check for FASTQ files
            if ! compgen -G "${input_folder}"/*_R1.fastq > /dev/null && \
            ! compgen -G "${input_folder}"/*_R1.fastq.gz > /dev/null; then
                echo "ERROR: No FASTQ files (*_R1.fastq or *_R1.fastq.gz) found in input folder: ${input_folder}" >&2
                exit 1
            fi
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

        -t|--adapters)
            adapters="${2:?ERROR: missing value for --adapters}"
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
            
        -p|--momps)
            momps="${2:?ERROR: missing value for --momps}"
            shift 2
            ;;
            
        -f|--snpeff)
            snpeff_other="${2:?ERROR: missing value for --snpeff}"
            shift 2
            ;;
            
        -a|--meta)
            metadata_user="${2:?ERROR: missing value for --meta}"
            shift 2
            ;;
            
        -r|--part)
            partition_folder="${2:?ERROR: missing value for --part}"
            shift 2
            ;;
            
        -z|--zoom)
            zoom="${2:?ERROR: missing value for --zoom}"

            if [[ "$zoom" =~ ^(none|all|analyse|[A-Za-z0-9-]{6,100})$ ]]; then
                :
            else
                echo "ERROR: invalid value for --zoom (allowed: analyse, none, all or sample ID list)" >&2
                exit 1
            fi
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
    echo "ERROR: --seq_id is Required"
    usage
	echo ""
    exit 1
fi

## Generate folder names
if [[ -z "${input_folder}" ]]; then
    input_folder="/srv/net/cluqumngs/BDD_COMMUN/Illumina/FASTQ/Legionella-${sequencing_id}"
fi
output_folder="${output_folder_prefix}/${sequencing_id}/${analyse_id}_Assembly-MLST"
save_folder="${save_folder_prefix}/${sequencing_id}"
tmp_folder="${tmp_folder_prefix}/${sequencing_id}"
work_folder="${work_folder_prefix}/${sequencing_id}/${analyse_id}_Assembly-MLST/work"
result_folder="${work_folder_prefix}/${sequencing_id}/${analyse_id}_Assembly-MLST"

## Input folder content check
if [[ ! -d "${input_folder}" ]] || \
   ! find "${input_folder}" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) -print -quit | grep -q .; then
    echo "ERROR: No FASTQ or FASTQ.GZ files found in ${input_folder}" >&2
    echo "     Please verify the path and consider re-downloading the data." >&2
    exit 1
fi

## Metadata file content check
if [[ -z "${metadata_user}" ]]; then
    metadata_user="${metadata_user_prefix}${sequencing_id}.txt"
fi

if [[ ! -f "${metadata_user}" ]]; then
    echo "ERROR: file not found ${metadata_user}" >&2
    exit 1
fi

### Check if TSV
header=$(head -n1 "${metadata_user}")
[[ "${header}" == *$'\t'* ]] || {
    echo "ERROR: ${metadata_user} is not a TSV file." >&2
    exit 1
}

### Check columns needed
for col in ID Year Origin; do
    grep -qw "${col}" <<< "$(tr '\t' '\n' <<< "${header}")" || {
        echo "ERROR: Missing column '${col}' in ${metadata_user}" >&2
        exit 1
    }
done

### Check FASTQs against metadata
valid_ids=$(tail -n +2 "${metadata_user}" | cut -f1 | tr '\n' '|' | sed 's/|$//')
find "${input_folder}" -maxdepth 1 -type f \( -name "*_R1.fastq" -o -name "*_R1.fastq.gz" \) | while read -r f; do
    base=$(basename "${f}")
    sample=${base%%_*}

    echo "${valid_ids}" | grep -qw "${sample}" || {
        echo "ERROR: FASTQ sample '${sample}' not found in Metadata file" >&2
        exit 1
    }
done


################################################################################
# start script
## Variables for launching nextflow
pipeline_file="${script_dir}/workflow_assembly_mlst.nf"
nf_exec="${script_dir}/../nextflow_25.10.4"

echo "START ---------------------------------------------------------------------------------------------------------------------"
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

## Start Assembly + MLST analysis
echo "--- ASSEMBLY + MLST ANALYSIS STARTING -------------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "L'analyse ASSEMBLY + MLST du run Legionella-${sequencing_id} est en cours" \
# | mail -s "Analyse Assembly + MLST Legionella-${sequencing_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

if ! k5start -U -f /home/chu-lyon.fr/ginevrach/login.kt \
    -- "${nf_exec}" \
    -C "${config_file}" \
    run "${pipeline_file}" \
    --suffix "${sequencing_id}" \
    --input_dir "${tmp_folder}" \
    -w "${work_folder}" \
    -profile local \
    --result "${result_folder}" \
    --adapters "${adapters}" \
    --decontamination "${deconta}" \
    --downsampling "${downsampling}" \
    --bbtools_downsampled "${down_to}" \
    --momps "${momps}" \
    --snpeff_other "${snpeff_other}" \
    --rep_zoom "${zoom}" \
    --rep_metadata "${metadata_user}" \
    --rep_partition "${partition_folder}" \
    -with-trace "${result_folder}/LOGS/nextflow_${sequencing_id}_${analyse_id}.txt" \
    -with-report "${result_folder}/LOGS/nextflow_${sequencing_id}_${analyse_id}.html"
then
    LOG="error"
fi

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
    "${result_folder}/" "${output_folder}/"
echo ""
### Synchronising every files/folders, dev and work not needed

## Backup and replace changed files only (ReporTree DB)
timestamp=$(date +"%Y%m%d-%H%M")
shopt -s nullglob
for src in "${result_folder}"/dev/Rsync/*.tsv; do
    [[ -s "${src}" ]] || continue
    dst="${partition_folder}/$(basename "${src}")"

    [[ -f "${dst}" ]] && cmp -s "${src}" "${dst}" && continue

    if [[ -f "${dst}" ]]; then
        mkdir -p "${partition_folder}/OLD/${timestamp}"
        mv "${dst}" "${partition_folder}/OLD/${timestamp}/"
    fi

    cp "${src}" "${dst}"
done
echo ""

echo "--- FINISHED - to SAVE FOLDER ---------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

## Remove results from calculation engine
echo "Deleting... ${work_folder}"
rm -r "${work_folder}"
rm -r "${result_folder}/dev/0-1_Trimmed" # Warning: Delete Trimmed Reads for space
rm -r "${result_folder}/dev/Rsync"       # Warning: Delete partitions.tsv / alleles.tsv after Rsync
echo "Deleting... ${tmp_folder}"
rm -r "${tmp_folder}"

echo "--- FINISHED - to DELETE --------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "L'analyse ASSEMBLY + MLST du run Legionella-${sequencing_id} est disponible ici : ${output_folder}" \
# | mail -s "Analyse Assembly + MLST Legionella-${sequencing_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

echo "END -----------------------------------------------------------------------------------------------------------------------"
