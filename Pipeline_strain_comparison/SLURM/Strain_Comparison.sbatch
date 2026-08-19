#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --gres=submitter:1
#SBATCH --output=/srv/net/cluqumngs/TMP_IAI/redmine/65198/Strain-Comparison_%j.out
#SBATCH --error=/srv/net/cluqumngs/TMP_IAI/redmine/65198/Strain-Comparison_%j.err
#SBATCH --comment="Strain Comparison pipeline, SLURM out in redmine/65198"
#SBATCH --job-name=Strain-Comparison
#SBATCH --partition=diag_iai
#SBATCH --nodelist=di6182su
#SBATCH --account=iai
#SBATCH --constraint=public,dump

set -e

### include/nextflow-conf.bash ###

# In local storage allocated by SLURM 
cd ${SCRATCH}

################################################################################
#                                                                              #
# Strain_Comparison.sbatch version 1                                           #
#                                                                              #
# Aurelie PETICCA, last update: 2026-08                                        #
# Christophe GINEVRA                                                           #
#                                                                              #
# Aim: SLURM launch for SNP + previous MLST nextflow pipeline                  #
#                                                                              #
# Usage:  sbatch Strain_Comparison.sbatch -s Comparison_Name [options]         #
# WARNING : not suitable for NGS-Web launching yet                             #
#                                                                              #
################################################################################

# Help screen
display_help() { 
    echo "Usage: sbatch Strain_Comparison.sbatch -s Comparison_Name [options]" >&2
    echo >&2
    echo >&2
    echo "   -s             [str]           name assigned to the comparison, corresponding to the suffix of the Strains_XXXX.txt metadata file. (required)" >&2
    echo "   -c             [path]          nextflow config file to use, by default : 
                                                <nextflow_folder>/config/strain_comparison.config" >&2
    echo >&2
    echo "   -h                             print this help message and exit" >&2
    echo >&2
}

usage() {
    echo "Usage: sbatch Strain_Comparison.sbatch -s <name> [options]" >&2
    echo >&2
    echo "Options:" >&2
    echo "  -s             Comparison name (required)" >&2
    echo "  -c             Config file" >&2
    echo "  -h             Show help" >&2
}

################################################################################
# Configuration
## Variables init
comparison_id=""
script_dir=""
config_file=""
pipeline_nf=""
nf_exec=""

## Nextflow default values
script_dir="${DEPOT}/pipelines/nextflow_pipeline_legio/Pipeline_strain_comparison/SLURM"

config_file="${script_dir}/config/strain_comparison.config"
pipeline_nf="${script_dir}/../workflow_strain_comparison.nf"
dump_nf="${script_dir}/dump_params.nf"
nf_exec="${script_dir}/../../nextflow_25.10.4"

### Dynamic Nextflow configuration
nextflow-conf > "${SCRATCH}/SLURM.conf"

## User values
### Check args presence
if [ $# -eq 0 ]; then
    echo "ERROR: Comparison Name is required" >&2
    usage >&2
    echo >&2
    exit 1
fi

### Parsing loop
while getopts ":s:hc:" opt; do
    case ${opt} in
    s)
        # Comparison name
        comparison_id="${OPTARG}"
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
        echo "ERROR: Option -${OPTARG} requires an argument" >&2
        exit 1
        ;;
    \?)
        echo "ERROR: Unknown option -${OPTARG}" >&2
        usage
        exit 1
        ;;
    esac
done

### Required argument check
if [[ -z "${comparison_id}" ]]; then
    echo "ERROR: Comparison Name is required" >&2
    usage
    echo >&2
    exit 1
fi

## Retrieving paths from nextflow config, contains "="
mapfile -t paths < <(
    "${nf_exec}" run "${dump_nf}" \
        -c "${config_file}" \
        --suffix "${comparison_id}" \
        -profile local \
        2>&1 | grep -E '^.*=' | cut -d= -f2-
)
output_folder="${paths[0]}"
#tmp_folder="${paths[1]}"
#work_folder="${paths[2]}"
result_folder="${paths[3]}"
analyse_id="${paths[4]}"
metadata_user="${paths[5]}"
db_snippy="${paths[6]}"

tmp_folder=${TMPDIR}/${comparison_id}/${analyse_id}_Strain-Comparison/data
work_folder=${TMPDIR}/${comparison_id}/${analyse_id}_Strain-Comparison/work

### If path empty error
for d in "${output_folder}" "${tmp_folder}" "${work_folder}" "${result_folder}"; do
    if [[ -z "$d" ]]; then
        echo "ERROR: One of the required paths is empty" >&2
        echo "     output_folder/output_dir : ${output_folder}" >&2
        echo "     tmp_folder/tmp_dir       : ${tmp_folder}" >&2
        echo "     work_folder/workDir      : ${work_folder}" >&2
        echo "     result_folder/result     : ${result_folder}" >&2
        echo >&2
        exit 1
    fi
done

### Create output paths
for d in "${output_folder}" "${tmp_folder}" "${work_folder}" "${result_folder}"; do
    mkdir -p "$d"
done

## Check Metadata file
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
for col in ID ST Year Origin Comparison; do
    grep -qw "${col}" <<< "$(tr '\t' '\n' <<< "${header}")" || {
        echo "ERROR: Missing column '${col}' in ${metadata_user}" >&2
        exit 1
    }
done


################################################################################
# start script
echo "START -----------------------------------------------------------------------------------------------------------------"
echo

## Start SNP and cgMLST Comparison
echo "--- STRAIN COMPARISON ANALYSIS STARTING -------------------------------------------------------------------------------"
echo "Start: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

### Notification
# echo "L'analyse SNP + MLST de la comparaison ${comparison_id} est en cours" \
# | mail -s "Analyse SNP + MLST comparaison ${comparison_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

"${nf_exec}" run \
    "${pipeline_nf}" \
    -c "${SCRATCH}/SLURM.conf" \
    -c "${config_file}" \
    --suffix "${comparison_id}" \
    -profile slurm \
    -with-trace "${result_folder}/LOGS/nextflow_${comparison_id}.txt" \
    -with-report "${result_folder}/LOGS/nextflow_${comparison_id}.html"

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
    "${result_folder}/" "${output_folder}/"
echo
### Synchronising every files/folders, dev and work not needed

### Add new Snippy results vs ST to Database
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

echo "--- FINISHED - to SAVE FOLDER ---------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo

## Remove results from calculation engine
echo "Deleting... ${work_folder}"
rm -r "${work_folder}"
rm -r "${result_folder}/dev/0-1_Trimmed"            # Warning: Delete Trimmed Reads for space
rm -f "${result_folder}"/dev/1_Snippy-*/*/*/*.bam*  # Warning: Delete Bam files for space
echo "Deleting... ${tmp_folder}"
rm -r "${tmp_folder}"

echo "--- FINISHED - to DELETE --------------------------------------------------------------------------------------------------"
echo "End: $(date '+%d/%m/%Y %H:%M:%S')"
echo

### Notification
# echo "L'analyse SNP + MLST du run de la comparaison ${comparison_id} est disponible ici : ${output_folder}" \
# | mail -s "Analyse SNP + MLST comparaison ${comparison_id}" christophe.ginevra@chu-lyon.fr GHE.CNR-LEGIO@chu-lyon.fr

echo "END -------------------------------------------------------------------------------------------------------------------"
