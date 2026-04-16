#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 20 ]; then
    echo "ERROR: 20 arguments expected, got $#"
    exit 1
fi

# -------------------------
# Arguments
# -------------------------
input_dir="$1"
result_dir="$2"
suffix="$3"

paired_end="$4"
all_in_one="$5"
trimming="$6"

min_quality="$7"
min_length="$8"

trim_left_f="$9"
trim_left_r="${10}"
trunc_len_f="${11}"
trunc_len_r="${12}"
n_threads="${13}"
reads_learn="${14}"
fold_parents="${15}"

db="${16}"
reads="${17}"
taxa="${18}"

confidence="${19}"
n_jobs="${20}"

software_track_file="pipeline_${suffix}.txt"

# -------------------------
# File content
# -------------------------
{
echo "QIIME2 - AMPLICONS ANALYSIS CONFIGURATION"
echo ""

echo "Generated: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

echo "GENERAL SETTINGS"
echo "Input folder  : ${input_dir}"
echo "Output folder : ${result_dir}"
echo "Suffix        : ${suffix}"
echo ""

echo "ANALYSIS STRATEGY"

if [ "${paired_end}" = true ]; then
    echo "Sequencing type : Paired-end (PE)"
else
    echo "Sequencing type : Single-end (SE)"
fi

if [ "${all_in_one}" = true ]; then
    echo "Sample handling : All samples processed together"
else
    echo "Sample handling : Samples processed separately"
fi

if [ "${trimming}" = true ]; then
    echo "Trimming        : Enabled"
else
    echo "Trimming        : Disabled"
fi

echo ""

echo "CUTADAPT DATA TRIMMING"
echo "Phred Score Qual. : ${min_quality}"
echo "Length min        : ${min_length}"
echo ""

echo "DADA2 DENOISING PARAMETERS"
echo "Trim left forward : ${trim_left_f} (not used if 0)"
echo "Trim left reverse : ${trim_left_r} (not used if 0)"
echo "Trunc length F    : ${trunc_len_f}"
echo "Trunc length R    : ${trunc_len_r}"
echo "Threads           : ${n_threads}"
echo "Reads for model   : ${reads_learn}"
echo "Fold parents      : ${fold_parents}"
echo ""

echo "CLASSIFIER TRAINING"
echo "Database          : ${db}"
echo "Reference reads   : ${reads}"
echo "Taxonomy file     : ${taxa}"
echo ""

echo "TAXONOMIC CLASSIFICATION"
echo "Confidence threshold : ${confidence}"
echo "Number of jobs       : ${n_jobs}"
echo ""

echo "CONFIGURATION COMPLETE"
echo ""
echo "--------------------------------------------------------------------------------"
echo "SOFTWARES VERSION"
echo ""

} > "$software_track_file"