#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 25 ]; then
    echo "ERROR: 25 arguments expected, got $#"
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
adapters="$6"

min_quality="$7"
min_length="$8"

trim_left_f="$9"
trim_left_r="${10}"
trunc_len_f="${11}"
trunc_len_r="${12}"
reads_learn="${13}"
fold_parents="${14}"

db="${15}"
reads="${16}"
taxa="${17}"

sklearn_confidence="${18}"
blast_identity="${19}"
blast_maxaccepts="${20}"
blast_query_cov="${21}"
vsearch_identity="${22}"
vsearch_maxaccepts="${23}"
vsearch_query_cov="${24}"
classifier="${25}"

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
    echo "Sequencing type           : Paired-end (PE)"
else
    echo "Sequencing type           : Single-end (SE)"
fi

if [ "${all_in_one}" = true ]; then
    echo "Sample handling           : All samples processed together"
else
    echo "Sample handling           : Samples processed separately"
fi

if [ "${adapters}" = true ]; then
    echo "Adapters                  : Enabled"
else
    echo "Adapters                  : Disabled"
fi

echo "Classifier used           : ${classifier}"

echo ""

echo "FASTP FILTERING - trimming"
echo "Phred Score Qual. : ${min_quality}"
echo "Length min        : ${min_length}"
echo ""

echo "DADA2 DENOISING"
echo "Trim left forward : ${trim_left_f} (not used if 0)"
echo "Trim left reverse : ${trim_left_r} (not used if 0)"
echo "Trunc length F    : ${trunc_len_f}"
echo "Trunc length R    : ${trunc_len_r}"
echo "Reads for model   : ${reads_learn}"
echo "Fold parents      : ${fold_parents}"
echo ""

echo "CLASSIFIER TRAINING"
echo "Database          : ${db}"
echo "Reference reads   : ${reads}"
echo "Taxonomy file     : ${taxa}"
echo ""

echo "SKLEARN CLASSIFICATION - taxonomic classification (if sklearn classifier)"
echo "sklearn confidence threshold : ${sklearn_confidence}"
echo ""

echo "BLAST CLASSIFICATION - taxonomic classification (if blast classifier)"
echo "Min. Identity percent : ${blast_identity}"
echo "Max. number of hits   : ${blast_maxaccepts}"
echo "Min. query Coverage   : ${blast_query_cov}"
echo ""

echo "VSEARCH CLASSIFICATION - taxonomic classification (if vsearch classifier)"
echo "Min. Identity percent : ${vsearch_identity}"
echo "Max. number of hits   : ${vsearch_maxaccepts}"
echo "Min. query Coverage   : ${vsearch_query_cov}"
echo ""

echo "CONFIGURATION COMPLETE"
echo ""
echo "--------------------------------------------------------------------------------"
echo "SOFTWARES VERSION"
echo ""

} > "$software_track_file"