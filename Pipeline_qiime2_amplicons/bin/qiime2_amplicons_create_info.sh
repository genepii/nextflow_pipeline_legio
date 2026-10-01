#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 31 ]; then
    echo "ERROR: 31 arguments expected, got $#"
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
denoising="$7"

min_quality="$8"
min_length="$9"

trim_left_f="${10}"
trim_left_r="${11}"
trunc_len_f="${12}"
trunc_len_r="${13}"
reads_learn="${14}"
fold_parents="${15}"

trunc_qual="${16}"
min_overlap="${17}"
max_diffs="${18}"
min_mergelenght="${19}"

db="${20}"
reads="${21}"
taxa="${22}"

sklearn_confidence="${23}"
blast_identity="${24}"
blast_maxaccepts="${25}"
blast_query_cov="${26}"
vsearch_identity="${27}"
vsearch_maxaccepts="${28}"
vsearch_query_cov="${29}"
classifier="${30}"

kraken_db="${31}"

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

if [ "${denoising}" = true ]; then
    echo "Data Processing           : Denoising (DADA2)"
else
    echo "Data Processing           : Deduplicating (VSearch)"
fi

echo "Classifier used           : ${classifier}"

echo ""

echo "FASTP FILTERING - trimming"
echo "Phred Score Qual. : ${min_quality}"
echo "Length min        : ${min_length}"
echo ""

echo "KRAKEN - identification"
echo "Database          : ${kraken_db}"
echo ""

echo "DADA2 DENOISING"
echo "Trim left forward : ${trim_left_f} (not used if 0)"
echo "Trim left reverse : ${trim_left_r} (not used if 0)"
echo "Trunc length F    : ${trunc_len_f}"
echo "Trunc length R    : ${trunc_len_r}"
echo "Reads for model   : ${reads_learn}"
echo "Fold parents      : ${fold_parents}"
echo ""

echo "VSEARCH DEDUPLICATING"
echo "Trim base quality : ${trunc_qual} (not used if 0)"
echo "Min. overlapping  : ${min_overlap}"
echo "Max. base diff    : ${max_diffs}"
echo "Min. final length : ${min_mergelenght}"
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