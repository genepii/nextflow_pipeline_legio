#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 33 ]; then
    echo "ERROR: 33 arguments expected, got $#"
    exit 1
fi

# -------------------------
# Arguments
# -------------------------
suffix="${1}"

input_dir="${2}"
result_dir="${3}"

paired_end="${4}"
adapter="${5}"
decontamination="${6}"
downsampling="${7}"
kraken2_assign="${8}"

min_quality="${9}"
min_length="${10}"

bbwrap_ref="${11}"
bbwrap_min_id="${12}"
bbwrap_max_indel="${13}"
bbwrap_bwr="${14}"
bbwrap_bw="${15}"
bbwrap_min_hits="${16}"
bbwrap_qtrim="${17}"
bbwrap_trimq="${18}"
bbwrap_qin="${19}"
bbwrap_path="${20}"

bbtools_downsampled="${21}"

kraken2_db="${22}"

min_overlap="${23}"
max_overlap="${24}"
dovetail_overlap="${25}"

blast_db="${26}"
perc_id="${27}"
loose_id="${28}"
query_cov="${29}"
loose_cov="${30}"
min_qlen="${31}"
loose_qlen="${32}"
delta="${33}"

software_track_file="pipeline_${suffix}.txt"

# -------------------------
# File content
# -------------------------
{
echo "BLAST AMPLICONS PIPELINE CONFIGURATION"
echo ""

echo "Generated: $(date '+%d/%m/%Y %H:%M:%S')"
echo ""

echo "GENERAL SETTINGS"
echo "Sequencing ID : ${suffix}"
echo "Input folder  : ${input_dir}"
echo "Output folder : ${result_dir}"
echo ""

echo "ANALYSIS OPTIONS"
echo "Paired-end        : ${paired_end}"
echo "Adapter trimming  : ${adapter}"
echo "Decontamination   : ${decontamination}"
echo "Downsampling      : ${downsampling}"
echo "Kraken2 assign    : ${kraken2_assign}"
echo ""

echo "FASTP FILTERING"
echo "Min quality : ${min_quality}"
echo "Min length  : ${min_length}"
echo ""

echo "BBWRAP MAPPING"
echo "Reference     : ${bbwrap_ref}         //if empty, DB path used"
echo "Min identity  : ${bbwrap_min_id}"
echo "Max indel     : ${bbwrap_max_indel}"
echo "BWR           : ${bbwrap_bwr}"
echo "BW            : ${bbwrap_bw}"
echo "Min hits      : ${bbwrap_min_hits}"
echo "Qtrim         : ${bbwrap_qtrim}"
echo "TrimQ         : ${bbwrap_trimq}"
echo "Qin           : ${bbwrap_qin}"
echo "DB path       : ${bbwrap_path}"
echo ""

echo "DOWNSAMPLING"
echo "Fraction : ${bbtools_downsampled}"
echo ""

echo "KRAKEN2"
echo "Database : ${kraken2_db}"
echo ""

echo "MERGED READS"
echo "Min. overlapping      : ${min_overlap}"
echo "Max. overlapping      : ${max_overlap}"
echo "Dovetail overlapping  : ${dovetail_overlap}"
echo ""

echo "BLAST"
echo "Database              : ${blast_db}"
echo "%identity strict      : ${perc_id}"
echo "%identity loose       : ${loose_id}"
echo "%coverage strict      : ${query_cov}"
echo "%coverage loose       : ${loose_cov}"
echo "Min. qlength strict   : ${min_qlen}"
echo "Min. qlength loose    : ${loose_qlen}"
echo "Max bitscore          : ${delta}"
echo ""

echo "--------------------------------------------------------------------------------"
echo "CONFIGURATION COMPLETE"
echo ""

} > "$software_track_file"