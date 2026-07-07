#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 80 ]; then
    echo "ERROR: 80 arguments expected, got $#"
    exit 1
fi

# -------------------------
# Core paths
# -------------------------
suffix="${1}"
input_dir="${2}"
result_dir="${3}"

# -------------------------
# Flags
# -------------------------
paired_end="${4}"
adapters="${5}"
decontamination="${6}"
downsampling="${7}"
momps="${8}"
snpeff_other="${9}"

bbtools_downsampled="${10}"

# -------------------------
# FASTP
# -------------------------
min_quality="${11}"
min_length="${12}"

# -------------------------
# BBWRAP
# -------------------------
bbwrap_ref="${13}"
bbwrap_path="${14}"
bbwrap_min_id="${15}"
bbwrap_max_indel="${16}"
bbwrap_bwr="${17}"
bbwrap_bw="${18}"
bbwrap_min_hits="${19}"
bbwrap_qtrim="${20}"
bbwrap_trimq="${21}"
bbwrap_qin="${22}"

# -------------------------
# KRAKEN2
# -------------------------
kraken2_db="${23}"
format_mpa="${24}"

# -------------------------
# ELGATO
# -------------------------
spp_target="${25}"
min_ratio_target="${26}"
min_ratio_legio="${27}"
min_ratio_legia="${28}"
elgato_depth="${29}"

# -------------------------
# MINIMAP2
# -------------------------
minimap_ref="${30}"
minimap_frag="${31}"
minimap_optF="${32}"
minimap_optk="${33}"
minimap_optw="${34}"
minimap_optA="${35}"
minimap_optB="${36}"
minimap_optO="${37}"
minimap_optE="${38}"
minimap_optr="${39}"
minimap_optp="${40}"
minimap_optN="${41}"
minimap_optf="${42}"
minimap_optn="${43}"
minimap_optm="${44}"
minimap_opts="${45}"
minimap_optg="${46}"
minimap_optheap="${47}"
minimap_optsec="${48}"

# -------------------------
# FREEBAYES
# -------------------------
freeb_targets="${49}"
freeb_theta="${50}"
freeb_ploidy="${51}"
freeb_best_n="${52}"
freeb_haplo_len="${53}"
freeb_max_it="${54}"
freeb_max_dep="${55}"
freeb_min_mapqual="${56}"
freeb_min_basequal="${57}"
freeb_min_var="${58}"
freeb_min_dep="${59}"

# -------------------------
# BCFTOOLS
# -------------------------
bcf_min_freq="${60}"
bcf_qa="${61}"

# -------------------------
# SNPEFF
# -------------------------
snpeff_amr_config="${62}"
snpeff_amr_scheme="${63}"
snpeff_other_config="${64}"
snpeff_other_scheme="${65}"

# -------------------------
# ASSEMBLY
# -------------------------
min_length_contig="${66}"

# -------------------------
# FASTANI
# -------------------------
fastani_genomes="${67}"
fastani_min="${68}"

# -------------------------
# REPORTREE
# -------------------------
rep_metadata="${69}"
rep_partition="${70}"
rep_interest="${71}"
rep_zoom="${72}"
rep_site_inclusion="${72}"
rep_min_allele="${74}"
rep_max_allele="${75}"
rep_loci_called="${76}"
rep_col_metadata="${77}"

# -------------------------
# CHEWBBACA / SPECIFIC ALLELES
# -------------------------
lb_set_json="${78}"
lp_set_json="${79}"
alleles_set_json="${80}"

software_track_file="pipeline_${suffix}.txt"


# -------------------------
# File content
# -------------------------
{
echo "LEGIONELLA ASSEMBLY + MLST CONFIGURATION"
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
    echo "Sequencing type     : Paired-end (PE)"
else
    echo "Sequencing type     : Single-end (SE)"
fi

if [ "${adapters}" = true ]; then
    echo "Adapters trimming   : Enabled"
else
    echo "Adapters trimming   : Disabled"
fi

if [ "${decontamination}" = true ]; then
    echo "Decontamination     : Enabled"
else
    echo "Decontamination     : Disabled"
fi

if [ "${downsampling}" = true ]; then
    echo "Downsampling        : Enabled"
    echo "Reads kept fraction : ${bbtools_downsampled}"
else
    echo "Downsampling        : Disabled"
fi

if [ "${snpeff_other}" = true ]; then
    echo "SnpEff other genes  : Enabled"
else
    echo "SnpEff other genes  : Disabled"
fi

if [ "${momps}" = true ]; then
    echo "mompS typing        : Enabled"
else
    echo "mompS typing        : Disabled"
fi

echo ""

echo "FASTP FILTERING - trimming"
echo "Minimum quality     : ${min_quality}"
echo "Minimum length      : ${min_length}"
echo ""

echo "BBWRAP DECONTAMINATION"
if [ -n "${bbwrap_ref}" ]; then
    echo "Reference fasta     : ${bbwrap_ref}"
else
    echo "Reference path      : ${bbwrap_path}"
fi

echo "Minimum identity    : ${bbwrap_min_id}"
echo "Maximum indels      : ${bbwrap_max_indel}"
echo "BWR                 : ${bbwrap_bwr}"
echo "BW                  : ${bbwrap_bw}"
echo "Minimum hits        : ${bbwrap_min_hits}"
echo "Quality trim mode   : ${bbwrap_qtrim}"
echo "Trim quality        : ${bbwrap_trimq}"
echo "Input quality ASCII : ${bbwrap_qin}"
echo ""

echo "KRAKEN2 CLASSIFICATION"
echo "Database            : ${kraken2_db}"
echo "MPA format          : ${format_mpa}"
echo ""

echo "EL GATO FILTERING"
echo "Target species      : ${spp_target}"
echo "Min target ratio    : ${min_ratio_target}"
echo "Min Legionella ratio: ${min_ratio_legio}"
echo "Min Legionellaceae  : ${min_ratio_legia}"
echo "Depth threshold     : ${elgato_depth}"
echo ""

echo "MINIMAP2 MAPPING"
echo "Reference           : ${minimap_ref}"
echo "Fragment mode       : ${minimap_frag}"
echo "Option -F           : ${minimap_optF}"
echo "Option -k           : ${minimap_optk}"
echo "Option -w           : ${minimap_optw}"
echo "Option -A           : ${minimap_optA}"
echo "Option -B           : ${minimap_optB}"
echo "Option -O           : ${minimap_optO}"
echo "Option -E           : ${minimap_optE}"
echo "Option -r           : ${minimap_optr}"
echo "Option -p           : ${minimap_optp}"
echo "Option -N           : ${minimap_optN}"
echo "Option -f           : ${minimap_optf}"
echo "Option -n           : ${minimap_optn}"
echo "Option -m           : ${minimap_optm}"
echo "Option -s           : ${minimap_opts}"
echo "Option -g           : ${minimap_optg}"
echo "Heap sort           : ${minimap_optheap}"
echo "Secondary alignments: ${minimap_optsec}"
echo ""

echo "FREEBAYES VARIANT CALLING"
echo "Targets BED         : ${freeb_targets}"
echo "Theta               : ${freeb_theta}"
echo "Ploidy              : ${freeb_ploidy}"
echo "Best alleles        : ${freeb_best_n}"
echo "Haplotype length    : ${freeb_haplo_len}"
echo "Max iterations      : ${freeb_max_it}"
echo "Max depth fraction  : ${freeb_max_dep}"
echo "Min mapping quality : ${freeb_min_mapqual}"
echo "Min base quality    : ${freeb_min_basequal}"
echo "Min variant count   : ${freeb_min_var}"
echo "Min coverage depth  : ${freeb_min_dep}"
echo ""

echo "BCFTOOLS FILTERING"
echo "Minimum frequency   : ${bcf_min_freq}"
echo "Minimum QA          : ${bcf_qa}"
echo ""

echo "SNPEFF"
echo "Config AMR          : ${snpeff_amr_config}"
echo "Scheme AMR          : ${snpeff_amr_scheme}"
echo "Config Other        : ${snpeff_other_config}"
echo "Scheme Other        : ${snpeff_other_scheme}"
echo ""

echo "ASSEMBLY"
echo "Min contig length   : ${min_length_contig}"
echo ""

echo "FASTANI"
echo "Genome list         : ${fastani_genomes}"
echo "Minimal value       : ${fastani_min}"
echo ""

echo "REPORTREE"
echo "Metadata file       : ${rep_metadata}"
echo "Partition folder    : ${rep_partition}"
echo "Interest threshold  : ${rep_interest}"
echo "Zoom view           : ${rep_zoom}"
echo "Threshold zoom loci : ${rep_site_inclusion}"
echo "Min allele distance : ${rep_min_allele}"
echo "Max allele distance : ${rep_max_allele}"
echo "Loci called cutoff  : ${rep_loci_called}"
echo "Metadata column     : ${rep_col_metadata}"
echo ""

echo "CHEWBBACA LB SET"
echo "${lb_set_json}"
echo ""

echo "CHEWBBACA LP SET"
echo "${lp_set_json}"
echo ""

echo "SPECIFIC ALLELES EXTRACTED"
echo "${alleles_set_json}"
echo ""

echo "CONFIGURATION COMPLETE"
echo ""

echo "--------------------------------------------------------------------------------"
echo "SOFTWARES VERSION"
echo ""

} > "${software_track_file}"