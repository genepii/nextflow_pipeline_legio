#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Processes for workflow_strain_comparison.nf


// -----------------------------------------------------------------------------
/*
* Find FASTA assemblies
* Input   : metadata table
* Output  : metadata with FASTA file path
* Purpose : recursively locate the most recent FASTA assembly for each sample
*/
process FIND_FASTA {
    label 'python'
    publishDir "${params.result}/ASSETS", mode: 'copy',
        pattern: "RealPath_metadata_fasta.tsv"

    input:
        path(metadata)

    output:
        path("metadata_fasta.tsv"), emit: newpath
        path("RealPath_metadata_fasta.tsv"), emit: oldpath

    script:
    def fasta_prefix = params.fasta_prefix ?: false
    def fasta_suffix = params.fasta_suffix ?: false
    def fasta_ext = params.fasta_ext ?: false

    def prefix_arg = (
        fasta_prefix instanceof List
        ? fasta_prefix.join(" ")
        : fasta_prefix.toString()
    )

    def suffix_arg = (
        fasta_suffix instanceof List
        ? fasta_suffix.join(" ")
        : fasta_suffix.toString()
    )

    def ext_arg = (
        fasta_ext instanceof List
        ? fasta_ext.join(" ")
        : fasta_ext.toString()
    )

    """
    strain_comparison_find_fasta.py \
        --metadata ${metadata} \
        --output metadata_fasta.tsv \
        --input ${params.input} \
        --directories ${params.path_fasta.join(' ')} \
        --prefix ${prefix_arg} \
        --suffix ${suffix_arg} \
        --extensions ${ext_arg}
    """
}

/*
* Find paired FASTQ files
* Input   : metadata table
* Output  : metadata table with READS_1 and READS_2 paths
* Purpose : recursively locate the most recent FASTQ R1/R2 files for each sample
*/
process FIND_FASTQ {
    label 'python'
    publishDir "${params.result}/ASSETS", mode: 'copy',
        pattern: "RealPath_metadata_fastq.tsv"

    input:
        path(metadata)

    output:
        path("metadata_fastq.tsv"), emit: newpath
        path("RealPath_metadata_fastq.tsv"), emit: oldpath

    script:
    def fastq_prefix = params.fastq_prefix ?: false
    def fastq_suffix = params.fastq_suffix ?: false
    def fastq_ext = params.fastq_ext ?: false

    def prefix_arg = (
        fastq_prefix instanceof List
        ? fastq_prefix.join(" ")
        : fastq_prefix.toString()
    )

    def suffix_arg = (
        fastq_suffix instanceof List
        ? fastq_suffix.join(" ")
        : fastq_suffix.toString()
    )

    def ext_arg = (
        fastq_ext instanceof List
        ? fastq_ext.join(" ")
        : fastq_ext.toString()
    )

    """
    strain_comparison_find_fastq.py \
        --metadata ${metadata} \
        --output metadata_fastq.tsv \
        --input ${params.input} \
        --directories ${params.path_fastq.join(' ')} \
        --prefix ${prefix_arg} \
        --suffix ${suffix_arg} \
        --extensions ${ext_arg}
    """
}

/*
* Build complete metadata table
* Input   : metadata table with FASTA paths + metadata table with FASTQ paths
* Output  : complete metadata table
* Purpose : merge FASTA and FASTQ search results into a single metadata table
*/
process BUILD_METADATA {
    label 'python'
    publishDir "${params.result}/ASSETS", mode: 'copy'

    input:
        path(fasta_metadata)
        path(fastq_metadata)

    output:
        path("Metadata_${params.suffix}.tsv")

    script:
    """
    strain_comparison_build_metadata.py \
        --fasta ${fasta_metadata} \
        --fastq ${fastq_metadata} \
        --output Metadata_${params.suffix}.tsv
    """
}

/*
* Check input files availability
* Input   : metadata table with FASTA and FASTQ paths
* Output  : missing files table + FASTA-ready metadata + FASTQ-ready metadata
* Purpose : identify samples missing sequencing files and split metadata
*/
process CHECK_METADATA {
    label 'python'

    publishDir "${params.result}/LOGS", mode: 'copy',
        pattern: "*.log"

    input:
        path(metadata)

    output:
        path("${params.suffix}_data_MISSING.log")
        path("Metadata_fasta.tsv"), emit: fasta
        path("Metadata_fastq.tsv"), emit: fastq

    script:
    """
    strain_comparison_check_metadata.py \
        --metadata ${metadata} \
        --missing ${params.suffix}_data_MISSING.log \
        --fasta Metadata_fasta.tsv \
        --fastq Metadata_fastq.tsv
    """
}


// -----------------------------------------------------------------------------
/*
* Paired-end reads filtering/trimming step
* Input   : R1 and R2 FASTQ files per sample
* Output  : trimmed paired FASTQ files
* Purpose : adapter trimming + quality filtering + length filtering
*/
process TRIM_FASTP {
    label 'fastp'
    publishDir "${params.result}/dev/0-1_Trimmed", mode: 'copy'

    input:
        tuple val(sample_id),
            val(comparison),
            path(r1),
            path(r2)

    output:
        tuple val(sample_id),
            val(comparison),
            path("${comparison}/${sample_id}_trimR1.fastq.gz"),
            path("${comparison}/${sample_id}_trimR2.fastq.gz")
        
    script:
    def adapters = params.adapters ? 
        "--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" : 
        "--disable_adapter_trimming"

    """
    mkdir "${comparison}"
    fastp \
        -i ${r1} \
        -I ${r2} \
        -o "${comparison}/${sample_id}_trimR1.fastq.gz" \
        -O "${comparison}/${sample_id}_trimR2.fastq.gz" \
        --qualified_quality_phred ${params.min_quality} \
        --length_required ${params.min_length} \
        --thread ${task.cpus} \
        ${adapters}
    """
}


// -----------------------------------------------------------------------------
/*
* Run Snippy variant calling
* Input   : sample FASTQ pair + comparison reference
* Output  : Snippy result directory
* Purpose : call variants for each sample against its Comparison reference
*/
process SNP_SNIPPY {
    label 'snippy'

    publishDir "${params.result}/dev", mode: "copy",
        pattern: "1_Snippy-${type}/*${comparison}/${sample_id}*"

    input:
        val(type)
        tuple val(sample_id),
            val(comparison),
            path(reference),
            path(r1),
            path(r2)

    output:
        tuple val(sample_id),
            val(comparison),
            val(type),
            path(reference),
            path("1_Snippy-${type}/*${comparison}/${sample_id}")

    script:
    """
    if [[ "${type}" == "ST" ]]; then
        out="1_Snippy-${type}/${type}${comparison}"
    else
        out="1_Snippy-${type}/${comparison}"
    fi

    snippy \
        --cpu ${task.cpus} \
        --reference ${reference} \
        --pe1 ${r1} \
        --pe2 ${r2} \
        --outdir \${out}/${sample_id} \
        ${params.snippy_opts}
    """
}

/*
* Check whether a Snippy Core analysis should be performed
* Input   : comparison identifier + reference + list of Snippy result directories
* Output  : comparison information + run_core flag file + comparison log
* Purpose : skip Snippy Core if only one sample or all VariantTotal = 0
*/
process CHECK_SNIPPY_CORE {
    label 'snippy'

    input:
    tuple val(comparison),
          val(type),
          path(reference),
          val(snippy_dir_list)

    output:
    tuple val(comparison),
          val(type),
          path(reference),
          val(snippy_dir_list),
          path("run_core.txt"),
          emit: list
    path("snippy_${comparison}.log"), emit: log

    script:
    def n = snippy_dir_list.size()

    def has_variant = snippy_dir_list.any { dir ->
        def snps = file("${dir}/snps.txt")
        snps.exists() && !snps.text.contains("VariantTotal\t0")
    }

    def run_core = (n > 1 && has_variant) ? 1 : 0

    def message = ""
    if (n == 1) {
        message = "ONLY ONE SAMPLE"
    } else if (!has_variant) {
        message = "NO VARIANT DETECTED"
    }

    """
    echo ${run_core} > run_core.txt
    # Note: Not used but could be to remove tuple from channel later

    if [[ ${run_core} -eq 0 ]]; then
        echo -e "${comparison}\t${message}" > "snippy_${comparison}.log"
    else
        touch "snippy_${comparison}.log"
    fi
    """
}

/*
* Run Snippy-core comparison
* Input   : reference FASTA + Snippy result directories
* Output  : core SNP alignment
* Purpose : generate core SNP alignment between strains of the same Comparison group
*/
process SNP_SNIPPY_CORE {
    label 'snippy'

    publishDir "${params.result}", mode: "copy",
        pattern: "1_Snippy-${type}/*${comparison}/*"

    input:
        tuple val(comparison),
            val(type),
            path(reference),
            path(snippy_dir)

    output:
        tuple val(comparison),
            path("1_Snippy-${type}/*/*"),
            emit: snippy_core
        tuple val(comparison),
            path("1_Snippy-${type}/*/core.clean.aln"), 
            emit: alignment
        tuple val(comparison),
            val(type),
            path("1_Snippy-${type}/*/core.aln"), 
            emit: core
        path("snippy-core_${comparison}.log"), emit: log

    script:
    """
    if [[ "${type}" == "ST" ]]; then
        out="1_Snippy-${type}/${type}${comparison}"
    else
        out="1_Snippy-${type}/${comparison}"
    fi

    mkdir -p "\${out}"

    touch "\${out}/core.aln"
    touch "\${out}/core.clean.aln"
    touch "snippy-core_${comparison}.log"

    snippy-core \
        --ref ${reference} \
        ${snippy_dir.join(' ')} > \${out}/core.log 2>&1 || true

    mv core.* "\${out}/."

    if awk -F '\\t' 'NR>1 && \$5 != 0 {found=1} END {exit !found}' "\${out}/core.txt"; then
        snippy-clean_full_aln "\${out}/core.full.aln" > "\${out}/core.clean.aln"
    elif [[ ! -s "\${out}/core.aln" ]]; then
        echo -e "${comparison}\tNOT ENOUGH VARIANT" > "snippy-core_${comparison}.log"
    fi
    """
}
// NB : snippy-core produces useful outputs on failure, so capture errors and use || true to avoid failure.


/*
* Compute pairwise SNP distances from a core SNP alignment
* Input   : comparison group identifier + output directory + core SNP alignment FASTA
* Output  : pairwise SNP distance matrix
* Purpose : calculate pairwise distances between strains within the same Comparison group
*/
process SNP_DIST_CORE {
    label 'snippy'
    publishDir "${params.result}", mode: "copy"

    input:
        tuple val(comparison),
            val(type),
            path(core)

    output:
        path("1_Snippy-${type}/*/distances_w_recomb.tsv")

    script:
    """
    if [[ "${type}" == "ST" ]]; then
        out="1_Snippy-${type}/${type}${comparison}"
    else
        out="1_Snippy-${type}/${comparison}"
    fi

    mkdir -p \${out}

    if [ ! -s "${core}" ]; then
        touch "\${out}/distances_w_recomb.tsv"
    else
        afa-pairwise.pl ${core} > "\${out}/distances_w_recomb.tsv"
    fi
    """
}

/*
* Write Snippy-core comparison warnings
* Input   : list of skipped comparisons and associated warning messages
* Output  : Snippy comparison log file
* Purpose : generate a log file listing comparisons not processed by Snippy-core
*/
process SNIPPY_CORE_LOG {
    label 'python'
    publishDir "${params.result}/LOGS", mode: "copy"

    input:
        val(type)
        path(log_entries)

    output:
        path("Snippy-${type}.log")

    script:
    """
    echo -e "Comparison\tINFO" > "Snippy-${type}.log"
    cat ${log_entries} >> "Snippy-${type}.log"
    """
}


// -----------------------------------------------------------------------------
/*
* Run IQ-TREE phylogenetic inference
* Input   : core SNP alignment generated by Snippy-core
* Output  : maximum likelihood tree files
* Purpose : generate a phylogenetic tree from the core SNP alignment
*/
process PHYLOTREE_IQTREE {
    label 'iqtree'
    publishDir "${params.result}", mode: "copy"

    input:
        val(type)
        tuple val(comparison),
            path(core_alignment)

    output:
        tuple val(type),
            val(comparison),
            path("2_Phylogeny/${type}${comparison}/${type}${comparison}.treefile"),
            emit: tree
        path("2_Phylogeny/${type}${comparison}")
        tuple val("IQtree-${type}"),
            path("2_Phylogeny/${type}${comparison}/${type}${comparison}.tree.log"),
            emit: log

    script:
    """
    out="2_Phylogeny/${type}${comparison}"

    mkdir -p "\${out}"

    if [ ! -s "${core_alignment}" ]; then
        touch "\${out}/${type}${comparison}.treefile"
        echo -e "${comparison}\tEMPTY ALIGNMENT" > "\${out}/${type}${comparison}.tree.log"
        exit 0
    else
        #TODO: ajout de bootstrap ultrarapide + support SH-aLRT (robustesse)
        nseq=\$(grep -c "^>" ${core_alignment})
        
        if [[ \${nseq} -le 1 ]]; then
            touch "\${out}/${type}${comparison}.treefile"
            echo -e "${comparison}\tSKIPPING IQTREE (\${nseq} seq)" > ${type}${comparison}.tree.log
            exit 0
        elif [[ \${nseq} -ge 4 ]]; then
            opt_bb="-bb ${params.iqtree_bootstrap}"
        else
            opt_bb=""
            echo -e "${comparison}\tNO BOOTSTRAP (\${nseq} seq)" > ${type}${comparison}.tree.log
        fi
    fi

    iqtree \
        -T ${task.cpus} \
        -st DNA \
        -alrt ${params.iqtree_alrt} \
        -m ${params.iqtree_model} \
        -s ${core_alignment} \
        --prefix ${type}${comparison} \
        \${opt_bb}

    touch ${type}${comparison}.tree.log

    mv ${type}${comparison}* \${out}/.
    """
}

/*
* Run Gubbins recombination detection
* Input   : core SNP alignment generated by Snippy-core
* Output  : Gubbins result directory
* Purpose : produce a recombination-filtered SNP alignment
*/
process FILTER_GUBBINS {
    label 'gubbins'
    publishDir "${params.result}/dev", mode: "copy"

    input:
        val(type)
        tuple val(comparison),
            path(core_alignment)

    output:
        tuple val(comparison),
            path("2_Phylogeny/*${comparison}/Gubbins/${type}${comparison}.filtered_polymorphic_sites.fasta"),
            emit: filt
        path("2_Phylogeny/*${comparison}/Gubbins")
        tuple val("Gubbins_${type}"),
            path("2_Phylogeny/*${comparison}/Gubbins/gubbins_${comparison}.err"),
            emit: log

    script:
    """
    if [[ "${type}" == *"ST"* ]]; then
        out="2_Phylogeny/ST${comparison}/Gubbins"
    else
        out="2_Phylogeny/${comparison}/Gubbins"
    fi

    mkdir -p \${out}

    if [[ ! -s "${core_alignment}" ]]; then
        echo -e "${comparison}\tEMPTY" > "\${out}/gubbins_${comparison}.err"
        touch "\${out}/${type}${comparison}.filtered_polymorphic_sites.fasta"
        exit 0
    fi

    if [[ \$(grep -c '^>' "${core_alignment}") -lt 4 ]]; then
        echo -e "${comparison}\tLESS THAN 4 SEQ" > "\${out}/gubbins_${comparison}.err"
        touch "\${out}/${type}${comparison}.filtered_polymorphic_sites.fasta"
        exit 0
    fi

    run_gubbins.py \
        --filter-percentage ${params.gubbins_percent} \
        --threads ${task.cpus} \
        --use-time-stamp \
        -p ${type}${comparison} \
        ${core_alignment} \
        > "\${out}/gubbins.log" 2>&1 || true

    touch "\${out}/gubbins_${comparison}.err"

    mv ${type}${comparison}* "\${out}/." 2>/dev/null || true
    """
}

/*
* Write Phylogenetic tree warnings
* Input   : list of skipped trees and associated warning messages
* Output  : tree log file
* Purpose : generate a log file listing visualisation not processed
*/
process TREE_LOG {
    label 'python'
    publishDir "${params.result}/LOGS", mode: "copy"

    input:
        tuple val(type),
            path(log_entries)

    output:
        path("${type}.log")

    script:
    """
    echo -e "Comparison\tINFO" > "${type}.log"
    cat ${log_entries} >> "${type}.log"
    """
}

/*
* Generate colored circular phylogenetic tree visualization
* Input   : Newick tree file and metadata file containing strain IDs
* Output  : SVG tree visualization with colored branches and legend
* Purpose : create phylogenetic tree annotated by sequence type, highlighting strains
*/
process PLOT_IQTREE {
    label 'python'
    publishDir "${params.result}", mode: "copy"

    input:
        tuple val(type),
            val(comparison),
            path(treefile)
        path(metadata)

    output:
        path("2_Phylogeny/${type}${comparison}/${type}${comparison}*")

    script:
    """
    out="2_Phylogeny/${type}${comparison}"

    mkdir -p \${out}

    if [[ ! -s ${treefile} ]]; then
        echo "WARNING: Empty or missing Newick tree. Visualization skipped." \
            > 2_Phylogeny/${type}${comparison}/${type}${comparison}.ete3.log 
    else
        strain_comparison_plot_tree.py \
            "${treefile}" \
            "${metadata}" \
            ${params.iqtree_color_by} \
            "${type}${comparison}"
        
        mv "${type}${comparison}"* \${out}/.
    fi
    """
}

// -----------------------------------------------------------------------------
/*
* Prepare allele table for GrapeTree
* Input   : sample metadata, cgMLST gene list and global chewBBACA summary table
* Output  : TSV containing Sample_ID and selected allele columns
* Purpose : extract the requested cgMLST scheme for the samples of one comparison
*/
process PREPARE_ALLELES_TSV {
    label 'python'

    publishDir "${params.result}/dev/3_cgMLST", mode: "copy",
        pattern: "${comparison}/cgMLST_${nb}genes.tsv"

    input:
        tuple path(metadata), 
            val(comparison), 
            path(genes_list), 
            val(nb),
            path(previous)

    output:
        tuple path("${comparison}/cgMLST_${nb}genes.tsv"), 
            val(nb), 
            val(comparison),
            path(genes_list),
            path(metadata)

    script:
    """
    mkdir ${comparison}

    strain_comparison_extract_cgMLST.sh \
        ${previous} \
        ${metadata} \
        ${genes_list} \
        ${comparison} \
        cgMLST_${nb}genes.tsv

    mv cgMLST_${nb}genes.tsv ${comparison}/.
    """
}

/*
* Run GrapeTree minimum spanning tree inference
* Input   : cgMLST allele table
* Output  : Newick minimum spanning tree
* Purpose : infer an MSTreeV2 tree from allele profiles
*/
process VISU_GRAPETREE {
    label 'grapetree'

    publishDir "${params.result}/3_cgMLST", mode: 'copy',
        pattern: "${comparison}/${nb}genes_GrapeTree/Lp_${nb}genes.nwk"

    input:
        tuple path(allele_tsv), 
            val(nb), 
            val(comparison),
            path(genes_list),
            path(metadata)

    output:
        tuple path("${comparison}/${nb}genes_GrapeTree/Lp_${nb}genes.nwk"), 
            val(nb), 
            val(comparison),
            path(genes_list),
            emit: tree
        tuple val("${nb}genes_GrapeTree"),
            path("${comparison}_${nb}genes_GrapeTree.log"),
            emit: log

    script:
    """
    mkdir -p "${comparison}/${nb}genes_GrapeTree"
    touch "${comparison}_${nb}genes_GrapeTree.log"
    touch "${comparison}/${nb}genes_GrapeTree/Lp_${nb}genes.nwk"

    n_profiles=\$(tail -n +2 ${allele_tsv} | cut -f2- | sort -u | wc -l)
    if [ "\$(wc -l < "${allele_tsv}")" -le 1 ]; then
        echo -e "${comparison}\tNO ALLELIC PROFILE\tNA" > "${comparison}_${nb}genes_GrapeTree.log"
        exit 0
    elif [ "\$n_profiles" -lt 2 ]; then
        echo -e "${comparison}\tONLY ONE ALLELIC PROFILE\tNA" > ${comparison}_${nb}genes_GrapeTree.log
        exit 0
    fi

    grapetree \
        -p ${allele_tsv} \
        -m ${params.grape_model} \
        > "${comparison}/${nb}genes_GrapeTree/Lp_${nb}genes.nwk"
    """
}


/*
* Prepare allele table for ReporTree
* Input   : sample metadata, cgMLST gene list and global chewBBACA summary table
* Output  : TSV containing Sample_ID and selected allele columns
* Purpose : extract the requested cgMLST scheme of one comparison
*/
process CHECK_ALLELES_TSV {
    label 'python'

    publishDir "${params.result}/dev/3_cgMLST", mode: "copy",
        pattern: "${comparison}/cgMLST_${nb}genes.tsv"

    input:
        tuple path(metadata), 
            val(comparison), 
            path(genes_list), 
            val(nb),
            path(previous)

    output:
        tuple path("${comparison}/cgMLST_${nb}genes.tsv"), 
            val(nb), 
            val(comparison),
            path(genes_list),
            path(metadata),
            path("${comparison}_${nb}genes_cgMLST.error"),
            emit: alleles
        tuple val("${nb}genes_cgMLST"),
            path("${comparison}_${nb}genes_cgMLST.log"),
            emit: log
        

    script:
    """
    mkdir ${comparison}

    strain_comparison_check_cgMLST.sh \
        ${previous} \
        ${metadata} \
        ${genes_list} \
        ${comparison} \
        "cgMLST_${nb}genes.tsv" \
        "${comparison}_${nb}genes_cgMLST.log" \
        "${comparison}_${nb}genes_cgMLST.error"

    mv "cgMLST_${nb}genes.tsv" ${comparison}/.
    """
}

/*
* Generate ReporTree visualisation
* Input   : nb genes + cgMLST allele tsv + metadata tsv
* Output  : ReporTree outputs
* Purpose : MLST newick tree and strain clustering
*/
process VISU_REPORTREE {
    label 'reportree'

    publishDir "${params.result}/3_cgMLST", mode: 'copy',
        pattern: "${comparison}/${nb}genes_ReporTree/Lp${nb}genes*"

    input:
        tuple path(allele_tsv), 
            val(nb), 
            val(comparison),
            path(genes_list),
            path(metadata),
            path(nomenclature)

    output:
        path("${comparison}/${nb}genes_ReporTree/*")
        tuple val("${nb}genes_ReporTree"),
            path("${comparison}_${nb}genes_ReporTree.log"),
            emit: log

    script:
    """
    mkdir -p ${comparison}/${nb}genes_ReporTree
    touch "${comparison}_${nb}genes_ReporTree.log"

    n_profiles=\$(tail -n +2 ${allele_tsv} | cut -f2- | sort -u | wc -l)
    if [ "\$(wc -l < "${allele_tsv}")" -le 1 ]; then
        echo -e "${comparison}\tNO ALLELIC PROFILE\tNA" > "${comparison}_${nb}genes_ReporTree.log"
        exit 0
    elif [ "\$n_profiles" -lt 2 ]; then
        echo -e "${comparison}\tONLY ONE ALLELIC PROFILE\tNA" > "${comparison}_${nb}genes_ReporTree.log"
        exit 0
    fi

    zoom_opts=""
    case "${params.rep_zoom}" in
        analyse)
            samples_of_comparison=\$(awk -F'\t' -v comp="${comparison}" '
                NR==1 {
                    for (i=1; i<=NF; i++) {
                        if (\$i=="Comparison")
                            comp_col=i
                    }
                    next
                }
                \$comp_col==comp {
                    print \$1
                }
            ' ${metadata} | paste -sd "," -)
            zoom_opts="--sample_of_interest \${samples_of_comparison} --zoom-cluster-of-interest ${params.rep_interest} --site-inclusion ${params.rep_site_inclusion}"
            ;;
        none)
            zoom_opts=""
            ;;
        *)
            zoom_opts="--sample_of_interest ${params.rep_zoom} --zoom-cluster-of-interest ${params.rep_interest} --site-inclusion ${params.rep_site_inclusion}"
            ;;
    esac

    nomenclature_opts=""
    if [ -f "${nomenclature}" ]; then
        nomenclature_opts="--nomenclature-file ${nomenclature}"
    fi

    threshold_opts=""
    if [ ${params.rep_min_allele} != "none" ]; then
        if [ ${params.rep_max_allele} != "none" ]; then
            threshold_opts="--threshold ${params.rep_min_allele}-${params.rep_max_allele}"
        else
            threshold_opts="--threshold ${params.rep_min_allele}"
        fi
    fi

    reportree.py \
        --n_proc ${task.cpus} \
        -m ${metadata} \
        -a ${allele_tsv} \
        -out Lp${nb}genes \
        -l ${genes_list} \
        --loci-called ${params.rep_loci_called} \
        --analysis ${params.rep_analysis} \
        --method ${params.rep_model} \
        --columns_summary_report "ST,Year,Origin" \
        --partitions2report 'all' \
        --metadata2report ${params.rep_col_metadata} \
        \$nomenclature_opts \
        \$zoom_opts \
        \$threshold_opts

    mv Lp${nb}genes* ${comparison}/${nb}genes_ReporTree/.
    """
}

/*
* Write Visualisation warnings
* Input   : list of skipped visualisation and associated warning messages
* Output  : visualisation log file
* Purpose : generate a log file listing visualisation not processed
*/
process VISU_LOG {
    label 'python'
    publishDir "${params.result}/LOGS", mode: "copy"

    input:
        tuple val(type),
            path(log_entries)

    output:
        path("${type}.log")

    script:
    """
    echo -e "Comparison\tINFO\tID" > "${type}.log"
    cat ${log_entries} >> "${type}.log"
    """
}