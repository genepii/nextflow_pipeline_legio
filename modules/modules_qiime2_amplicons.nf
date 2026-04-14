#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Processes for workflow_qiime2_amplicons.nf


// -----------------------------------------------------------------------------
/*
* Custom 23S–5S reference database
* Input   : FASTA sequences and taxonomy
* Output  : Naive Bayes classifier
* Purpose : adapt the classifier to the data (local issue)
*/
process import_refseq {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        val(reads_file)

    output:
        path "${params.db}.qza"

    script:
    """
    qiime tools import \
        --type FeatureData[Sequence] \
        --input-path ${reads_file} \
        --output-path "${params.db}.qza"
    """
}

process import_taxa {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        val(taxa_file)

    output:
        path "${params.db}_tax.qza"

    script:
    """
    qiime tools import \
        --type FeatureData[Taxonomy] \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path ${taxa_file} \
        --output-path "${params.db}_tax.qza"
    """
}

process generate_classifier_bayes {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        path classifier_reads
        path classifier_taxa

    output:
        path "${params.db}_classifier.qza"

    script:
    """
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads ${classifier_reads} \
        --i-reference-taxonomy ${classifier_taxa} \
        --o-classifier "${params.db}_classifier.qza"
    """
}

// -----------------------------------------------------------------------------
/*
* Remove Illumina adaptaters and bad quality reads
* Input   : tuples with sample_id, Illumina R1 and R2 (optionnal)
* Output  : tuples with sample_id, Illumina R1 and R2 (optionnal)
* Purpose : generate cleaned dataset for Qiime2 analysis
*/
process trimming_cutadapt {
    label 'cutadapt'
    publishDir "${params.result}", mode: 'copy'

    input:
        tuple val(sample_id), val(r1), val(r2)

    output:
        tuple val(sample_id),
            path("${sample_id}_trimR1.fastq.gz"),
            path(params.paired_end ? "${sample_id}_trimR2.fastq.gz" : null)

    script:
    def is_paired = params.paired_end

    def adapter = is_paired ?
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" :
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

    def quality = is_paired ?
        "-q ${params.min_quality},${params.min_quality}" :
        "-q ${params.min_quality}"

    def minlen = "--minimum-length ${params.min_length}"

    def r1_out = "${sample_id}_trimR1.fastq.gz"
    def r2_out = is_paired ? "${sample_id}_trimR2.fastq.gz" : "null"

    def inputs = is_paired ?
        "${r1} ${r2}" :
        "${r1}"

    def outputs = is_paired ?
        "-o ${r1_out} -p ${r2_out}" :
        "-o ${r1_out}"

    """
    cutadapt \
        ${adapter} \
        ${quality} \
        ${minlen} \
        ${outputs} \
        ${inputs}
    """
}

// -----------------------------------------------------------------------------
/*
* Generate TSV manifest based on params.paired_end = True or not
* Input   : tuples with sample_id, Illumina R1 and R2 (optionnal)
* Output  : TSV manifest for this sample_id
* Purpose : generate input for QIIME2 import
* Note    : FastqManifestPhred33V2 manifest
*/
process generate_manifest {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        tuple val(sample_id), val(r1), val(r2)

    output:
        tuple val(sample_id), path("${sample_id}_manifest.tsv")

    script:
    """
    if [ "${params.paired_end}" = "true" ]; then
        printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > ${sample_id}_manifest.tsv
        printf "%s\t%s\t%s\n" "${sample_id}" "${r1}" "${r2}" >> ${sample_id}_manifest.tsv
    else
        printf "sample-id\tabsolute-filepath\n" > ${sample_id}_manifest.tsv
        printf "%s\t%s\n" "${sample_id}" "${r1}" >> ${sample_id}_manifest.tsv
    fi
    """
}

/*
* Generate TSV manifest for SE or PE data (params.paired_end)
* Input   : list of TSV manifests
* Output  : TSV manifest for all sample_id
* Purpose : generate input for QIIME2 import
* Note    : FastqManifestPhred33V2 manifest
*/
process generate_manifest_all {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        tuple val(sample_id), path(manifests)

    output:
        tuple val("All-samples"), path("All-samples_manifest.tsv")

    script:
    def header = params.paired_end ? 
        "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" : 
        "sample-id\tabsolute-filepath"
    """
    echo -e "${header}" > All-samples_manifest.tsv

    for f in ${manifests}; do
        tail -n +2 "\$f" >> All-samples_manifest.tsv
    done
    """
}

// -----------------------------------------------------------------------------
/*
* Import data into QIIME2 for SE or PE data (params.paired_end)
* Input   : TSV file with path to Illumina R1 and R2 files (manifest)
* Output  : QIIME2 artifact (.qza)
* Purpose : prepare data for denoising
*/
process import_manifest {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        tuple val(sample_id), path(manifest)

    output:
        tuple val(sample_id), path("${sample_id}_demux.qza")

    script:
    def type = params.paired_end ? 
        "SampleData[PairedEndSequencesWithQuality]" : 
        "SampleData[SequencesWithQuality]"
    def format = params.paired_end ? 
        "PairedEndFastqManifestPhred33V2" : 
        "SingleEndFastqManifestPhred33V2"
    """
    qiime tools import \
        --input-path ${manifest} \
        --output-path ${sample_id}_demux.qza \
        --type "${type}" \
        --input-format "${format}"
    """
}

/*
* Generate demultiplexing summary
* Input   : imported QIIME2 artifact (.qza)
* Output  : plot .qzv with information about Qiime2 demultiplexing
* Purpose : generate quality plots (QC analysis)
*/
process qc_demux {
    label 'qiime'
    publishDir params.result, mode: 'copy'

    input:
        tuple val(sample_id), path(demux)

    output:
        tuple val(sample_id), path("${sample_id}_demux.qzv")

    script:
    """
    qiime demux summarize \
        --i-data ${demux} \
        --o-visualization ${sample_id}_demux.qzv
    """
}

// -----------------------------------------------------------------------------
/*
* Denoising by DADA2 for 23S–5S amplicons
* Input   : imported QIIME2 artifact (.qza)
* Output  : feature table and representative sequences (.qza)
* Purpose : denoise and generate ASVs
*/
process denoise_dada2 {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        tuple val(sample_id), path(demux)

    output:
        tuple val(sample_id), 
            path("${sample_id}_table-dada2_${params.suffix}.qza"), 
            emit: table_dada2
        tuple val(sample_id), 
            path("${sample_id}_stats-dada2_${params.suffix}.qza"), 
            emit: stats_dada2
        tuple val(sample_id), 
            path("${sample_id}_rep-seqs-dada2_${params.suffix}.qza"), 
            emit: rep_dada2

    script:
    def type = params.paired_end ? 
        "denoise-paired" : 
        "denoise-single"
    def trunc = params.paired_end ?
        "--p-trunc-len-f ${params.trunc_len_f} --p-trunc-len-r ${params.trunc_len_r}" :
        "--p-trunc-len ${params.trunc_len_f}"
    """
    qiime dada2 ${type} \
        --i-demultiplexed-seqs ${demux} \
        ${trunc} \
        --p-min-fold-parent-over-abundance ${params.fold_parents} \
        --p-n-reads-learn ${params.reads_learn} \
        --p-n-threads ${params.n_threads} \
        --o-representative-sequences ${sample_id}_rep-seqs-dada2_${params.suffix}.qza \
        --o-table ${sample_id}_table-dada2_${params.suffix}.qza \
        --o-denoising-stats ${sample_id}_stats-dada2_${params.suffix}.qza \
        --verbose
    """
}

/*
* DADA2 metrics (QC)
* Input   : feature table and representative sequences (.qza)
* Output  : three types of visualisation (qzv)
* Purpose : generate quality plots (QC analysis)
*/
process qc_dada2_meta {
    label 'qiime'
    publishDir params.result, mode: 'copy'

    input:
        tuple val(sample_id), path(stats_dada2)

    output:
        tuple val(sample_id), path("${sample_id}_stats-dada2_${params.suffix}.qzv")

    script:
    """
    qiime metadata tabulate \
        --m-input-file ${stats_dada2} \
        --o-visualization ${sample_id}_stats-dada2_${params.suffix}.qzv
    """
}

process qc_dada2_table {
    label 'qiime'
    publishDir params.result, mode: 'copy'

    input:
        tuple val(sample_id), path(table_dada2)

    output:
        tuple val(sample_id), path("${sample_id}_table-dada2_${params.suffix}.qzv")

    script:
    """
    qiime feature-table summarize \
        --i-table ${table_dada2} \
        --o-visualization ${sample_id}_table-dada2_${params.suffix}.qzv
    """
}

process qc_dada2_rep {
    label 'qiime'
    publishDir params.result, mode: 'copy'

    input:
        tuple val(sample_id), path(rep_dada2)

    output:
        tuple val(sample_id), path("${sample_id}_rep-seqs-dada2_${params.suffix}.qzv")

    script:
    """
    qiime feature-table tabulate-seqs \
        --i-data ${rep_dada2} \
        --o-visualization ${sample_id}_rep-seqs-dada2_${params.suffix}.qzv
    """
}

// -----------------------------------------------------------------------------
/*
* Taxonomic classification
* Input   : representative sequences
* Output  : taxonomy artifact
* Purpose : assign taxonomy using pre-trained classifier
*/
process taxa_classification {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        path classifier
        tuple val(sample_id), path(rep_dada2)

    output:
        tuple val(sample_id), path("${sample_id}_taxonomy_${params.suffix}.qza")

    script:
    """
    qiime feature-classifier classify-sklearn \
        --p-n-jobs ${params.n_jobs} \
        --p-confidence ${params.confidence} \
        --i-classifier ${classifier} \
        --i-reads ${rep_dada2} \
        --o-classification "${sample_id}_taxonomy_${params.suffix}.qza"
    """
}

/*
* Taxonomic classification filtering
* Input   : taxonomy artifact
* Output  : taxonomy artifact
* Purpose : focus on the taxa of interest
*/
process taxa_filtering {
    label 'qiime'
    publishDir "${params.result}/dev", mode: 'copy'

    input:
        tuple val(sample_id), path(taxa_classified), path(table_dada2)

    output:
       tuple val(sample_id), path(taxa_classified), path("${sample_id}_filtered-table_${params.suffix}.qza")

    script:
    """
    qiime taxa filter-table \
        --i-table ${table_dada2} \
        --i-taxonomy ${taxa_classified} \
        --p-include _ \
        --o-filtered-table "${sample_id}_filtered-table_${params.suffix}.qza"
    """
}

/*
* Taxonomic classification metrics
* Input   : taxonomy artifact
* Output  : visualisation (qzv), barplot or krona
* Purpose : generate quality plots (QC analysis)
*/
process qc_init_classification {
    label 'qiime'
    publishDir params.result, mode: 'copy'

    input:
        tuple val(sample_id), path(taxa_classified), path(table_dada2)

    output:
        tuple val(sample_id), 
            path(taxa_classified), 
            path("${sample_id}_rartaxa_barplot_${params.suffix}.qzv")

    script:
    """
    qiime taxa barplot \
        --i-table ${table_dada2} \
        --i-taxonomy ${taxa_classified} \
        --o-visualization "${sample_id}_rartaxa_barplot_${params.suffix}.qzv"
    """
}

process qc_filtered_classification {
    label 'qiime'
    publishDir params.result, mode: 'copy'

    input:
        tuple val(sample_id), path(taxa_classified), path(table_filtered)

    output:
        tuple val(sample_id), 
            path(taxa_classified), 
            path("${sample_id}_filtered_barplot_${params.suffix}.qzv")

    script:
    """
    qiime taxa barplot \
        --i-table ${table_filtered} \
        --i-taxonomy ${taxa_classified} \
        --o-visualization "${sample_id}_filtered_barplot_${params.suffix}.qzv"
    """
}

process krona_init_classification {
    label 'qiime'
    publishDir params.result, mode: 'copy'

    input:
        tuple val(sample_id), path(taxa_classified), path(table_dada2)

    output:
        tuple val(sample_id), 
            path(taxa_classified), 
            path("${sample_id}_rartaxa_krona_${params.suffix}.qzv")

    script:
    """
    qiime krona collapse-and-plot \
        --i-table ${table_dada2} \
        --i-taxonomy ${taxa_classified} \
        --o-krona-plot "${sample_id}_rartaxa_krona_${params.suffix}.qzv"
    """
}

process krona_filtered_classification {
    label 'qiime'
    publishDir params.result, mode: 'copy'

    input:
        tuple val(sample_id), path(taxa_classified), path(table_filtered)

    output:
        tuple val(sample_id), 
            path(taxa_classified), 
            path("${sample_id}_filtered_krona_${params.suffix}.qzv")

    script:
    """
    qiime krona collapse-and-plot \
        --i-table ${table_filtered} \
        --i-taxonomy ${taxa_classified} \
        --o-krona-plot "${sample_id}_filtered_krona_${params.suffix}.qzv"
    """
}

// -----------------------------------------------------------------------------
/*
* Information about Software used for analysis
* Input   : none
* Output  : software version (txt)
* Purpose : save the information about software version
*/
process create_info {
    label 'qiime'

    input:
        val input_dir
        val result_dir
        val suffix

        val paired_end
        val all_in_one
        val trimming

        val min_quality
        val min_length

        val trim_left_f
        val trim_left_r
        val trunc_len_f
        val trunc_len_r
        val n_threads
        val reads_learn
        val fold_parents

        val db
        val reads
        val taxa

        val confidence
        val n_jobs

    output:
        path "pipeline_${params.suffix}.txt"

    script:
    """
    sofrware_track_file="pipeline_${params.suffix}.txt"

    echo "QIIME2 - AMPLICONS ANALYSIS CONFIGURATION" > \$sofrware_track_file
    echo "" >> \$sofrware_track_file

    echo "Generated: \$(date '+%d/%m/%Y %H:%M:%S')" >> \$sofrware_track_file
    echo "" >> \$sofrware_track_file

    # -------------------------
    # General settings
    # -------------------------
    echo "GENERAL SETTINGS" >> \$sofrware_track_file
    echo "Input folder  : ${input_dir}" >> \$sofrware_track_file
    echo "Output folder : ${result_dir}" >> \$sofrware_track_file
    echo "Suffix        : ${suffix}" >> \$sofrware_track_file
    echo "" >> \$sofrware_track_file

    # -------------------------
    # Data type
    # -------------------------
    echo "ANALYSIS STRATEGY" >> \$sofrware_track_file

    if [ "${paired_end}" = true ]; then
        echo "Sequencing type : Paired-end (PE)" >> \$sofrware_track_file
    else
        echo "Sequencing type : Single-end (SE)" >> \$sofrware_track_file
    fi

    if [ "${all_in_one}" = true ]; then
        echo "Sample handling : All samples processed together" >> \$sofrware_track_file
    else
        echo "Sample handling : Samples processed separately" >> \$sofrware_track_file
    fi

    if [ "${trimming}" = true ]; then
        echo "Trimming        : Enabled" >> \$sofrware_track_file
    else
        echo "Trimming        : Disabled" >> \$sofrware_track_file
    fi

    echo "" >> \$sofrware_track_file

    # -------------------------
    # Trimming parameters
    # -------------------------
    echo "CUTADAPT DATA TRIMMING" >> \$sofrware_track_file
    echo "Phred Score Qual. : ${min_quality}" >> \$sofrware_track_file
    echo "Length min        : ${min_length}" >> \$sofrware_track_file

    echo "" >> \$sofrware_track_file

    # -------------------------
    # DADA2 parameters
    # -------------------------
    echo "DADA2 DENOISING PARAMETERS" >> \$sofrware_track_file
    echo "Trim left forward : ${trim_left_f} (not used if 0)" >> \$sofrware_track_file
    echo "Trim left reverse : ${trim_left_r} (not used if 0)" >> \$sofrware_track_file
    echo "Trunc length F    : ${trunc_len_f}" >> \$sofrware_track_file
    echo "Trunc length R    : ${trunc_len_r}" >> \$sofrware_track_file
    echo "Threads           : ${n_threads}" >> \$sofrware_track_file
    echo "Reads for model   : ${reads_learn}" >> \$sofrware_track_file
    echo "Fold parents      : ${fold_parents}" >> \$sofrware_track_file

    echo "" >> \$sofrware_track_file

    # -------------------------
    # Classifier training
    # -------------------------
    echo "CLASSIFIER TRAINING" >> \$sofrware_track_file
    echo "Database          : ${db}" >> \$sofrware_track_file
    echo "Reference reads   : ${reads}" >> \$sofrware_track_file
    echo "Taxonomy file     : ${taxa}" >> \$sofrware_track_file

    echo "" >> \$sofrware_track_file

    # -------------------------
    # Taxonomic classification
    # -------------------------
    echo "TAXONOMIC CLASSIFICATION" >> \$sofrware_track_file
    echo "Confidence threshold : ${confidence}" >> \$sofrware_track_file
    echo "Number of jobs       : ${n_jobs}" >> \$sofrware_track_file

    echo "" >> \$sofrware_track_file
    echo "CONFIGURATION COMPLETE" >> \$sofrware_track_file

    echo "" >> \$sofrware_track_file
    echo "--------------------------------------------------------------------------------" \
        >> \$sofrware_track_file

    echo "SOFTWARES VERSION" >> \$sofrware_track_file
    echo "" >> \$sofrware_track_file
    """
}

process qiime_info {
    label 'qiime'

    input:
        path file 

    output: 
        path "qiime_${params.suffix}.txt"

    script:
    """
    sofrware_track_file="qiime_${params.suffix}.txt"
    cat $file > \$sofrware_track_file

    echo "QIIME2" >> \$sofrware_track_file
    qiime info >> \$sofrware_track_file || true

    echo "" >> \$sofrware_track_file

    echo "KRONA VERSION" >> \$sofrware_track_file
    qiime krona --version >> \$sofrware_track_file
    """
}

process cutadapt_info {
    label 'cutadapt'

    input:
        path file 

    output: 
        path "cutadapt_${params.suffix}.txt"

    script:
    """
    sofrware_track_file="cutadapt_${params.suffix}.txt"
    cat $file > \$sofrware_track_file

    if [ "${params.trimming}" = "true" ]; then
        echo "" >> \$sofrware_track_file
        echo "CUTADAPT" >> \$sofrware_track_file
        cutadapt --version >> \$sofrware_track_file
    fi
    """
}

process publish_info {
    label 'qiime'
    publishDir "${params.result}", mode: 'copy'

    input:
        path file 

    output: 
        path "config_${params.suffix}.txt"

    script:
    """
    sofrware_track_file="config_${params.suffix}.txt"
    cat $file > \$sofrware_track_file
    """
}