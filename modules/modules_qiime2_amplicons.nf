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
process IMPORT_REFSEQ {
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

process IMPORT_TAXA {
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

process GENERATE_CLASSIFIER_BAYES {
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
process TRIMMING_CUTADAPT {
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
        -j ${task.cpus} \
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
process GENERATE_MANIFEST {
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
process GENERATE_MANIFEST_ALL {
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
process IMPORT_MANIFEST {
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
process QC_DEMUX {
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
process DENOISE_DADA2 {
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
process QC_DADA2_META {
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

process QC_DADA2_TABLE {
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

process QC_DADA2_REP {
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
process TAXA_CLASSIFICATION {
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
        --p-n-jobs ${task.cpus} \
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
process TAXA_FILTERING {
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
process QC_INIT_CLASSIFICATION {
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

process QC_FILTERED_CLASSIFICATION {
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

process KRONA_INIT_CLASSIFICATION {
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

process KRONA_FILT_CLASSIFICATION {
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
* Input   : all params.values
* Output  : software version (txt)
* Purpose : save the information about software version
*/
process CREATE_INFO {
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
        path "pipeline_${suffix}.txt"

    script:
    """
    qiime2_amplicons_create_info.sh \
        "${input_dir}" \
        "${result_dir}" \
        "${suffix}" \
        "${paired_end}" \
        "${all_in_one}" \
        "${trimming}" \
        "${min_quality}" \
        "${min_length}" \
        "${trim_left_f}" \
        "${trim_left_r}" \
        "${trunc_len_f}" \
        "${trunc_len_r}" \
        "${n_threads}" \
        "${reads_learn}" \
        "${fold_parents}" \
        "${db}" \
        "${reads}" \
        "${taxa}" \
        "${confidence}" \
        "${n_jobs}"
    """
}

process QIIME_INFO {
    label 'qiime'

    input:
        path file 

    output: 
        path "qiime_${params.suffix}.txt"

    script:
    """
    software_track_file="qiime_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "QIIME2" >> \$software_track_file
    qiime info >> \$software_track_file || true

    echo "" >> \$software_track_file

    echo "KRONA VERSION" >> \$software_track_file
    qiime krona --version >> \$software_track_file
    """
}

process CUTADAPT_INFO {
    label 'cutadapt'

    input:
        path file 

    output: 
        path "cutadapt_${params.suffix}.txt"

    script:
    """
    software_track_file="cutadapt_${params.suffix}.txt"
    cat $file > \$software_track_file

    if [ "${params.trimming}" = "true" ]; then
        echo "" >> \$software_track_file
        echo "CUTADAPT" >> \$software_track_file
        cutadapt --version >> \$software_track_file
    fi
    """
}

process PUBLISH_INFO {
    label 'qiime'
    publishDir "${params.result}", mode: 'copy'

    input:
        path file 

    output: 
        path "config_${params.suffix}.txt"

    script:
    """
    software_track_file="config_${params.suffix}.txt"
    cat $file > \$software_track_file
    """
}