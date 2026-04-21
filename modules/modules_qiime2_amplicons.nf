#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Processes for workflow_qiime2_amplicons.nf


// -----------------------------------------------------------------------------
/*
* Custom 23S–5S reference database
* Input   : FASTA sequences and associated taxonomy file
* Output  : trained Naive Bayes classifier (QIIME2 artifact)
* Purpose : adapt the classifier to local dataset specificity
*/
process IMPORT_REFSEQ {
    label 'qiime'
    publishDir "${params.result}/dev/0_Classifier", mode: 'copy'

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

/*
* Taxonomy import step
* Input   : taxonomy file in TSV format
* Output  : QIIME2 taxonomy artifact (.qza)
* Purpose : prepare taxonomic annotations for classifier training
*/
process IMPORT_TAXA {
    label 'qiime'
    publishDir "${params.result}/dev/0_Classifier", mode: 'copy'

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

/*
* Naive Bayes classifier training
* Input   : reference sequences + taxonomy artifact
* Output  : trained QIIME2 classifier (.qza)
* Purpose : enable taxonomic assignment of ASVs/reads
*/
process GENERATE_CLASSIFIER_BAYES {
    label 'qiime'
    publishDir "${params.result}/dev/0_Classifier", mode: 'copy'

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
* Reads quality control
* Input   : FASTQ files (R1 or R1/R2)
* Output  : FastQC reports (HTML + ZIP)
* Purpose : assess sequencing quality before or after downstream processing
*/
process QC_FASTQC {
    label 'fastqc'
    publishDir "${params.result}/0_FastQC/${read_type}", mode: 'copy'

    input:
        val(read_type)
        tuple val(sample_id), val(r1), val(r2)

    output:
        tuple val(sample_id), path("*.zip"), emit: zip_files
        tuple val(sample_id), path("*.html"), emit: html_files

    script:
    def inputs = params.paired_end ?
        "${r1} ${r2}" :
        "${r1}"

    """
    fastqc \
        ${inputs} \
        --threads ${task.cpus} \
        --outdir ./       
    """
}

/*
* Aggregated quality control report
* Input   : collection of FastQC reports (ZIP files)
* Output  : MultiQC HTML report
* Purpose : provide a global overview of sequencing quality
*/
process QC_MULTIQC {
    label 'multiqc'
    publishDir "${params.result}/0_FastQC/${read_type}", mode: 'copy'

    input:
        val(read_type)
        path fastqc_zip

    output:
        path "General_multiQC_report.html"

    script:
    """
    multiqc ${fastqc_zip} \
        --filename "General_multiQC_report.html" \
        --force
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
        tuple val(sample_id), val(r1), val(r2)

    output:
        tuple val(sample_id),
            path("${sample_id}_trimR1.fastq.gz"),
            path("${sample_id}_trimR2.fastq.gz")
        
    script:
    def adapters = params.adapters ? "" : "--disable_adapter_trimming"
    
    def r1_out = "${sample_id}_trimR1.fastq.gz"
    def r2_out = "${sample_id}_trimR2.fastq.gz"
    def io_opts = params.paired_end ?
        "-i ${r1} -I ${r2} -o ${r1_out} -O ${r2_out}" :
        "-i ${r1} -o ${r1_out}"

    """
    fastp \
        ${io_opts} \
        --qualified_quality_phred ${params.min_quality} \
        --length_required ${params.min_length} \
        --thread ${task.cpus} \
        ${adapters}
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
    publishDir "${params.result}/dev/1_Qiime2", mode: 'copy'

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
    publishDir "${params.result}/dev/1_Qiime2", mode: 'copy'

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
    publishDir "${params.result}/dev/1_Qiime2", mode: 'copy'

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
    publishDir "${params.result}/1_Qiime2", mode: 'copy'

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
    publishDir "${params.result}/dev/2_Dada2", mode: 'copy'

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
* DADA2 denoising metadata QC report
* Input   : DADA2 summary statistics artifact
* Output  : QIIME2 visualization (.qzv) of denoising stats
* Purpose : assess read filtering, error correction, and denoising performance
*/
process QC_DADA2_META {
    label 'qiime'
    publishDir "${params.result}/2_Dada2", mode: 'copy'

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

/*
* DADA2 feature table summary QC
* Input   : feature table generated by DADA2
* Output  : QIIME2 visualization of feature table statistics
* Purpose : evaluate sequencing depth distribution and sample composition
*/
process QC_DADA2_TABLE {
    label 'qiime'
    publishDir "${params.result}/2_Dada2", mode: 'copy'

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

/*
* Representative sequences QC visualization
* Input   : representative sequences from DADA2
* Output  : QIIME2 visualization of sequence distribution
* Purpose : inspect sequence diversity and representative ASVs
*/
process QC_DADA2_REP {
    label 'qiime'
    publishDir "${params.result}/2_Dada2", mode: 'copy'

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
    publishDir "${params.result}/dev/3_Classification", mode: 'copy'

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
    publishDir "${params.result}/dev/3_Classification", mode: 'copy'

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
* Initial taxonomic classification overview (rarefied/unfiltered data)
* Input   : DADA2 feature table + taxonomic assignments
* Output  : taxonomic barplot visualization (.qzv)
* Purpose : evaluate initial community composition before filtering
*/
process QC_INIT_CLASSIFICATION {
    label 'qiime'
    publishDir "${params.result}/3_Classification", mode: 'copy'

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

/*
* Taxonomic classification after feature filtering
* Input   : filtered feature table + taxonomic assignments
* Output  : filtered taxonomic barplot visualization (.qzv)
* Purpose : assess how filtering impacts community structure representation
*/
process QC_FILTERED_CLASSIFICATION {
    label 'qiime'
    publishDir "${params.result}/3_Classification", mode: 'copy'

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

/*
* Initial taxonomic composition visualization using Krona (unfiltered data)
* Input   : DADA2 feature table + taxonomic assignments
* Output  : interactive Krona plot (.qzv)
* Purpose : explore hierarchical taxonomic abundance before filtering
*/
process KRONA_INIT_CLASSIFICATION {
    label 'qiime'
    publishDir "${params.result}/3_Classification", mode: 'copy'

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

/*
* Taxonomic composition visualization using Krona after filtering
* Input   : filtered feature table + taxonomic assignments
* Output  : interactive Krona plot (.qzv)
* Purpose : evaluate hierarchical taxonomy structure after quality filtering
*/
process KRONA_FILT_CLASSIFICATION {
    label 'qiime'
    publishDir "${params.result}/3_Classification", mode: 'copy'

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
        val adapters

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
        "${adapters}" \
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

process FASTQC_INFO {
    label 'fastqc'

    input:
        path file 

    output: 
        path "fastqc_${params.suffix}.txt"

    script:
    """
    software_track_file="fastqc_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "FASTQC VERSION" >> \$software_track_file
    fastqc --version >> \$software_track_file || true
    """
}

process MULTIQC_INFO {
    label 'multiqc'

    input:
        path file 

    output: 
        path "multiqc_${params.suffix}.txt"

    script:
    """
    software_track_file="multiqc_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "MULTIQC VERSION" >> \$software_track_file
    multiqc --version >> \$software_track_file || true
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

process FASTP_INFO {
    label 'fastp'

    input:
        path file 

    output: 
        path "fastp_${params.suffix}.txt"

    script:
    """
    software_track_file="fastp_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "FASTP VERSION" >> \$software_track_file
    fastp --version >> \$software_track_file || true
    """
}

process PUBLISH_INFO {
    label 'qiime'
    publishDir "${params.result}", mode: 'copy'

    input:
        path file 

    output: 
        path "softwaresTrackfile_${params.suffix}.txt"

    script:
    """
    software_track_file="softwaresTrackfile_${params.suffix}.txt"
    cat $file > \$software_track_file
    """
}