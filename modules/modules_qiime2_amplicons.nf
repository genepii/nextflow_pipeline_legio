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
        path(reads_file)

    output:
        path("${params.db}.qza")

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
        path(taxa_file)

    output:
        path("${params.db}_tax.qza")

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
        path(classifier_reads)
        path(classifier_taxa)

    output:
        path("${params.db}_classifier.qza")

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
        path(fastqc_zip)

    output:
        path("General_multiQC_report.html")

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
    def adapters = params.adapters ? 
        "--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" : 
        "--disable_adapter_trimming"
    
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
        tuple val(sample_id), val(r1), val(r2), val(reads_learn)

    output:
        tuple val(sample_id), path("${sample_id}_manifest.tsv"), val(reads_learn)

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
        path(manifests)
        val(min_reads_learn)

    output:
        tuple val("All-samples"), 
            path("All-samples_manifest.tsv"), 
            val(min_reads_learn)

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
        tuple val(sample_id), path(manifest), val(reads_learn)

    output:
        tuple val(sample_id), path("${sample_id}_demux.qza"), val(reads_learn)

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
        tuple val(sample_id), path(demux), val(reads_learn)

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
* Total sequences count recovery from trimmed FASTQ
* Input   : trimmed R1 FASTQ
* Output  : number of reads
* Purpose : stop the analysis for sample with not enough reads inside
*/
process COUNT_READS {
    label 'qiime'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_totalseq.txt"), path(r1), path(r2)

    script:
    """
    nb_reads=\$(zcat ${r1} | wc -l)
    nb_reads=\$((nb_reads / 4))
    echo \$nb_reads > "${sample_id}_totalseq.txt"
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
    errorStrategy 'ignore'

    input:
        tuple val(sample_id), path(demux), val(reads_learn)

    output:
        tuple val(sample_id), 
            path("${sample_id}_table-dada2.qza"), 
            emit: table_dada2
        tuple val(sample_id), 
            path("${sample_id}_stats-dada2.qza"), 
            emit: stats_dada2
        tuple val(sample_id), 
            path("${sample_id}_rep-seqs-dada2.qza"), 
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
        --p-n-reads-learn ${reads_learn} \
        --p-n-threads ${task.cpus} \
        --o-representative-sequences ${sample_id}_rep-seqs-dada2.qza \
        --o-table ${sample_id}_table-dada2.qza \
        --o-denoising-stats ${sample_id}_stats-dada2.qza \
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
        tuple val(sample_id), path("${sample_id}_stats-dada2.qzv")

    script:
    """
    qiime metadata tabulate \
        --m-input-file ${stats_dada2} \
        --o-visualization ${sample_id}_stats-dada2.qzv
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
        tuple val(sample_id), path("${sample_id}_table-dada2.qzv")

    script:
    """
    qiime feature-table summarize \
        --i-table ${table_dada2} \
        --o-visualization ${sample_id}_table-dada2.qzv
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
        tuple val(sample_id), path("${sample_id}_rep-seqs-dada2.qzv")

    script:
    """
    qiime feature-table tabulate-seqs \
        --i-data ${rep_dada2} \
        --o-visualization ${sample_id}_rep-seqs-dada2.qzv
    """
}

// -----------------------------------------------------------------------------
/*
* Taxonomic classification with Bayes classifier
* Input   : representative sequences
* Output  : taxonomy artifact + status file
* Purpose : assign taxonomy using pre-trained classifier
*/
process SKLEARN_CLASSIFIER {
    label 'qiime'
    publishDir "${params.result}/dev/3_Classification", mode: 'copy'

    input:
        path(classifier)
        tuple val(sample_id), path(rep_dada2)

    output:
        tuple val(sample_id), 
            path("${sample_id}_taxonomySklearn.qza"), 
            path("${sample_id}_statusSklearn.txt"), 
            emit: taxonomy

    script:
    """
    set +e

    qiime feature-classifier classify-sklearn \
        --p-n-jobs ${task.cpus} \
        --p-confidence ${params.sklearn_confidence} \
        --i-classifier ${classifier} \
        --i-reads ${rep_dada2} \
        --o-classification "${sample_id}_taxonomySklearn.qza"
        
    status=\$?

    if [ \$status -ne 0 ]; then
        echo "TECH_FAIL" > ${sample_id}_statusSklearn.txt
        touch "${sample_id}_taxonomySklearn.qza"
        exit 0
    fi

    # Content check
    qiime tools export \
        --input-path "${sample_id}_taxonomySklearn.qza" \
        --output-path tmp_export

    if [ ! -s tmp_export/taxonomy.tsv ] || [ \$(wc -l < tmp_export/taxonomy.tsv) -le 1 ]; then
        echo "NO_HIT" > ${sample_id}_statusSklearn.txt
    else
        echo "OK" > ${sample_id}_statusSklearn.txt
    fi
    """
}

/*
* Taxonomic classification with Blast
* Input   : representative sequences + reference database
* Output  : taxonomy assignment (qza) + status file
* Purpose : assign taxonomy using Blast consensus classifier
*/
process BLAST_CLASSIFIER {
    label 'qiime'
    publishDir "${params.result}/dev/3_Classification", mode: 'copy'

    input:
        path(reads)
        path(taxa)
        tuple val(sample_id), path(rep_dada2)

    output:
        tuple val(sample_id), 
            path("${sample_id}_taxonomyBlast.qza"), 
            path("${sample_id}_statusBlast.txt"), 
            emit: taxonomy
        path("${sample_id}_resultsBlast.qza")

    script:
    """
    set +e

    qiime feature-classifier classify-consensus-blast \
        --i-query ${rep_dada2} \
        --i-reference-reads ${reads} \
        --i-reference-taxonomy ${taxa} \
        --p-perc-identity ${params.blast_identity} \
        --p-maxaccepts ${params.blast_maxaccepts} \
        --p-query-cov ${params.blast_query_cov} \
        --p-strand both \
        --p-num-threads ${task.cpus} \
        --o-classification "${sample_id}_taxonomyBlast.qza" \
        --o-search-results "${sample_id}_resultsBlast.qza"

    status=\$?

    if [ \$status -ne 0 ]; then
        echo "TECH_FAIL" > ${sample_id}_statusBlast.txt
        touch "${sample_id}_taxonomyBlast.qza"
        touch "${sample_id}_resultsBlast.qza"
        exit 0
    fi

    # Content check
    qiime tools export \
        --input-path "${sample_id}_taxonomyBlast.qza" \
        --output-path tmp_export

    if [ ! -s tmp_export/taxonomy.tsv ] || [ \$(wc -l < tmp_export/taxonomy.tsv) -le 1 ]; then
        echo "NO_HIT" > ${sample_id}_statusBlast.txt
    else
        echo "OK" > ${sample_id}_statusBlast.txt
    fi
    """
}

/*
* Taxonomic classification with Vsearch
* Input   : representative sequences + reference database
* Output  : taxonomy assignment (qza) + status file
* Purpose : assign taxonomy using vsearch consensus classifier
*/
process VSEARCH_CLASSIFIER {
    label 'qiime'
    publishDir "${params.result}/dev/3_Classification", mode: 'copy'

    input:
        path(reads)
        path(taxa)
        tuple val(sample_id), path(rep_dada2)

    output:
        tuple val(sample_id),
            path("${sample_id}_taxonomyVsearch.qza"),
            path("${sample_id}_statusVsearch.txt"), 
            emit: taxonomy
        path("${sample_id}_resultsVsearch.qza")

    script:
    """
    set +e

    qiime feature-classifier classify-consensus-vsearch \
        --i-query ${rep_dada2} \
        --i-reference-reads ${reads} \
        --i-reference-taxonomy ${taxa} \
        --p-perc-identity ${params.vsearch_identity} \
        --p-maxaccepts ${params.vsearch_maxaccepts} \
        --p-query-cov ${params.vsearch_query_cov} \
        --p-strand both \
        --p-threads ${task.cpus} \
        --o-classification "${sample_id}_taxonomyVsearch.qza" \
        --o-search-results "${sample_id}_resultsVsearch.qza"

    status=\$?

    if [ \$status -ne 0 ]; then
        echo "TECH_FAIL" > ${sample_id}_statusVsearch.txt
        touch "${sample_id}_taxonomyVsearch.qza"
        touch "${sample_id}_resultsVsearch.qza"
        exit 0
    fi

    # Export check
    qiime tools export \
        --input-path "${sample_id}_taxonomyVsearch.qza" \
        --output-path tmp_export

    if [ ! -s tmp_export/taxonomy.tsv ] || [ \$(wc -l < tmp_export/taxonomy.tsv) -le 1 ]; then
        echo "NO_HIT" > ${sample_id}_statusVsearch.txt
    else
        echo "OK" > ${sample_id}_statusVsearch.txt
    fi
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
        tuple val(sample_id),
            path(taxa_classified),
            path("${sample_id}_filtTable.qza"),
            path("${sample_id}_statusFilter.txt")

    script:
    """
    set +e

    qiime taxa filter-table \
        --i-table ${table_dada2} \
        --i-taxonomy ${taxa_classified} \
        --p-include '_' \
        --o-filtered-table ${sample_id}_filtTable.qza

    status=\$?

    if [ \$status -ne 0 ]; then
        echo "EMPTY" > ${sample_id}_statusFilter.txt
        touch "${sample_id}_filtTable.qza"
        exit 0
    else
        echo "OK" > ${sample_id}_statusFilter.txt
    fi
    """
}

/*
* Taxonomic classification overview
* Input   : DADA2 feature table + taxonomic assignments
* Output  : taxonomic barplot visualization (.qzv)
* Purpose : evaluate community composition
*/
process QC_CLASSIFICATION {
    label 'qiime'
    publishDir "${params.result}/3_Classification", mode: 'copy'

    input:
        val(results_type)
        tuple val(sample_id), path(taxa_classified), path(table_dada2)

    output:
        tuple val(sample_id), path("${sample_id}_${results_type}Barplot.qzv")

    script:
    """
    qiime taxa barplot \
        --i-table ${table_dada2} \
        --i-taxonomy ${taxa_classified} \
        --o-visualization "${sample_id}_${results_type}Barplot.qzv"
    """
}

/*
* Export QIIME2 taxonomy artifact for downstream processing
* Input   : taxonomic assignments
* Output  : max level taxonomy (.txt)
* Purpose : extract taxonomy and compute the maximum depth
*           to safely parameterize Krona (--p-level)
*/
process KRONA_TAXA_LEVEL{
    label 'qiime'

    input:
        tuple val(sample_id), path(taxa_classified)

    output:
        tuple val(sample_id), path("${sample_id}_maxLevel.txt")

    script:
    """
    set -euo pipefail

    EXPORT_DIR=${sample_id}_export
    mkdir -p \$EXPORT_DIR

    # Export taxonomy from QIIME2 artifact
    qiime tools export \
        --input-path ${taxa_classified} \
        --output-path \$EXPORT_DIR

    # Direct path (known structure for FeatureData[Taxonomy])
    TAX_FILE=\$EXPORT_DIR/taxonomy.tsv

    # Compute max depth
    awk -F '\\t' '
    NR > 1 {
        n = split(\$2, a, ";")
        count = 0
        for (i = 1; i <= n; i++) {
            if (a[i] != "") count++
        }
        if (count > max) max = count
    }
    END {
        if (max < 1) max = 1
        print max
    }' \$TAX_FILE > ${sample_id}_maxLevel.txt
    """
}

/*
* Taxonomic classification overview
* Input   : DADA2 feature table + taxonomic assignments
* Output  : interactive Krona plot (.qzv)
* Purpose : explore hierarchical taxonomic abundance
*/
process KRONA_CLASSIFICATION {
    label 'qiime'
    publishDir "${params.result}/3_Classification", mode: 'copy'

    input:
        val(results_type)
        tuple val(sample_id), path(taxa_classified), path(table_dada2), path(maxlevel)

    output:
        tuple val(sample_id), path("${sample_id}_${results_type}Krona.qzv")

    script:
    """
    max=\$(awk 'NR==1{print \$1}' ${maxlevel})

    qiime krona collapse-and-plot \
        --i-table ${table_dada2} \
        --i-taxonomy ${taxa_classified} \
        --o-krona-plot ${sample_id}_${results_type}Krona.qzv \
        --p-level \$max
    """
}

/*
* Taxonomic classification overview in HTML
* Input   : interactive Krona plot (.qzv)
* Output  : interactive Krona plot (.html)
* Purpose : explore hierarchical taxonomic abundance
*/
process KRONA_TO_HTML {
    label 'qiime'
    publishDir "${params.result}/3_Classification", mode: 'copy'

    input:
        val(results_type)
        tuple val(sample_id), path(krona_qzv)

    output:
        tuple val(sample_id),
            path("${sample_id}/${results_type}Krona/*")

    script:
    """
    qiime tools export \
        --input-path ${krona_qzv} \
        --output-path "${sample_id}/${results_type}Krona"
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
        val(input_dir)
        val(result_dir)
        val(suffix)

        val(paired_end)
        val(all_in_one)
        val(adapters)

        val(min_quality)
        val(min_length)

        val(trim_left_f)
        val(trim_left_r)
        val(trunc_len_f)
        val(trunc_len_r)
        val(reads_learn)
        val(fold_parents)

        val(db)
        val(reads)
        val(taxa)

        val(sklearn_confidence)
        val(blast_identity)
        val(blast_maxaccepts)
        val(blast_query_cov)
        val(vsearch_identity)
        val(vsearch_maxaccepts)
        val(vsearch_query_cov)
        val(classifier)

    output:
        path("pipeline_${suffix}.txt")

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
        "${reads_learn}" \
        "${fold_parents}" \
        "${db}" \
        "${reads}" \
        "${taxa}" \
        "${sklearn_confidence}" \
        "${blast_identity}" \
        "${blast_maxaccepts}" \
        "${blast_query_cov}" \
        "${vsearch_identity}" \
        "${vsearch_maxaccepts}" \
        "${vsearch_query_cov}" \
        "${classifier}"
    """
}

process FASTQC_INFO {
    label 'fastqc'

    input:
        path(file) 

    output: 
        path("fastqc_${params.suffix}.txt")

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
        path(file)

    output: 
        path("multiqc_${params.suffix}.txt")

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
        path(file)

    output: 
        path("qiime_${params.suffix}.txt")

    script:
    """
    software_track_file="qiime_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

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
        path(file)

    output: 
        path("fastp_${params.suffix}.txt")

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
        path(file)

    output: 
        path("softwaresTrackfile_${params.suffix}.txt")

    script:
    """
    software_track_file="softwaresTrackfile_${params.suffix}.txt"
    cat $file > \$software_track_file
    """
}