#!/usr/bin/env nextflow
    
// enable dsl2
nextflow.enable.dsl=2


// -----------------------------------------------------------------------------
// MAIN WORFLOW FOR QIIME2 AMPLICONS ANALYSIS
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// INPUT CHANNEL
// -----------------------------------------------------------------------------

// collect fastq files in tuple [sample_id, R1, R2] or [sample_id, R1]
if (params.paired_end) {
    inputs_ch = Channel
        .fromFilePairs("${params.input_dir}/*_{R1,R2}{.fastq*,_*.fastq*}")
        .map { id, reads ->

            // Sample ID = everything before the first '_'
            def sample_id = reads[0].baseName.split('_')[0]

            assert sample_id
            assert reads.size() == 2

            tuple(sample_id, reads[0], reads[1])
        }
} else {
    inputs_ch = Channel
        .fromPath("${params.input_dir}/*_R1{.fastq*,_*.fastq*}")
        .map { r1 ->

            // Sample ID = everything before the first '_'
            def sample_id = r1.baseName.split('_')[0]

            assert sample_id

            tuple(sample_id, r1, null)
        }
}


// get files for classifier training
reads_ch = Channel.value(file(params.reads))
taxa_ch = Channel.value(file(params.taxa))


// -----------------------------------------------------------------------------
// INCLUDE MODULES
// -----------------------------------------------------------------------------

include {
    TRIM_FASTP
    GENERATE_MANIFEST
    GENERATE_MANIFEST_ALL
    IMPORT_MANIFEST
    ASSIGN_KRAKEN2
    MPA_MODIF
    MPA_TO_KRONA
    COUNT_FASTQ_READS
    MPA_FAMILY_BARPLOT
    QC_DEMUX
    DENOISE_DADA2
    QC_DADA2_META
    QC_DADA2_TABLE
    QC_DADA2_REP
    IMPORT_REFSEQ
    IMPORT_TAXA
    GENERATE_CLASSIFIER_BAYES
    SKLEARN_CLASSIFIER
    BLAST_CLASSIFIER
    VSEARCH_CLASSIFIER
    TAXA_FILTERING
    CREATE_INFO
    FASTQC_INFO
    MULTIQC_INFO
    SEQKIT_INFO
    KRAKEN2_INFO
    PYTHON_INFO
    KRONA_INFO
    QIIME_INFO
    FASTP_INFO
    PUBLISH_INFO
} from "./modules/modules_qiime2_amplicons.nf"

include { 
    QC_FASTQC as QC_FASTQC_RAW 
    QC_MULTIQC as QC_MULTIQC_RAW
    QC_CLASSIFICATION as QC_INIT_CLASSIFICATION
    KRONA_TAXA_LEVEL as KRONA_INIT_LEVEL
    KRONA_CLASSIFICATION as KRONA_INIT_CLASSIFICATION
    KRONA_TO_HTML as KRONA_INIT_TO_HTML
} from './modules/modules_qiime2_amplicons.nf'

include { 
    QC_FASTQC as QC_FASTQC_TRIM 
    QC_MULTIQC as QC_MULTIQC_TRIM
    QC_CLASSIFICATION as QC_FILT_CLASSIFICATION
    KRONA_TAXA_LEVEL as KRONA_FILT_LEVEL
    KRONA_CLASSIFICATION as KRONA_FILT_CLASSIFICATION
    KRONA_TO_HTML as KRONA_FILT_TO_HTML
} from './modules/modules_qiime2_amplicons.nf'


// -----------------------------------------------------------------------------
// WORKFLOW
// -----------------------------------------------------------------------------

workflow {

    // ---------------------------
    // IMPORT CLASSIFIER
    // ---------------------------
    classifier_path = file("${params.path_db}/${params.db}_classifier.qza")
    def use_reference = params.classifier in ['blast', 'vsearch']

    if (use_reference) {
        IMPORT_REFSEQ(reads_ch)
        IMPORT_TAXA(taxa_ch)
        refseq_ch = IMPORT_REFSEQ.out
        taxa_out_ch = IMPORT_TAXA.out

    } else if (classifier_path?.exists()) {
        classifier_ch = Channel.value(file(classifier_path))
        log.info "Classifier already in files, skipping training steps"

    } else {
        IMPORT_REFSEQ(reads_ch)
        IMPORT_TAXA(taxa_ch)
        refseq_ch = IMPORT_REFSEQ.out
        taxa_out_ch = IMPORT_TAXA.out

        classifier_ch = GENERATE_CLASSIFIER_BAYES(refseq_ch, taxa_out_ch)
    }


    // ---------------------------
    // raw QC PLOTS
    // ---------------------------
    read_type_raw = "0_Raw"

    qc_raw = QC_FASTQC_RAW(read_type_raw, inputs_ch)

    raw_fastqc_zips = qc_raw.zip_files
        .map { sample_id, zip -> zip }
        .flatten()
        .collect()

    QC_MULTIQC_RAW(read_type_raw, raw_fastqc_zips)


    // ---------------------------
    // TRIMMING
    // ---------------------------
    samples_ch = TRIM_FASTP(inputs_ch)


    // ---------------------------
    // trim QC PLOTS
    // ---------------------------
    read_type_trim = "0-1_Trimmed"

    qc_trimmed = QC_FASTQC_TRIM(read_type_trim, samples_ch)

    trimmed_fastqc_zips = qc_trimmed.zip_files
        .map { sample_id, zip -> zip }
        .flatten()
        .collect()

    QC_MULTIQC_TRIM(read_type_trim, trimmed_fastqc_zips)


    // ---------------------------
    // ORGANISM IDENTIFICATION - Kraken2
    // ---------------------------
    ASSIGN_KRAKEN2(samples_ch)
    MPA_MODIF(ASSIGN_KRAKEN2.out)
    MPA_TO_KRONA(MPA_MODIF.out)

    COUNT_FASTQ_READS(samples_ch)
    joined_mpamodif_total_ch = MPA_MODIF.out.join(COUNT_FASTQ_READS.out.kraken)
    MPA_FAMILY_BARPLOT(joined_mpamodif_total_ch)


    // ---------------------------
    // COUNT READS + FILTER
    // ---------------------------
    counted_parsed_ch = COUNT_FASTQ_READS.out.qiime
        .map { sample_id, file_nb, r1, r2 ->
            def nb_reads = file_nb.text.trim().toInteger()
            
            def learn_reads = Math.min(
                nb_reads,
                params.reads_learn
            )

            tuple(sample_id, learn_reads, nb_reads, r1, r2)
        }

    branched = counted_parsed_ch.branch {
        sufficient: it[1] >= params.min_reads
        insufficient: it[1] < params.min_reads
    }

    // Samples OK, enough reads for DADA2
    samples_ok_ch = branched.sufficient.map { sample_id, nb_reads, learn_reads, r1, r2  ->
        tuple(sample_id, r1, r2, learn_reads)
    }

    // Failed samples, not enough reads for DADA2
    failed_trim_ch = branched.insufficient.map { sample_id, nb_reads, learn_reads, r1, r2  ->
        "${sample_id}\t${nb_reads}\tLOW_READ_COUNT\n"
    }

    failed_trim_ch
        .collectFile(
            name: "Failed_TRIM.tsv",
            storeDir: "${params.result}/LOGS"
        )


    // ---------------------------
    // QIIME2 MANIFEST
    // ---------------------------
    GENERATE_MANIFEST(samples_ok_ch)

    def manifests_list = GENERATE_MANIFEST.out.map { it[1] }.collect()
    def min_reads_learn = params.reads_learn

    def manifests_ch = params.all_in_one \
        ? GENERATE_MANIFEST_ALL(manifests_list, min_reads_learn) \
        : GENERATE_MANIFEST.out


    // ---------------------------
    // IMPORT DATA
    // ---------------------------
    IMPORT_MANIFEST(manifests_ch)
    QC_DEMUX(IMPORT_MANIFEST.out)


    // ---------------------------
    // DADA2
    // ---------------------------
    DENOISE_DADA2(IMPORT_MANIFEST.out)
    QC_DADA2_META(DENOISE_DADA2.out.stats_dada2)
    QC_DADA2_TABLE(DENOISE_DADA2.out.table_dada2)
    QC_DADA2_REP(DENOISE_DADA2.out.rep_dada2)


    // ---------------------------
    // COUNT DADA2 SEQS + FILTER
    // ---------------------------
    rep_seqs_ch = DENOISE_DADA2.out.rep_dada2

    rep_checked_ch = rep_seqs_ch.map { sample_id, rep_file ->
        def size_ok = file(rep_file).exists() && file(rep_file).size() > 0
        tuple(sample_id, rep_file, size_ok)
    }

    rep_branched = rep_checked_ch.branch {
        sufficient: it[2] == true
        insufficient: it[2] == false
    }

    // files with DADA2 results
    rep_ok_ch = rep_branched.sufficient.map { sample_id, rep_file, _ ->
        tuple(sample_id, rep_file)
    }

    // files without any DADA2 results
    failed_dada2_ch = rep_branched.insufficient.map { sample_id, rep_file, _ ->
        "${sample_id}\tNO_PASS\n"
    }

    failed_dada2_ch
        .collectFile(
            name: "Failed_DADA2.tsv",
            storeDir: "${params.result}/LOGS"
        )


    // ---------------------------
    // CLASSIFICATION
    // ---------------------------
    if (params.classifier == 'blast') {
        BLAST_CLASSIFIER(refseq_ch, taxa_out_ch, rep_ok_ch)
        taxa_classified_ch = BLAST_CLASSIFIER.out.taxonomy
    
    } else if (params.classifier == 'sklearn') {
        SKLEARN_CLASSIFIER(classifier_ch, rep_ok_ch)
        taxa_classified_ch = SKLEARN_CLASSIFIER.out.taxonomy

    } else if (params.classifier == 'vsearch') {
        VSEARCH_CLASSIFIER(refseq_ch, taxa_out_ch, rep_ok_ch)
        taxa_classified_ch = VSEARCH_CLASSIFIER.out.taxonomy
    }


    // ---------------------------
    // COUNT RESULTS + FILTER
    // ---------------------------
    // files with taxa_classified results
    taxa_classified_ok_ch = taxa_classified_ch
        .filter { sample_id, init_qza, status ->
                file(status).text.trim() == "OK"
            }
        .map { sample_id, init_qza, status -> tuple(sample_id, init_qza) }

    // files with no taxa_classified results
    taxa_classified_failed_ch = taxa_classified_ch
        .filter { sample_id, init_qza, status ->
            file(status).text.trim() != "OK"
        }
        .map { sample_id, init_qza, status ->
            "${sample_id}\t${file(status).text.trim()}\n"
        }

    taxa_classified_failed_ch
        .collectFile(
            name: "Failed_CLASSIFY.tsv",
            storeDir: "${params.result}/LOGS"
        )


    // ---------------------------
    // INIT QC
    // ---------------------------
    results_type = "init"

    KRONA_INIT_LEVEL(taxa_classified_ok_ch)

    joined_taxa_table_ch = taxa_classified_ok_ch.join(DENOISE_DADA2.out.table_dada2)
    joined_taxa_table_max_ch = joined_taxa_table_ch.join(KRONA_INIT_LEVEL.out)

    QC_INIT_CLASSIFICATION(results_type, joined_taxa_table_ch)
    KRONA_INIT_CLASSIFICATION(results_type, joined_taxa_table_max_ch)
    KRONA_INIT_TO_HTML(results_type, KRONA_INIT_CLASSIFICATION.out)


    // ---------------------------
    // FILTERING
    // ---------------------------
    TAXA_FILTERING(joined_taxa_table_ch)


    // ---------------------------
    // FILTERING STATUS CHECK
    // ---------------------------
    taxa_filtered_ok_ch = TAXA_FILTERING.out
        .filter { sample_id, taxa_classified, filtered_qza, status ->
            file(status).text.trim() == "OK"
        }
        .map { sample_id, taxa_classified, filtered_qza, status ->
            tuple(sample_id, taxa_classified, filtered_qza)
        }

    taxa_filtered_failed_ch = TAXA_FILTERING.out
        .filter { sample_id, taxa_classified, filtered_qza, status ->
            file(status).text.trim() != "OK"
        }
        .map { sample_id, taxa_classified, filtered_qza, status ->
            "${sample_id}\t${file(status).text.trim()}\n"
        }

    taxa_filtered_failed_ch
        .collectFile(
            name: "Failed_TAXAFILTER.tsv",
            storeDir: "${params.result}/LOGS"
        )
    

    // ---------------------------
    // FILTERING QC
    // ---------------------------
    results_type = "filt"

    taxa_filtered_ok_krona_ch = taxa_filtered_ok_ch
    .map { sample_id, taxa_classified, filtered_qza ->
        tuple(sample_id, taxa_classified)
    }

    KRONA_FILT_LEVEL(taxa_filtered_ok_krona_ch)

    joined_filt_table_max_ch = taxa_filtered_ok_ch.join(KRONA_FILT_LEVEL.out)

    QC_FILT_CLASSIFICATION(results_type, taxa_filtered_ok_ch)
    KRONA_FILT_CLASSIFICATION(results_type, joined_filt_table_max_ch)
    KRONA_FILT_TO_HTML(results_type, KRONA_FILT_CLASSIFICATION.out)

    // ---------------------------
    // TRACKING CONFIG
    // ---------------------------
    CREATE_INFO(
        params.input_dir,
        params.result,
        params.suffix,

        params.paired_end,
        params.all_in_one,
        params.adapters,

        params.min_quality,
        params.min_length,

        params.trim_left_f,
        params.trim_left_r,
        params.trunc_len_f,
        params.trunc_len_r,
        params.reads_learn,
        params.fold_parents,

        "${params.path_db}/${params.db}_classifier.qza",
        params.reads,
        params.taxa,

        params.sklearn_confidence,
        params.blast_identity,
        params.blast_maxaccepts,
        params.blast_query_cov,
        params.vsearch_identity,
        params.vsearch_maxaccepts,
        params.vsearch_query_cov,
        params.classifier,

        params.kraken2_db
    )
    FASTQC_INFO(CREATE_INFO.out)
    MULTIQC_INFO(FASTQC_INFO.out)
    SEQKIT_INFO(MULTIQC_INFO.out)
    KRAKEN2_INFO(SEQKIT_INFO.out)
    PYTHON_INFO(KRAKEN2_INFO.out)
    KRONA_INFO(PYTHON_INFO.out)
    FASTP_INFO(KRONA_INFO.out)
    QIIME_INFO(FASTP_INFO.out)
    PUBLISH_INFO(QIIME_INFO.out)
}