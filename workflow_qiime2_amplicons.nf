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
        .fromFilePairs("${params.input_dir}/*_{R1,R2}.fastq*")
        .map { id, reads ->
            def sample_id = id.replaceFirst(/_R1.*/, '')
            def r1 = reads[0]
            def r2 = reads[1]

            assert sample_id 
            assert r1 
            assert r2

            tuple(sample_id, r1, r2)
        }

} else {
    inputs_ch = Channel
        .fromPath("${params.input_dir}/*_R1.fastq*")
        .map { r1 ->
            def sample_id = r1.baseName.replaceFirst(/_R1.*/, '')

            assert sample_id 
            assert r1 

            tuple(sample_id, r1, null)
        }
}

// get files for classifier training
reads_ch = Channel.fromPath(params.reads)
taxa_ch = Channel.fromPath(params.taxa)


// -----------------------------------------------------------------------------
// INCLUDE MODULES
// -----------------------------------------------------------------------------

include {
    TRIMMING_CUTADAPT
    GENERATE_MANIFEST
    GENERATE_MANIFEST_ALL
    IMPORT_MANIFEST
    QC_DEMUX
    DENOISE_DADA2
    QC_DADA2_META
    QC_DADA2_TABLE
    QC_DADA2_REP
    IMPORT_REFSEQ
    IMPORT_TAXA
    GENERATE_CLASSIFIER_BAYES
    TAXA_CLASSIFICATION
    TAXA_FILTERING
    QC_INIT_CLASSIFICATION
    QC_FILTERED_CLASSIFICATION
    KRONA_INIT_CLASSIFICATION
    KRONA_FILT_CLASSIFICATION
    CREATE_INFO
    QIIME_INFO
    CUTADAPT_INFO
    PUBLISH_INFO
} from "./modules/modules_qiime2_amplicons.nf"


// -----------------------------------------------------------------------------
// WORKFLOW
// -----------------------------------------------------------------------------

workflow {

    // ---------------------------
    // BAYESIAN CLASSIFIER
    // ---------------------------
    IMPORT_REFSEQ(reads_ch)
    IMPORT_TAXA(taxa_ch)
    GENERATE_CLASSIFIER_BAYES(IMPORT_REFSEQ.out, IMPORT_TAXA.out)

    // ---------------------------
    // TRIMMING
    // ---------------------------
    def samples_ch = params.trimming \
        ? TRIMMING_CUTADAPT(inputs_ch) \
        : inputs_ch

    // ---------------------------
    // QIIME2 MANIFEST
    // ---------------------------
    GENERATE_MANIFEST(samples_ch)

    def manifests_ch = params.all_in_one \
        ? GENERATE_MANIFEST_ALL(GENERATE_MANIFEST.out.collect()) \
        : GENERATE_MANIFEST.out.collect()

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
    // CLASSIFICATION + QC
    // ---------------------------
    TAXA_CLASSIFICATION(GENERATE_CLASSIFIER_BAYES.out, DENOISE_DADA2.out.rep_dada2)

    joined_taxa_table_ch = TAXA_CLASSIFICATION.out.join(DENOISE_DADA2.out.table_dada2)
    QC_INIT_CLASSIFICATION(joined_taxa_table_ch)
    KRONA_INIT_CLASSIFICATION(joined_taxa_table_ch)

    // ---------------------------
    // FILTERING + QC
    // ---------------------------
    TAXA_FILTERING(joined_taxa_table_ch)
    QC_FILTERED_CLASSIFICATION(TAXA_FILTERING.out)
    KRONA_FILT_CLASSIFICATION(TAXA_FILTERING.out)

    // ---------------------------
    // TRACKING CONFIG
    // ---------------------------
    CREATE_INFO(
        params.input_dir,
        params.result,
        params.suffix,

        params.paired_end,
        params.all_in_one,
        params.trimming,

        params.min_quality,
        params.min_length,

        params.trim_left_f,
        params.trim_left_r,
        params.trunc_len_f,
        params.trunc_len_r,
        params.n_threads,
        params.reads_learn,
        params.fold_parents,

        params.db,
        params.reads,
        params.taxa,

        params.confidence,
        params.n_jobs
    )
    QIIME_INFO(CREATE_INFO.out)
    CUTADAPT_INFO(QIIME_INFO.out)
    PUBLISH_INFO(CUTADAPT_INFO.out)
}