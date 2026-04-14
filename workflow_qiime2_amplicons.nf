#!/usr/bin/env nextflow

// Test
// ./nextflow_25.10.4 workflow_qiime2_amplicons.nf -c config/qiime2_amplicons.config -profile test
    
// enable dsl2
nextflow.enable.dsl=2

/*
* Main workflow for Qiime2 Amplicons analysis
*/


// -----------------------------------------------------------------------------
// Channels
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
// Import modules
include {
    trimming_cutadapt
    generate_manifest
    generate_manifest_all
    import_manifest
    qc_demux
    denoise_dada2
    qc_dada2_meta
    qc_dada2_table
    qc_dada2_rep
    import_refseq
    import_taxa
    generate_classifier_bayes
    taxa_classification
    taxa_filtering
    qc_init_classification
    qc_filtered_classification
    krona_init_classification
    krona_filtered_classification
    create_info
    qiime_info
    cutadapt_info
    publish_info
} from "./modules/modules_qiime2_amplicons.nf"


// -----------------------------------------------------------------------------
// Workflow execution
workflow {
    // Debug
    // inputs_ch.view()
    // reads_ch.view()
    // taxa_ch.view()

    // Create the classifier for Qiime2
    import_refseq(reads_ch)
    import_taxa(taxa_ch)
    generate_classifier_bayes(import_refseq.out, import_taxa.out)

    // Data trimming if asked
    def samples_ch = params.trimming \
        ? trimming_cutadapt(inputs_ch) \
        : inputs_ch

    // Create the manifest for Qiime2
    generate_manifest(samples_ch)

    // Create a big manifest if all results in 1 file
    def manifests_ch = params.all_in_one \
        ? generate_manifest_all(generate_manifest.out.collect()) \
        : generate_manifest.out.collect()

    // Qiime2 import
    import_manifest(manifests_ch)
    qc_demux(import_manifest.out)

    // Denoising by DADA2 + QC
    denoise_dada2(import_manifest.out)
    qc_dada2_meta(denoise_dada2.out.stats_dada2)
    qc_dada2_table(denoise_dada2.out.table_dada2)
    qc_dada2_rep(denoise_dada2.out.rep_dada2)

    // Reads classification
    taxa_classification(generate_classifier_bayes.out, denoise_dada2.out.rep_dada2)

    // QC on classification without filtering    
    joined_taxa_table_ch = taxa_classification.out.join(denoise_dada2.out.table_dada2)
    qc_init_classification(joined_taxa_table_ch)
    krona_init_classification(joined_taxa_table_ch)

    // QC on classification with filtering
    taxa_filtering(joined_taxa_table_ch)
    qc_filtered_classification(taxa_filtering.out)
    krona_filtered_classification(taxa_filtering.out)

    //General information for tracking
    create_info(
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
    qiime_info(create_info.out)
    cutadapt_info(qiime_info.out)
    publish_info(cutadapt_info.out)
}