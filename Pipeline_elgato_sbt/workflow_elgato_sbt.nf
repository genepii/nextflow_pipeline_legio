#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2


// -----------------------------------------------------------------------------
// MAIN WORFLOW FOR EL GATO NESTED SBT ANALYSIS
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// INPUT CHANNEL
// -----------------------------------------------------------------------------

// collect fastq files in tuple [sample_id, R1, R2] or [sample_id, R1]
if (params.paired_end) {
    inputs_ch = Channel
        .fromFilePairs("${params.input_dir}/*_{R1,R2}.fastq*")
        .map { id, reads ->

            def base = id

            def sample_id = base.contains('_') ?
                base.split('_')[0] :
                base.replaceFirst(/_R1.*/, '')

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

            def base = r1.baseName

            def sample_id = base.contains('_') ?
                base.split('_')[0] :
                base.replaceFirst(/_R1.*/, '')

            assert sample_id
            assert r1

            tuple(sample_id, r1, null)
        }
}


// -----------------------------------------------------------------------------
// INCLUDE MODULES
// -----------------------------------------------------------------------------

include {
    TRIM_FASTP
    DECONTA_BBWRAP
    DOWNSAMPLE_BBTOOLS
    QC_SEQKIT
    ASSIGN_KRAKEN2
    MPA_MODIF
    MPA_TO_KRONA
    COUNT_FASTQ_READS
    MPA_FAMILY_BARPLOT
    MLST_ELGATO
    MERGE_ELGATO
    CREATE_INFO
    FASTQC_INFO
    MULTIQC_INFO
    FASTP_INFO
    BBTOOLS_INFO
    SEQKIT_INFO
    KRAKEN2_INFO
    PYTHON_INFO
    KRONA_INFO
    ELGATO_INFO
    PUBLISH_INFO
} from './modules/modules_elgato_sbt.nf'

include { 
    QC_FASTQC as QC_FASTQC_RAW 
    QC_MULTIQC as QC_MULTIQC_RAW
} from './modules/modules_elgato_sbt.nf'

include { 
    QC_FASTQC as QC_FASTQC_TRIM 
    QC_MULTIQC as QC_MULTIQC_TRIM
} from './modules/modules_elgato_sbt.nf'

include { 
    QC_FASTQC as QC_FASTQC_PROC 
    QC_MULTIQC as QC_MULTIQC_PROC
} from './modules/modules_elgato_sbt.nf'


// -----------------------------------------------------------------------------
// WORKFLOW
// -----------------------------------------------------------------------------

workflow {

    // ---------------------------
    // ONLY PAIRED-END DATA
    // ---------------------------
    if (!params.paired_end.toBoolean()) {
        error "Only paired-end (PE) data are supported. Single-end mode is not allowed."
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
    // DECONTAMINATION (optionnal)
    // ---------------------------
    read_type_proc = "1_Processed"

    if (params.decontamination) {
        samples_ch = DECONTA_BBWRAP(samples_ch)
    }

    // ---------------------------
    // DOWNSAMPLING (optionnal)
    // ---------------------------
    if (params.downsampling) {
        samples_ch = DOWNSAMPLE_BBTOOLS(samples_ch)
        QC_SEQKIT(DOWNSAMPLE_BBTOOLS.out)
    }

    // ---------------------------
    // preprocessing QC PLOTS (optionnal)
    // ---------------------------
    if (params.downsampling || params.decontamination) {
        qc_processed = QC_FASTQC_PROC(read_type_proc, samples_ch)

        proc_fastqc_zips = qc_processed.zip_files
            .map { sample_id, zip -> zip }
            .flatten()
            .collect()

        QC_MULTIQC_PROC(read_type_proc, proc_fastqc_zips)
    }

    // ---------------------------
    // TAXONOMIC ASSIGNATION + QC (optional)
    // ---------------------------
    if (params.kraken2_assign) {
        ASSIGN_KRAKEN2(samples_ch)
        MPA_MODIF(ASSIGN_KRAKEN2.out)
        MPA_TO_KRONA(MPA_MODIF.out)

        COUNT_FASTQ_READS(samples_ch)
        joined_mpa_total_ch = MPA_MODIF.out.join(COUNT_FASTQ_READS.out)
        MPA_FAMILY_BARPLOT(joined_mpa_total_ch)
    }

    // ---------------------------
    // MLST PROFILE
    // ---------------------------
    MLST_ELGATO(samples_ch)
    mlst_files_ch = MLST_ELGATO.out.mlst
                    .map { sample_id, file -> file }
                    .collect()
    MERGE_ELGATO(mlst_files_ch)


    // ---------------------------
    // TRACKING CONFIG
    // ---------------------------
    CREATE_INFO(
        params.suffix,
        params.input_dir,
        params.result,

        params.paired_end,
        params.adapters,
        params.decontamination,
        params.downsampling,
        params.kraken2_assign,

        params.min_quality,
        params.min_length,

        params.bbwrap_ref,
        params.bbwrap_min_id,
        params.bbwrap_max_indel,
        params.bbwrap_bwr,
        params.bbwrap_bw,
        params.bbwrap_min_hits,
        params.bbwrap_qtrim,
        params.bbwrap_trimq,
        params.bbwrap_qin,
        params.bbwrap_path,

        params.bbtools_downsampled,

        params.kraken2_db,

        params.elgato_depth
    )

    FASTQC_INFO(CREATE_INFO.out)
    MULTIQC_INFO(FASTQC_INFO.out)
    tmp_out = FASTP_INFO(MULTIQC_INFO.out)
    if (params.downsampling || params.decontamination) {
        BBTOOLS_INFO(tmp_out)
        tmp_out = SEQKIT_INFO(BBTOOLS_INFO.out)
    }
    if (params.kraken2_assign) {
        KRAKEN2_INFO(tmp_out)
        PYTHON_INFO(KRAKEN2_INFO.out)
        tmp_out = KRONA_INFO(PYTHON_INFO.out)
    }
    ELGATO_INFO(tmp_out)

    PUBLISH_INFO(ELGATO_INFO.out)
}