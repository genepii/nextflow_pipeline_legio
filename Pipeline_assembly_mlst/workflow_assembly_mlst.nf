#!/usr/bin/env nextflow
    
// enable dsl2
nextflow.enable.dsl=2


// -----------------------------------------------------------------------------
// MAIN WORFLOW FOR ASSEMBLY + MLST ANALYSIS
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
    MERGE_TOTAL_BASES
    MPA_FAMILY_BARPLOT
    RATIO_KRAKEN2
    MERGE_RATIO_KRAKEN2
    MLST_ELGATO
    MERGE_ELGATO
    CLEAN_REFERENCE_FASTA
    INDEX_REFERENCE_FASTA
    ALIGN_MINIMAP2
    TO_BAM_SAMTOOLS
    VARCALL_FREEBAYES
    FILTER_BCFTOOLS
    NORM_BCFTOOLS
    ASSEMBLY_SPADES
    FILTER_CONTIGS
    QC_QUAST
    MERGE_QC_QUAST
    MLST_MOMPS
    MERGE_MOMPS
    STRAIN_FASTANI
    MERGE_STRAIN_FASTANI
    MLST_CHEWBBACA
    EXTRACT_ALLELES
    MERGE_EXTRACT_ALLELES
    CHEWBBACA_GRAPETREE
    MERGE_REPORTREE_TSV
    CHEWBBACA_REPORTREE
    LP_MERGE_METADATA
    VISU_REPORTREE
    ASSEMBLY_MLST_SUMMARY_TABLE 
    ASSEMBLY_MLST_SUMMARY_HTML
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
    MINIMAP_INFO
    FREEBAYES_INFO
    BCFTOOLS_INFO
    SNPEFF_INFO
    FASTANI_INFO
    SPADES_INFO
    SAMTOOLS_INFO
    QUAST_INFO
    MOMPS_INFO
    CHEWBBACA_INFO
    GRAPETREE_INFO
    REPORTREE_INFO
    PUBLISH_INFO
} from "./modules/modules_assembly_mlst.nf"

include { 
    QC_FASTQC as QC_FASTQC_RAW 
    QC_MULTIQC as QC_MULTIQC_RAW
    LP_GRAPETREE as LP_GRAPETREE_MOMPS
    IMPACT_SNPEFF as IMPACT_SNPEFF_AMR
    MERGE_IMPACT_SNPEFF as MERGE_IMPACT_SNPEFF_AMR
    PARSE_SNPEFF_GENES as PARSE_SNPEFF_GENES_AMR
} from './modules/modules_assembly_mlst.nf'

include { 
    QC_FASTQC as QC_FASTQC_TRIM 
    QC_MULTIQC as QC_MULTIQC_TRIM
    LP_GRAPETREE as LP_GRAPETREE_ELGATO
    IMPACT_SNPEFF as IMPACT_SNPEFF_OTHER
    MERGE_IMPACT_SNPEFF as MERGE_IMPACT_SNPEFF_OTHER
    PARSE_SNPEFF_GENES as PARSE_SNPEFF_GENES_OTHER
} from './modules/modules_assembly_mlst.nf'

include { 
    QC_FASTQC as QC_FASTQC_PROC 
    QC_MULTIQC as QC_MULTIQC_PROC
} from './modules/modules_assembly_mlst.nf'


// -----------------------------------------------------------------------------
// WORKFLOW
// -----------------------------------------------------------------------------

workflow {
    // ---------------------------
    // TRACKING CONFIG
    // ---------------------------
    // NB : We could use a config.json file to avoid this lengthy call, 
    // but it was decided to stick with a standard Nextflow config.txt file 
    CREATE_INFO(
        params.suffix,
        params.input_dir,
        params.result,

        params.paired_end,
        params.adapters,
        params.decontamination,
        params.downsampling,
        params.momps,
        params.snpeff_other,

        params.bbtools_downsampled,

        params.min_quality,
        params.min_length,

        params.bbwrap_ref,
        params.bbwrap_path,
        params.bbwrap_min_id,
        params.bbwrap_max_indel,
        params.bbwrap_bwr,
        params.bbwrap_bw,
        params.bbwrap_min_hits,
        params.bbwrap_qtrim,
        params.bbwrap_trimq,
        params.bbwrap_qin,

        params.kraken2_db,
        params.format_mpa,

        params.spp_target,
        params.min_ratio_target,
        params.min_ratio_legio,
        params.min_ratio_legia,
        params.elgato_depth,

        params.minimap_ref,
        params.minimap_frag,
        params.minimap_optF,
        params.minimap_optk,
        params.minimap_optw,
        params.minimap_optA,
        params.minimap_optB,
        params.minimap_optO,
        params.minimap_optE,
        params.minimap_optr,
        params.minimap_optp,
        params.minimap_optN,
        params.minimap_optf,
        params.minimap_optn,
        params.minimap_optm,
        params.minimap_opts,
        params.minimap_optg,
        params.minimap_optheap,
        params.minimap_optsec,

        params.freeb_targets,
        params.freeb_theta,
        params.freeb_ploidy,
        params.freeb_best_n,
        params.freeb_haplo_len,
        params.freeb_max_it,
        params.freeb_max_dep,
        params.freeb_min_mapqual,
        params.freeb_min_basequal,
        params.freeb_min_var,
        params.freeb_min_dep,

        params.bcf_min_freq,
        params.bcf_qa,

        params.snpeff_amr_config,
        params.snpeff_amr_scheme,
        params.snpeff_other_config,
        params.snpeff_other_scheme,

        params.min_length_contig,

        params.fastani_genomes,
        params.fastani_min,

        params.rep_metadata,
        params.rep_partition,
        params.rep_interest,
        params.rep_zoom,
        params.rep_site_inclusion,
        params.rep_min_allele,
        params.rep_max_allele,
        params.rep_loci_called,
        params.rep_col_metadata,

        // IMPORTANT: serialize complex objects
        groovy.json.JsonOutput.toJson(params.lb_set),
        groovy.json.JsonOutput.toJson(params.lp_set),
        groovy.json.JsonOutput.toJson(params.alleles_set)
    )

    FASTQC_INFO(CREATE_INFO.out)
    MULTIQC_INFO(FASTQC_INFO.out)
    tmp_out = FASTP_INFO(MULTIQC_INFO.out)
    if (params.downsampling || params.decontamination) {
        BBTOOLS_INFO(tmp_out)
        tmp_out = SEQKIT_INFO(BBTOOLS_INFO.out)
    }
    KRAKEN2_INFO(tmp_out)
    PYTHON_INFO(KRAKEN2_INFO.out)
    tmp_out = KRONA_INFO(PYTHON_INFO.out)
    ELGATO_INFO(tmp_out)
    MINIMAP_INFO(ELGATO_INFO.out)
    FREEBAYES_INFO(MINIMAP_INFO.out)
    BCFTOOLS_INFO(FREEBAYES_INFO.out)
    SNPEFF_INFO(BCFTOOLS_INFO.out)
    FASTANI_INFO(SNPEFF_INFO.out)
    SPADES_INFO(FASTANI_INFO.out)
    SAMTOOLS_INFO(SPADES_INFO.out)
    QUAST_INFO(SAMTOOLS_INFO.out)
    MOMPS_INFO(QUAST_INFO.out)
    CHEWBBACA_INFO(MOMPS_INFO.out)
    GRAPETREE_INFO(CHEWBBACA_INFO.out)
    REPORTREE_INFO(GRAPETREE_INFO.out)

    PUBLISH_INFO(REPORTREE_INFO.out)


    // ---------------------------
    // ANALYSIS - START
    // ---------------------------

    // only Paired End data
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
    ASSIGN_KRAKEN2(samples_ch)
    MPA_MODIF(ASSIGN_KRAKEN2.out)
    MPA_TO_KRONA(MPA_MODIF.out)

    COUNT_FASTQ_READS(samples_ch)
    joined_mpamodif_total_ch = MPA_MODIF.out.join(COUNT_FASTQ_READS.out.total)
    MPA_FAMILY_BARPLOT(joined_mpamodif_total_ch)

    total_bases_ch = COUNT_FASTQ_READS.out.base
        .collect()
    MERGE_TOTAL_BASES(total_bases_ch)

    
    // ---------------------------
    // ELGATO - MLST PROFILE (only if spp_target)
    // ---------------------------
    joined_mpa_total_ch = ASSIGN_KRAKEN2.out.join(COUNT_FASTQ_READS.out.total)
    ratio_ch = RATIO_KRAKEN2(joined_mpa_total_ch)

    merge_ratio_ch = ratio_ch
        .map { sample_id, ratio_file -> ratio_file }
        .collect()
    MERGE_RATIO_KRAKEN2(merge_ratio_ch)

    // Get only samples identified as params.spp_target 
    filtered_samples_ch =
        ratio_ch
            .map { sample_id, ratio_file ->

                // Read file in a Nextflow-safe staged context (process-safe execution)
                def content = ratio_file.text

                if (content.contains(params.spp_target)) {
                    tuple(sample_id, true)
                } else {
                    tuple(sample_id, false)
                }
            }
            .filter { sample_id, keep -> keep }
            .map { sample_id, keep -> sample_id }

    elgato_input_ch =
        samples_ch
            .map { sample_id, r1, r2 -> tuple(sample_id, r1, r2) }
            .join(filtered_samples_ch)
            .map { sample_id, r1, r2 ->
                tuple(sample_id, r1, r2)
            }

    MLST_ELGATO(elgato_input_ch)

    mlst_elgato_ch = MLST_ELGATO.out.mlst
        .map { sample_id, file -> file }
        .collect()
    fastfinder_elgato_ch = MLST_ELGATO.out.fastfinder
        .map { sample_id, file -> file }
        .collect()
    MERGE_ELGATO(mlst_elgato_ch, fastfinder_elgato_ch)

    // GrapeTree only if > 2 samples == Lp
    filtered_mlst_ch = MERGE_ELGATO.out.mlst.filter { tsv ->
        tsv.text.readLines().size() > 3
    }
    LP_GRAPETREE_ELGATO("ElGato", filtered_mlst_ch)

    
    // ---------------------------
    // FREEBAYES - AMR PROFILE (only if spp_target)
    // ---------------------------
    ref_file = CLEAN_REFERENCE_FASTA(params.minimap_ref)
    INDEX_REFERENCE_FASTA(ref_file)

    ALIGN_MINIMAP2(elgato_input_ch, ref_file)
    TO_BAM_SAMTOOLS(ALIGN_MINIMAP2.out)
    
    VARCALL_FREEBAYES(TO_BAM_SAMTOOLS.out, ref_file)
    FILTER_BCFTOOLS(VARCALL_FREEBAYES.out)

    NORM_BCFTOOLS(FILTER_BCFTOOLS.out, ref_file)
    
    IMPACT_SNPEFF_AMR(NORM_BCFTOOLS.out, "AMR")
    PARSE_SNPEFF_GENES_AMR(IMPACT_SNPEFF_AMR.out.vcf)

    snpeff_genes_ch = IMPACT_SNPEFF_AMR.out.genes
        .map { sample_id, type, file -> file }
        .collect()
    MERGE_IMPACT_SNPEFF_AMR(snpeff_genes_ch, "AMR")

    if (params.snpeff_other) {
        IMPACT_SNPEFF_OTHER(NORM_BCFTOOLS.out, "Other")
        PARSE_SNPEFF_GENES_OTHER(IMPACT_SNPEFF_OTHER.out.vcf)

        snpeff_genes_ch = IMPACT_SNPEFF_OTHER.out.genes
            .map { sample_id, type, file -> file }
            .collect()
        MERGE_IMPACT_SNPEFF_OTHER(snpeff_genes_ch, "Other")
    }

    // ---------------------------
    // ASSEMBLY + QC + FILTER
    // ---------------------------
    ASSEMBLY_SPADES(samples_ch)
    FILTER_CONTIGS(ASSEMBLY_SPADES.out.contigs)
    
    all_contigs_ch = ASSEMBLY_SPADES.out.contigs
        .mix(FILTER_CONTIGS.out)
        .map { sample_id, file -> file }
        .collect()

    QC_QUAST(all_contigs_ch)

    MERGE_QC_QUAST(QC_QUAST.out.extract)
    

    // ---------------------------
    // MOMPS - MLST PROFILE (optionnal)
    // ---------------------------
    if (params.momps) {
        joined_reads_contigs_ch = samples_ch.join(FILTER_CONTIGS.out.filtered)
        MLST_MOMPS(joined_reads_contigs_ch)
        fastfinder_momps_ch = MLST_MOMPS.out.fastfinder
                        .map { sample_id, file -> file }
                        .collect()
        mlst_momps_ch = MLST_MOMPS.out.mlst
                        .map { sample_id, file -> file }
                        .collect()
        MERGE_MOMPS(fastfinder_momps_ch, mlst_momps_ch)

        // GrapeTree only if > 2 samples == Lp 
        filtered_momps_mlst = MERGE_MOMPS.out.mlst.filter { tsv ->
            tsv.text.readLines().size() > 3
        }
        LP_GRAPETREE_MOMPS("mompS", filtered_momps_mlst)
    }


    // ---------------------------
    // STRAIN IDENTIFICATION
    // ---------------------------
    STRAIN_FASTANI(FILTER_CONTIGS.out)

    fastani_files_ch = STRAIN_FASTANI.out
        .map { sample_id, tsv -> tsv }
        .collect()
    MERGE_STRAIN_FASTANI(fastani_files_ch)

    // Get tuple [sample_id, tsv, strain, fastANI_value], 
    // if fastANI_value >= params.fastani_min
    strain_ch = STRAIN_FASTANI.out
        .map { sample_id, tsv ->
            def cols = tsv.text.readLines()[0].split('\t')

            def strain = cols[3].tokenize('/').last()
            def fastani = cols[4] as Double

            tuple(
                sample_id,
                tsv,
                strain,
                fastani
            )
        }
        .filter { sample_id, tsv, strain, fastani ->
            fastani >= params.fastani_min
        }

    // Get strain Lb or Lp + Join on sample_id with contig filt 
    strain_contigs_ch = strain_ch
        .map { sample_id, tsv, strain, fastani ->
            def species =
                strain?.contains("L_longbeachae") ? "Lb" :
                strain?.contains("L_pneumophila") ? "Lp" :
                null

            tuple(sample_id, species)
        }
        .filter { sample_id, species ->
            species != null
        }
        .join(FILTER_CONTIGS.out)
        .map { sample_id, species, contig ->
            tuple(species, contig)
        }
        .groupTuple()
        .map { strain, contigs_files ->
            // force list consistency
            tuple(strain, contigs_files as List)
        }
        

    // ---------------------------
    // CHEWBBACA - cgMLST PROFILE
    // ---------------------------
    // Chewbbaca inputs depending on strain
    chew_input_ch = strain_contigs_ch
        .flatMap { strain, contigs_files ->

            def sets = (strain == "Lb") ? (params.lb_set ?: []) : (params.lp_set ?: [])

            sets.collect { s ->
                tuple(
                    strain,
                    contigs_files,
                    file(s.genes),
                    file(s.ptf),
                    s.nb,
                    "${strain}_${s.nb}"   // unique id
                )
            }
        }

    MLST_CHEWBBACA(chew_input_ch)
    EXTRACT_ALLELES(MLST_CHEWBBACA.out.mlst)

    allele_files_ch = EXTRACT_ALLELES.out
        .map { strain, nb, file -> file }
        .collect()
    MERGE_EXTRACT_ALLELES(allele_files_ch)

    // Build a lookup table: "strain:nb" -> flags
    def treeCfg = [:]
    (params.lp_set ?: []).each { s ->
        treeCfg["Lp:${s.nb}"] = [
            grapetree: s.grapetree,
            reportree: s.reportree
        ]
    }
    (params.lb_set ?: []).each { s ->
        treeCfg["Lb:${s.nb}"] = [
            grapetree: s.grapetree,
            reportree: s.reportree
        ]
    }

    // Merge both sources
    all_mlst_ch = EXTRACT_ALLELES.out.mix(MLST_CHEWBBACA.out.mlst)


    // ---------------------------
    // XXX TREE - cgMLST TREES
    // ---------------------------
    // Two channels, 1 for GrapeTree / 1 for ReporTree
    annotated_ch = all_mlst_ch.map { strain, nb, file ->
        def cfg = treeCfg["${strain}:${nb}"] ?: [grapetree:false, reportree:false]
        tuple(strain, nb, file, cfg.grapetree, cfg.reportree)
    }

    grapetree_ch = annotated_ch
        .filter { s, nb, f, g, r -> g }
        .map { s, nb, f, g, r ->
            tuple(s, nb, f)
        }
    reportree_ch = annotated_ch
        .filter { s, nb, f, g, r -> r }
        .map { s, nb, f, g, r ->
            tuple(s, nb, f)
        }

    // Keep only files containing >= 2 samples
    filtered_grapetree_ch = grapetree_ch
        .filter { s, nb, f ->
            def lines = f.text.readLines().findAll { it.trim() }

            def n_lines = lines.size()
            def n_cols = n_lines > 0 ? lines[0].split('\t', -1).size() : 0

            n_lines > 3 && n_cols >= 2
        }
    CHEWBBACA_GRAPETREE(filtered_grapetree_ch)

    // ReporTree not on .filt.tsv
    filtered_chewbbaca_ch = reportree_ch
        .filter { s, nb, f ->
            def lines = f.text.readLines().findAll { it.trim() }
            def n_lines = lines.size()
            def n_cols = n_lines > 0 ? lines[0].split('\t', -1).size() : 0

            !f.name.endsWith(".filt.tsv") &&
            n_cols >= 2
        }
    MERGE_REPORTREE_TSV(filtered_chewbbaca_ch)

    // ReporTree only if > 2 samples
    filtered_mergedchewbbaca_ch = MERGE_REPORTREE_TSV.out.mlst
        .filter { s, nb, f ->
            def lines = f.text.readLines().findAll { it.trim() }
            def n_lines = lines.size()
            n_lines > 3
        }
    CHEWBBACA_REPORTREE(filtered_mergedchewbbaca_ch)

    LP_MERGE_METADATA(LP_GRAPETREE_ELGATO.out.metadata)

    // Add cgMLST reference schema
    def lp_lit_map = params.lp_set.collectEntries { cfg ->
        [(cfg.nb): cfg.lit]
    }
    filtered_reportree_ch = CHEWBBACA_REPORTREE.out.mlst
        .map { strain, nb, allele_tsv ->
            def lit = lp_lit_map[nb]
            return tuple(strain, nb, allele_tsv, lit)
        }

    VISU_REPORTREE(filtered_reportree_ch, LP_MERGE_METADATA.out)
        

    // ---------------------------
    // SUMMARY - .tsv and .html
    // ---------------------------
    summary_ch = MERGE_RATIO_KRAKEN2.out
        .mix(MERGE_STRAIN_FASTANI.out)
        .mix(MERGE_ELGATO.out.mlst)
        .mix(MERGE_EXTRACT_ALLELES.out)
        .mix(MERGE_QC_QUAST.out)
        .mix(MERGE_IMPACT_SNPEFF_AMR.out)
        .mix(MERGE_TOTAL_BASES.out)

    if (params.momps) {
        summary_ch = summary_ch.mix(MERGE_MOMPS.out.mlst)
    }
    if (params.snpeff_other) {
        summary_ch = summary_ch.mix(MERGE_IMPACT_SNPEFF_OTHER.out)
    }
    ASSEMBLY_MLST_SUMMARY_TABLE(summary_ch.collect())

    ASSEMBLY_MLST_SUMMARY_HTML(ASSEMBLY_MLST_SUMMARY_TABLE.out, PUBLISH_INFO.out)
}
