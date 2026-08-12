#!/usr/bin/env nextflow
    
// enable dsl2
nextflow.enable.dsl=2


// -----------------------------------------------------------------------------
// MAIN WORFLOW FOR COMPARISON OF STRAINS
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// INCLUDE MODULES
// -----------------------------------------------------------------------------

include {
    FIND_FASTA
    FIND_FASTQ
    BUILD_METADATA
    CHECK_METADATA
    TRIM_FASTP
    PREPARE_ALLELES_TSV
    CHECK_ALLELES_TSV
    FILTER_GUBBINS
    VISU_GRAPETREE
    VISU_REPORTREE
} from "./modules/modules_strain_comparison.nf"

include {
    SNP_SNIPPY as INTRA_SNIPPY
    CHECK_SNIPPY_CORE as CHECK_INTRA_CORE
    SNP_SNIPPY_CORE as INTRA_CORE
    SNP_DIST_CORE as INTRA_DIST_CORE
    SNIPPY_CORE_LOG as INTRA_LOG
    PHYLOTREE_IQTREE as INIT_IQTREE
    TREE_LOG as INIT_LOG
    PLOT_IQTREE as INIT_PLOT
    VISU_LOG as GRAPETREE_LOG
} from "./modules/modules_strain_comparison.nf"

include {
    SNP_SNIPPY as ST_SNIPPY
    SNP_SNIPPY_CORE as ST_CORE
    SNP_DIST_CORE as ST_DIST_CORE
    SNIPPY_CORE_LOG as ST_LOG
    PHYLOTREE_IQTREE as FILT_IQTREE
    TREE_LOG as FILT_LOG
    PLOT_IQTREE as FILT_PLOT
    VISU_LOG as REPORTREE_LOG
} from "./modules/modules_strain_comparison.nf"

include {
    VISU_LOG as CHECK_LOG
    TREE_LOG as GUBBINS_LOG
} from "./modules/modules_strain_comparison.nf"


// -----------------------------------------------------------------------------
// WORKFLOW
// -----------------------------------------------------------------------------

workflow {
    // ---------------------------
    // METADATA FILE SCREENING
    // ---------------------------
    FIND_FASTA(params.metadata)
    FIND_FASTQ(params.metadata)

    BUILD_METADATA(FIND_FASTA.out.newpath, FIND_FASTQ.out.newpath)
    CHECK_METADATA(BUILD_METADATA.out)


    // ---------------------------
    // TRIMMING
    // ---------------------------
    fastq_for_trim_ch = CHECK_METADATA.out.fastq
        .splitCsv(
            header: true,
            sep: "\t"
        )
        .map { row ->
            tuple(
                row.ID,
                row.Comparison,
                file(row.READS_1),
                file(row.READS_2)
            )
        }
    TRIM_FASTP(fastq_for_trim_ch)


    // ---------------------------
    // CORE + SNP - 1st fasta vs fastq
    // ---------------------------
    // Get the 1st fasta found for each comparison value
    ref_by_comparison_ch = CHECK_METADATA.out.fasta
        .splitCsv(
            header: true,
            sep: "\t"
        )
        .map { row ->
            tuple(
                row.Comparison,
                file(row.FASTA)
            )
        }
        .groupTuple()
        .map { comparison, fasta_list ->
            tuple(
                fasta_list[0],
                comparison
            )
        }

    // SNP research for each sample by comparison
    snippy_intra_ch = TRIM_FASTP.out
        .combine(ref_by_comparison_ch, by: 1)
        .map { comparison, sample_id, r1, r2, reference ->
            tuple(
                sample_id,
                comparison,
                reference,
                r1,
                r2
            )
        }
    INTRA_SNIPPY("Intragroup", snippy_intra_ch)

    // SNP core research for each comparison
    all_intra_ch = INTRA_SNIPPY.out
        .map { sample_id, comparison, type, reference, snippy_dir ->
            tuple(
                comparison,
                type,
                reference,
                snippy_dir
            )
        }
        .groupTuple()
        .map { comparison, type_list, reference_list, snippy_dir_list ->
            tuple(
                comparison,
                type_list[0],
                reference_list[0],
                snippy_dir_list
            )
        }

    // Create a LOG if no Variant detected or only 1 sample
    CHECK_INTRA_CORE(all_intra_ch)

    INTRA_CORE(all_intra_ch)

    intra_logs_ch = CHECK_INTRA_CORE.out.log
        .mix(INTRA_CORE.out.log)
        .collect()
    INTRA_LOG("Intragroup", intra_logs_ch)

    INTRA_DIST_CORE(INTRA_CORE.out.core)


    // ---------------------------
    // SNP - ref ST vs fastq
    // ---------------------------
    // Create reference channel from config
    snippy_ref_ch = Channel
        .fromList(params.snippy_ref)
        .map { snippy_info ->
            tuple(
                snippy_info.st.toString(),
                file(snippy_info.ref)
            )
        }

    // Metadata containing sample ID and ST
    metadata_st_ch = CHECK_METADATA.out.fastq
        .splitCsv(
            header: true,
            sep: "\t"
        )
        .map { row ->
            tuple(
                row.ID,
                row.ST.toString()
            )
        }

    // Add ST information to trimmed reads
    reads_with_st_ch = TRIM_FASTP.out
        .join(metadata_st_ch)
        .map { sample_id, comparison, r1, r2, st ->
            tuple(
                sample_id,
                st,
                r1,
                r2
            )
        }

    // Match reads ST with reference ST
    snippy_st_ch = reads_with_st_ch
        .combine(snippy_ref_ch)
        .filter { sample_id, st, r1, r2, ref_st, reference ->
            st == ref_st
        }
        .map { sample_id, st, r1, r2, ref_st, reference ->
            tuple(
                sample_id,
                st,
                reference,
                r1,
                r2
            )
        }

    ST_SNIPPY("ST", snippy_st_ch)


    // ---------------------------
    // CORE SNP - ref ST vs fastq + previous
    // ---------------------------
    // Get basedir for previous Snippy results for each ST
    previous_snp_ch = Channel
        .fromList(params.snippy_ref)
        .map { snippy_info ->
            tuple(
                snippy_info.st.toString(),
                file(snippy_info.previous)
                    .listFiles()
                    .findAll { it.isDirectory() }
            )
        }

    // SNP core research for each ST reference
    snp_st_ch = ST_SNIPPY.out
        .map { sample_id, st, type, reference, snippy_dir ->
            tuple(
                st.toString(),
                type,
                reference,
                snippy_dir
            )
        }
        .groupTuple()
        .map { st, type_list, reference_list, snippy_dir_list ->
            tuple(
                st,
                type_list[0],
                reference_list[0],
                snippy_dir_list
            )
        }

    // Merge new results with previous results using ST/comparison, remove duplicates
    snp_with_previous_ch = snp_st_ch
        .join(previous_snp_ch)
        .map { st, type, reference, snippy_dirs, previous_dirs ->
            def snippy_ids = snippy_dirs.collect { it.baseName.toString() }

            def filtered_previous_dirs = previous_dirs.findAll { prev_dir ->
                !snippy_ids.contains(prev_dir.baseName.toString())
            }

            tuple(
                st,
                type,
                reference,
                snippy_dirs + filtered_previous_dirs
            )
        }

    ST_CORE(snp_with_previous_ch)

    ST_LOG("ST", ST_CORE.out.log.collect())

    ST_DIST_CORE(ST_CORE.out.core)


    // ---------------------------
    // PHYLOGENY - ref ST vs fastq + previous
    // ---------------------------
    INIT_IQTREE("Init_ST", ST_CORE.out.alignment)

    init_log_ch = INIT_IQTREE.out.log
        .map { type, log ->
            tuple(type, log)
        }
        .groupTuple()
    INIT_LOG(init_log_ch)

    INIT_PLOT(INIT_IQTREE.out.tree, BUILD_METADATA.out)

    FILTER_GUBBINS("ST", ST_CORE.out.alignment)

    gubbins_log_ch = FILTER_GUBBINS.out.log
        .map { type, log ->
            tuple(type, log)
        }
        .groupTuple()
    GUBBINS_LOG(gubbins_log_ch)

    FILT_IQTREE("Filt_ST", FILTER_GUBBINS.out.filt)

    filt_log_ch = FILT_IQTREE.out.log
        .map { type, log ->
            tuple(type, log)
        }
        .groupTuple()
    FILT_LOG(filt_log_ch)

    FILT_PLOT(FILT_IQTREE.out.tree, BUILD_METADATA.out)


    // ---------------------------
    // VISUALISATION - GrapeTree
    // ---------------------------
    // Extract unique comparison values
    comparisons_ch = BUILD_METADATA.out
        .splitCsv(
            header: true,
            sep: '\t'
        )
        .map { row ->
            row.Comparison
        }
        .unique()
    
    // Create cgMLST schemes enabled for GrapeTree
    grapetree_schema_ch = Channel
        .from(params.chewbbaca_set)
        .filter { scheme ->
            scheme.grapetree
        }
        .map { scheme ->
            tuple(
                file(scheme.genes),
                scheme.nb,
                file(scheme.previous)
            )
        }
    
    // Create input tuples for TSV preparation
    grapetree_prepare_ch = comparisons_ch
        .combine(BUILD_METADATA.out)
        .combine(grapetree_schema_ch)
        .map { comparison, metadata_file, genes, nb, previous ->
            tuple(
                metadata_file,
                comparison,
                genes,
                nb,
                previous
            )
        }
        
    PREPARE_ALLELES_TSV(grapetree_prepare_ch)

    VISU_GRAPETREE(PREPARE_ALLELES_TSV.out)

    grapetree_log_ch = VISU_GRAPETREE.out.log
        .map { type, log ->
            tuple(type, log)
        }
        .groupTuple()
    GRAPETREE_LOG(grapetree_log_ch)


    // ---------------------------
    // VISUALISATION - ReporTree
    // ---------------------------    
    // Create cgMLST schemes enabled for ReporTree
    reportree_schema_ch = Channel
        .from(params.chewbbaca_set)
        .filter { scheme ->
            scheme.reportree
        }
        .map { scheme ->
            tuple(
                file(scheme.genes),
                scheme.nb,
                file(scheme.previous),
                file(scheme.nomenclature)
            )
        }
    
    // Create input tuples for TSV preparation
    reportree_prepare_ch = comparisons_ch
        .combine(BUILD_METADATA.out)
        .combine(reportree_schema_ch)
        .map { comparison, metadata_file, genes, nb, previous, nomenclature ->
            tuple(
                metadata_file,
                comparison,
                genes,
                nb,
                previous
            )
        }

    CHECK_ALLELES_TSV(reportree_prepare_ch)

    check_log_ch = CHECK_ALLELES_TSV.out.log
        .map { type, log ->
            tuple(type, log)
        }
        .groupTuple()

    CHECK_LOG(check_log_ch)
    
    // Get tuples for ReporTree, join by chewbbaca_set nb
    reportree_visu_ch = CHECK_ALLELES_TSV.out.alleles
        .filter { alleles, nb, comparison, genes, metadata, error ->
            error.text.trim() == "0"
        }
        .combine(reportree_schema_ch, by: 1)
        .map { nb, alleles, comparison, genes, metadata, error, genes_2, previous, nomenclature ->
            tuple(
                alleles,
                nb,
                comparison,
                genes,
                metadata,
                nomenclature
            )
        }

    VISU_REPORTREE(reportree_visu_ch)

    reportree_log_ch = VISU_REPORTREE.out.log
        .map { type, log ->
            tuple(type, log)
        }
        .groupTuple()
    REPORTREE_LOG(reportree_log_ch)
}