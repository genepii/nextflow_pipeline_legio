#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Processes for workflow_assembly_mlst.nf


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
* Reference-based decontamination using BBMap
* Input   : paired-end FASTQ files (R1 and R2) per sample
* Output  : paired FASTQ files containing only unmapped reads (outu1/outu2)
* Purpose : remove reads mapping to a reference database (host, PhiX, contaminants)
*           and retain non-aligned reads for downstream microbial profiling
* Note    : pairing is preserved natively by BBMap (no manual reconstruction required)
*/
process DECONTA_BBWRAP {
    label 'bbtools'
    publishDir "${params.result}/dev/0-2_Decontamination", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id),
            path("${sample_id}_R1_deconta.fastq.gz"),
            path("${sample_id}_R2_deconta.fastq.gz")

    script:
    def ref_opts = params.bbwrap_ref ? 
        "ref=${params.bbwrap_ref}" : 
        "path=${params.bbwrap_path} append=t"

    """
    bbwrap.sh \
        -Xmx${task.memory.toGiga()}g \
        mapper=bbmap \
        ${ref_opts} \
        in1="${r1}" \
        in2="${r2}" \
        outu1="${sample_id}_R1_deconta.fastq" \
        outu2="${sample_id}_R2_deconta.fastq" \
        minid=${params.bbwrap_min_id} \
        maxindel=${params.bbwrap_max_indel} \
        bwr=${params.bbwrap_bwr} \
        bw=${params.bbwrap_bw} \
        minhits=${params.bbwrap_min_hits} \
        qtrim=${params.bbwrap_qtrim} \
        trimq=${params.bbwrap_trimq} \
        qin=${params.bbwrap_qin} \
        threads=${task.cpus} \
        fast=t \
        quickmatch=t \
        ow=t

    gzip "${sample_id}_R1_deconta.fastq" 
    gzip "${sample_id}_R2_deconta.fastq"
    """
}

// -----------------------------------------------------------------------------
/*
* Reads downsampling
* Input   : paired-end FASTQ files (R1, R2)
* Output  : downsampled Paired FASTQ files
* Purpose : reduce dataset size for testing or resource optimization
*/
process DOWNSAMPLE_BBTOOLS {
    label 'bbtools'
    publishDir "${params.result}/dev/0-3_Downsampled", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id),
            path("${sample_id}_R1_downsampled.fastq.gz"),
            path("${sample_id}_R2_downsampled.fastq.gz")

    script:
    """
    reformat.sh \
        in1="${r1}" in2="${r2}" \
        out1="${sample_id}_R1_downsampled.fastq.gz" \
        out2="${sample_id}_R2_downsampled.fastq.gz" \
        samplerate=${params.bbtools_downsampled} \
        -Xmx${task.memory.toGiga()}g
    """
}

/*
* QC of the reads downsampling
* Input   : downsampled paired FASTQ files (R1, R2)
* Output  : QC files with global stats + Paired Check .txt (empty if OK)
* Purpose : provide stats about downsampled data, and check if data not damaged
*/
process QC_SEQKIT {
    label 'seqkit'
    publishDir "${params.result}/dev/0-3_Downsampled/QC", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        path("${sample_id}_stats.txt")
        path("${sample_id}_paired_check.txt")

    script:
    """
    #!/usr/bin/env bash

    seqkit stats -Ta \
        "${r1}" \
        "${r2}" \
        > "${sample_id}_stats.txt"

    touch "${sample_id}_paired_check.txt"

    paste <(zcat "${r1}" | sed -n '1~4p') \
        <(zcat "${r2}" | sed -n '1~4p') \
    | awk '{if(\$1!=\$2) {print "Mismatch:", \$1, \$2; exit 1}}' \
    > "${sample_id}_paired_check.txt"
    """
}

// -----------------------------------------------------------------------------
/*
* Taxonomic classification using Kraken2
* Input   : paired-end FASTQ files (R1, R2)
* Output  : Kraken2 standard report file or MPA-style taxonomic profile
* Purpose : rapid k-mer based taxonomic assignment for quality control
*/
process ASSIGN_KRAKEN2 {
    label 'kraken2'
    publishDir "${params.result}/1_Kraken2", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}.${params.format_mpa ? 'mpa' : 'report'}")

    script:
    def output_format = params.format_mpa ? 
        "--use-mpa-style --report ${sample_id}.mpa" : 
        "--report ${sample_id}.report"

    """
    kraken2 \
        --db "${params.kraken2_db}" \
        --paired \
        --classified-out "${sample_id}#.fastq" \
        --output "${sample_id}.out" \
        ${output_format} \
        "${r1}" "${r2}" \
        --threads ${task.cpus}
    """
}

/*
* Format Kraken2 output for visualisation
* Input   : MPA-style taxonomic profile
* Output  : results formatted (.tsv), with no double count
* Purpose : explore full taxonomic composition from flattened taxonomy
*/
process MPA_MODIF {
    label 'python'
    publishDir "${params.result}/dev/1_Kraken2", mode: 'copy'

    input:
        tuple val(sample_id), path(mpa)

    output:
        tuple val(sample_id), path("${sample_id}_mpaModif.tsv")
    
    script:
    """
    assembly_mlst_mpa_modified.py ${mpa} ${sample_id}_mpaModif.tsv
    """
}

/*
* Krona visualisation from Kraken2 MPA-style report formatted
* Input   : Kraken2 results formatted (.tsv)
* Output  : interactive Krona HTML visualisation
* Purpose : explore full taxonomic composition from flattened taxonomy
*/
process MPA_TO_KRONA {
    label 'krona'
    publishDir "${params.result}/1_Kraken2", mode: 'copy'

    input:
        tuple val(sample_id), path(krona_input)

    output:
        tuple val(sample_id), path("${sample_id}_allKrona.html")

    script:
    """
    ktImportText ${krona_input} -o ${sample_id}_allKrona.html
    """
}

/*
* Total reads count recovery from FASTQ
* Input   : FASTQ files (R1, R2)
* Output  : total number of reads 
* Purpose : keep the information for plotting later
*/
process COUNT_FASTQ_READS {
    label 'seqkit'
    publishDir "${params.result}/dev/1_Kraken2", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_totalreads.txt"), emit : total
        path("${sample_id}_totalbases.txt"), emit : base

    script:
    """
    seqkit stats -T ${r1} > stats.tsv

    awk -F '\\t' 'NR==2 {gsub(/,/, "", \$4); print \$4}' stats.tsv \
        > ${sample_id}_totalreads.txt

    awk -F '\\t' 'NR==2 {gsub(/,/, "", \$5); print \$5}' stats.tsv \
    > ${sample_id}_totalbases.txt
    """
}

/*
* Merge total bases counts across samples
* Input   : list of *_totalbases.txt files (one value per sample)
* Output  : unified table with Sample_ID and total bases
* Purpose : aggregate sequencing depth information for downstream QC/plots
*/
process MERGE_TOTAL_BASES {
    label 'bash'
    publishDir "${params.result}/dev/Summary", mode: 'copy'

    input:
        path(base_files)

    output:
        path("Total_bases_summary.tsv")

    script:
    """
    echo -e "Sample_ID\tTotal_bases" > Total_bases_summary.tsv

    for f in ${base_files}; do
        sample=\$(basename \$f _totalbases.txt)
        value=\$(cat \$f)

        echo -e "\${sample}\t\${value}" >> Total_bases_summary.tsv
    done
    """
}

/*
* Family-level abundance barplot from Kraken2 MPA report
* Input   : MPA-style taxonomic profile
* Output  : TSV table + barplot (top 15 families)
* Purpose : clean visualisation of dominant bacterial families
*/
process MPA_FAMILY_BARPLOT {
    label 'python'

    publishDir "${params.result}/dev/1_Kraken2", mode: 'copy',
        pattern: "*_familyBarplot.tsv"
    publishDir "${params.result}/1_Kraken2", mode: 'copy',
        pattern: "*_familyBarplot.png"

    input:
        tuple val(sample_id), path(mpa_modif), path(total_reads)

    output:
        tuple val(sample_id),
            path("${sample_id}_familyBarplot.tsv"),
            path("${sample_id}_familyBarplot.png")

    script:
    """
    total=\$(cat ${total_reads} | head -n 1)

    assembly_mlst_mpa_family_barplot.py \
        ${mpa_modif} \
        \$total \
        ${sample_id}_familyBarplot.tsv \
        ${sample_id}_familyBarplot.png
    """
}

/*
* Calculate ratio nb_reads_interest / nb_reads_total
* Input   : Kraken2 MPA taxonomic profile
* Output  : TSV table (sample_ID, ratio_target, ratio_legio, commentary)
* Purpose : fast check if strain of interest or Legio in samples
*/
process RATIO_KRAKEN2 {
    label 'python'
    publishDir "${params.result}/dev/1_Kraken2", mode: 'copy'

    input:
        tuple val(sample_id), path(mpa), path(total_reads)

    output:
        tuple val(sample_id), path("${sample_id}_ratio.tsv")

    script:
    """
    total=\$(cat ${total_reads} | head -n 1)

    assembly_mlst_kraken2_ratio.py \
        ${sample_id} \
        \$total \
        ${mpa} \
        ${sample_id}_ratio.tsv \
        "${params.spp_target}" \
        ${params.min_ratio_target} \
        ${params.min_ratio_legio} \
        ${params.min_ratio_legia}
    """
}

/*
* Merge ratio Kraken2
* Input   : TSV table (sample_ID, ratio_target, ratio_legio, commentary)
* Output  : TSV table (sample_ID, ratio_target, ratio_legio, commentary)
* Purpose : create only one table with every samples
*/
process MERGE_RATIO_KRAKEN2 {
    label 'python'
    publishDir "${params.result}/dev/Summary", mode: 'copy'

    input:
        path(ratio_files)

    output:
        path("Summary_ratio.tsv")

    script:
    def target = params.spp_target.toString().replaceAll(/\s+/, '_')
    """
    printf "Sample_ID\t${target}_percent\tLegionella_spp_percent\tKraken2_results\n" \
        > Summary_ratio.tsv

    cat ${ratio_files} >> Summary_ratio.tsv
    """
}


// -----------------------------------------------------------------------------
/*
* Derive Legionella pneumophila MLST profile using El Gato
* Input   : paired-end reads (R1/R2 FASTQ)
* Output  : MLST results in TSV format (stdout redirected to file)
* Purpose : compute sequence type (ST) and allele calls per sample_id from reads
*/
process MLST_ELGATO {
    label 'elgato'

    publishDir "${params.result}/dev/2_ElGato", mode: 'copy',
        pattern: "*.csv"
    publishDir "${params.result}/2_ElGato", mode: 'copy',
        pattern: "*_*"

    input:
        tuple val(sample_id), val(r1), val(r2)

    output:
        tuple val(sample_id), path("${sample_id}_MLST.tsv"), emit: mlst
        tuple val(sample_id), path("${sample_id}.csv"), emit: fastfinder
        tuple val(sample_id), path("${sample_id}_reads")

    script:
    """
    el_gato.py \
        --read1 ${r1} \
        --read2 ${r2} \
        --sample ${sample_id} \
        --threads ${task.cpus} \
        --depth ${params.elgato_depth} \
        --out ${sample_id}_reads \
        -w \
        > ${sample_id}_MLST.tsv
    
    sed "s/${sample_id}/${sample_id},${params.fastfinder_value}/g" \
        ${sample_id}_MLST.tsv \
        | sed 's/\\t/,/g' \
        | sed 's/,,/,Nvx,/g' \
        | sed 's/,0,/,FAILED,/g' \
        | sed 's/,999,/,Nvx,/g' \
        > ${sample_id}.csv
    """
}

/*
* Merge MLST CSV or TSV outputs generated by El Gato
* Input   : multiple per-sample MLST CSV or TSV files
* Output  : single merged and sorted MLST table (CSV or TSV)
* Purpose : aggregate sequence type (ST) and allele calls across all sample_id
*/
process MERGE_ELGATO {
    label 'elgato'

    publishDir "${params.result}/dev/Summary", mode: 'copy',
        pattern: "*.tsv"
    publishDir "${params.result}", mode: 'copy',
        pattern: "*.csv"

    input:
        path(mlst_files)
        path(fastfinder_files)
    
    output:
        path("MLST_ElGatoResults_${params.suffix}.tsv"), emit: mlst
        path("Fastfinder_ElGatoResults_${params.suffix}.csv"), emit: fastfinder

    script:
    """
    printf "Sample_ID\tST\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA\n" \
        > MLST_ElGatoResults_${params.suffix}.tsv

    printf "Sample ID,${params.fastfinder_desc},ST,flaA,pilE,asd,mip,mompS,proA,neuA\n" \
        > Fastfinder_ElGatoResults_${params.suffix}.csv

    if [ ${mlst_files.size()} -ne 0 ]; then
        cat ${mlst_files} \
            | sort -t\$'\\t' -k1,1 \
            >> MLST_ElGatoResults_${params.suffix}.tsv

        cat ${fastfinder_files} \
            | sort -t',' -k1,1 \
            >> Fastfinder_ElGatoResults_${params.suffix}.csv
    fi
    """
}


// -----------------------------------------------------------------------------
/*
* Remove Windows carriage return characters from FASTA reference
* Input   : reference FASTA (CRLF possible)
* Output  : cleaned FASTA (Unix format)
* Purpose : ensure compatibility with downstream bioinformatics tools
*/
process CLEAN_REFERENCE_FASTA {
    label 'minimap2'

    input:
        path(minimap_ref)

    output:
        path("ref.fa")

    script:
    """
    sed 's/\r\$//' ${minimap_ref} > ref.fa
    """
}

/*
* Index reference genome using samtools faidx
* Input   : reference FASTA
* Output  : FASTA index (.fai)
* Purpose : enable fast random access to reference sequences for bcftools and alignment tools
*/
process INDEX_REFERENCE_FASTA {
    label 'samtools'

    input:
        path(ref_file)

    output:
        path("${ref_file}.fai")

    script:
    """
    samtools faidx ${ref_file}
    """
}

/*
* Align paired-end reads to a reference genome using Minimap2 with short-read optimized parameters
* Input   : paired-end FASTQ files (R1, R2) + reference FASTA genome
* Output  : SAM alignment file + cleaned reference FASTA file (Windows EOL removed)
* Purpose : generate read alignments against reference for downstream BAM conversion and variant calling
*/
process ALIGN_MINIMAP2 {
    label 'minimap2'
    // no PublishDir for .sam

    input:
        tuple val(sample_id), val(r1), val(r2)
        path(ref_file)

    output:
        tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    minimap2 \
        -t ${task.cpus} \
        -a \
        --sr \
        --frag=${params.minimap_frag} \
        -F${params.minimap_optF} \
        -k${params.minimap_optk} \
        -w${params.minimap_optw} \
        -A${params.minimap_optA} \
        -B${params.minimap_optB} \
        -O${params.minimap_optO} \
        -E${params.minimap_optE} \
        -r${params.minimap_optr} \
        -p${params.minimap_optp} \
        -N${params.minimap_optN} \
        -f${params.minimap_optf} \
        -n${params.minimap_optn} \
        -m${params.minimap_optm} \
        -s${params.minimap_opts} \
        -g${params.minimap_optg} \
        --heap-sort=${params.minimap_optheap} \
        --secondary=${params.minimap_optsec} \
        ${ref_file} \
        ${r1} \
        ${r2} \
        > ${sample_id}.sam
    """
}

/*
* Convert SAM alignments to sorted, indexed BAM format using Samtools
* Input   : SAM file + reference FASTA file
* Output  : sorted BAM file + BAM index (.bai) + reference FASTA file
* Purpose : prepare aligned reads in indexed binary format for efficient variant calling
*/
process TO_BAM_SAMTOOLS {
    label 'samtools'
    publishDir "${params.result}/dev/3_AMR", mode: 'copy'

    input:
        tuple val(sample_id), path(sam_file)
    
    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    samtools view -@ ${task.cpus} -b ${sam_file} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}

/*
* Call genomic variants using FreeBayes from aligned BAM reads
* Input   : sorted BAM file + BAM index + reference FASTA file
* Output  : raw VCF file (initial variant calls) + reference FASTA file
* Purpose : detect SNPs and small variants using haplotype-based Bayesian variant calling
*/
process VARCALL_FREEBAYES {
    label 'freebayes'
    publishDir "${params.result}/dev/3_AMR", mode: 'copy'

    input:
        tuple val(sample_id), path(bam_file), path(bai_file)
        path(ref_file)

    output:
        tuple val(sample_id), path("${sample_id}.init.vcf")

    script:
    """
    freebayes \
        --targets ${params.freeb_targets} \
        --theta ${params.freeb_theta} \
        --ploidy ${params.freeb_ploidy} \
        --report-all-haplotype-alleles \
        --pooled-continuous \
        --use-best-n-alleles ${params.freeb_best_n} \
        --allele-balance-priors-off \
        --haplotype-length ${params.freeb_haplo_len} \
        --use-duplicate-reads \
        --genotyping-max-iterations ${params.freeb_max_it} \
        --genotyping-max-banddepth ${params.freeb_max_dep} \
        --min-mapping-quality ${params.freeb_min_mapqual} \
        --min-base-quality ${params.freeb_min_basequal} \
        -F 0.05 \
        -C ${params.freeb_min_var} \
        --min-coverage ${params.freeb_min_dep} \
        -f ${ref_file} \
        -b ${bam_file} \
        > "${sample_id}.init.vcf"
    """
}

/*
* Filter variants using bcftools based on quality and allele frequency thresholds
* Input   : raw VCF
* Output  : filtered VCF (high-confidence variant set)
* Purpose : remove low-quality variants before normalization
*/
process FILTER_BCFTOOLS {
    label 'bcftools'
    publishDir "${params.result}/dev/3_AMR", mode: 'copy'

    input:
        tuple val(sample_id), path(init_vcf)

    output:
        tuple val(sample_id), path("${sample_id}.filt.vcf")

    script:
    """
    bcftools filter \
        -i 'AF >= ${params.bcf_min_freq} && \
            SAF >= ${params.bcf_saf} && \
            SAR >= ${params.bcf_sar} && \
            INFO/QA >= ${params.bcf_qa}' \
        -o ${sample_id}.filt.vcf \
        ${init_vcf}
    """
}

/*
* Normalize and decompose variants against reference genome using bcftools norm
* Input   : filtered VCF + reference genome
* Output  : final normalized VCF (biallelic, left-aligned, standardized)
* Purpose : ensure canonical variant representation for downstream AMR analysis
*/
process NORM_BCFTOOLS {
    label 'bcftools'
    publishDir "${params.result}/dev/3_AMR", mode: 'copy'

    input:
        tuple val(sample_id), path(filt_vcf)
        path(ref_file)

    output:
        tuple val(sample_id), path("${sample_id}.norm.vcf")

    script:
    """
    bcftools norm \
        -m -any \
        -f ${ref_file} \
        ${filt_vcf} \
        -o ${sample_id}.norm.vcf 
    """
}

/*
* Annotate genomic variants using SnpEff
* Input   : filtered VCF file (high-confidence variants)
* Output  : annotated VCF + HTML report + gene-level summary table
* Purpose : generate per-gene mutation impact statistics for downstream AMR analysis
*/
process IMPACT_SNPEFF {
    label 'snpeff'

    publishDir { type == 'AMR'
        ? "${params.result}/3_AMR"
        : "${params.result}/3_Other"
    }, mode: 'copy', pattern: "*.SnpEff.vcf"

    publishDir { type == 'AMR'
        ? "${params.result}/3_AMR"
        : "${params.result}/3_Other"
    }, mode: 'copy', pattern: "*.SnpEff.html"

    publishDir { type == 'AMR'
        ? "${params.result}/dev/3_AMR"
        : "${params.result}/dev/3_Other"
    }, mode: 'copy', pattern: "*.SnpEff.genes.txt"

    input:
        tuple val(sample_id), path(final_vcf)
        val(type)

    output:
        tuple val(sample_id), val(type), path("${sample_id}.SnpEff.vcf"), emit : vcf
        tuple val(sample_id), val(type), path("${sample_id}.SnpEff.html"), emit: html
        tuple val(sample_id), val(type), path("${sample_id}.SnpEff.genes.txt"), emit: genes

    script:
    """
    snpEff \
        -v \
        -c ${ type == 'AMR' ? params.snpeff_amr_config : params.snpeff_other_config } \
        ${ type == 'AMR' ? params.snpeff_amr_scheme : params.snpeff_other_scheme } \
        ${final_vcf} \
        > ${sample_id}.SnpEff.vcf \
        -stats ${sample_id}.SnpEff.html
    """
}

/*
* Merge SnpEff gene-level annotation tables across samples
* Input   : list of SnpEff genes.txt files (per sample)
* Output  : unified gene impact table
* Purpose : aggregate per-sample mutation impacts at gene level
*/
process MERGE_IMPACT_SNPEFF {
    label 'python'
    publishDir "${params.result}/dev/Summary", mode: 'copy'

    input:
        path(genes_impact)
        val(type)

    output:
        path("SnpEff_impact_merged_${type}.tsv")

    script:
    """
    if [ ${genes_impact.size()} -eq 0 ]; then
        printf "Sample_ID\tAMR_Nb_Mutated_Genes\tAMR_Nb_with_Impact\tAMR_Nb_Non_Coding" \
            > SnpEff_impact_merged_${type}.tsv
    else
        assembly_mlst_merge_impact_snpeff.py \
            -o SnpEff_impact_merged_${type}.tsv \
            -t ${type} \
            ${genes_impact.flatten().join(' ')}
    fi
    """
}

/*
* Extract gene-level functional annotations from SnpEff VCF outputs
* Input   : SnpEff annotated VCF (.vcf or .vcf.gz)
* Output  : normalized TSV containing gene, transcript and variant effect annotations
* Purpose : generate structured AMR-ready mutation tables for downstream analyses
*/
process PARSE_SNPEFF_GENES {
    label 'python'

    publishDir { type == 'AMR'
        ? "${params.result}/3_AMR"
        : "${params.result}/3_Other"
    }, mode: 'copy'

    input:
        tuple val(sample_id), val(type), path(vcf)

    output:
        tuple val(sample_id), path("${sample_id}.SnpEff.genes.parsed.tsv")

    script:
    """
    assembly_mlst_parse_snpeff_genes.py \
        --input ${vcf} \
        --sample ${sample_id} \
        --output ${sample_id}.SnpEff.genes.parsed.tsv
    """
}


// -----------------------------------------------------------------------------
/*
* Perform de novo genome assembly using SPAdes on paired-end reads
* Input   : paired-end reads (R1/R2 FASTQ)
* Output  : contigs + scaffolds + logs and assembly metadata
* Purpose : generate high-quality bacterial genome assemblies (careful mode, auto coverage cutoff)
*/
process ASSEMBLY_SPADES {
    label 'spades'
    publishDir "${params.result}/dev/4_SPAdes", mode: 'copy'

    input:
        tuple val(sample_id), val(r1), val(r2)
    
    output:
        tuple val(sample_id), path("${sample_id}/${sample_id}.contigs.fasta"), emit: contigs
        tuple val(sample_id), path("${sample_id}/${sample_id}.scaffolds.fasta"), emit: scaffolds
        path("${sample_id}/spades.log")
        path("${sample_id}/params.txt")
        path("${sample_id}/assembly_graph_with_scaffolds.gfa")

    script:
    """
    mkdir -p ${sample_id}
    
    spades.py \
        -1 ${r1} \
        -2 ${r2} \
        -t ${task.cpus} \
        -m ${task.memory.toGiga()} \
        --cov-cutoff auto \
        --careful \
        -o ${sample_id}
    
    mv ${sample_id}/contigs.fasta ${sample_id}/${sample_id}.contigs.fasta
    mv ${sample_id}/scaffolds.fasta ${sample_id}/${sample_id}.scaffolds.fasta
    """
}

/*
* Filter SPAdes contigs shorter than 500 bp
* Input  : output contigs by SPAdes
* Output : XXX.contigsXXXb.fasta stored per sample
* Purpose: remove short contigs for downstream QC and analysis
*/
process FILTER_CONTIGS {
    label 'seqkit'
    publishDir "${params.result}/4_SPAdes/${sample_id}", mode: 'copy'

    input:
        tuple val(sample_id), path(contigs)

    output:
        tuple val(sample_id), path("${sample_id}.contigs${params.min_length_contig}b.fasta")

    script:
    """
    seqkit seq \
        -m ${params.min_length_contig} \
        ${contigs} \
        > ${sample_id}.contigs${params.min_length_contig}b.fasta
    """
}

/*
* Assess assembly quality using QUAST
* Input   : list of contigs FASTA files
* Output  : QUAST reports and extracted metrics table
* Purpose : assess assembly quality metrics (N50, GC%, contigs count, etc.)
*/
process QC_QUAST {
    label 'quast'
    publishDir "${params.result}/4_SPAdes", mode: 'copy'

    input:
        path(contigs_files)

    output:
        path("QC-Quast/*")
        path("QC-Quast/extract_quast.tsv"), emit: extract

    script:
    """
    mkdir -p QC-Quast

    quast.py \
        -t ${task.cpus} \
        ${contigs_files.join(' ')} \
        -o QC-Quast

    tail -n +2 QC-Quast/transposed_report.tsv | \
        awk -v OFS="\\t" '{print \$1"\\t"\$16"\\t"\$14"\\t"\$17"\\t"\$18"\\t"\$20"\\t"\$22}' | \
        sed -E 's/\\.fasta//g' \
        > QC-Quast/extract_quast.tsv
    """
}

/*
* Filter and format QUAST extract results across samples
* Input   : extract_quast.tsv file
* Output  : unified QC summary table
* Purpose : consolidate assembly metrics per sample
*/
process MERGE_QC_QUAST {
    label 'quast'
    publishDir "${params.result}/dev/Summary", mode: 'copy'

    input:
        path(quast_file)

    output:
        path("QUAST_summary.tsv")

    script:
    """
    awk -v suffix=".contigs${params.min_length_contig}b" '
    BEGIN {
        FS=OFS="\\t"
        print "Sample_ID\\tTotal_length\\tNumber_of_contigs\\tGC\\tN50\\tauN\\tL90"
    }

    \$1 ~ suffix "\$" {

        sample = \$1
        sub(suffix "\$", "", sample)

        print sample, \$2, \$3, \$4, \$5, \$6, \$7
    }
    ' ${quast_file} > QUAST_summary.tsv
    """
}


// -----------------------------------------------------------------------------
/*
* Run mompS typing from paired-end reads and contigs
* Input   : sample_id, r1, r2, contigs
* Output  : FastFinder formatted CSV + mompS result directory
* Purpose : standardized Legionella typing workflow
*/
process MLST_MOMPS {
    label 'momps'

    publishDir "${params.result}/dev/Summary", mode: 'copy',
        pattern: "*.tsv"
    publishDir "${params.result}/dev/5_mompS", mode: 'copy',
        pattern: "*.csv"
    publishDir "${params.result}/5_mompS", mode: 'copy',
        pattern: "*/**"

    input:
        tuple val(sample_id), val(r1), val(r2), path(contigs)
    
    output:
        tuple val(sample_id), path("${sample_id}.csv"), emit: fastfinder
        tuple val(sample_id), path("${sample_id}.tsv"), emit: mlst
        tuple val(sample_id), path("${sample_id}/*"), emit: momps

    script:
    """
    perl /opt/mompS/momps.pl \
        -f ${r1} \
        -r ${r2} \
        -a ${contigs} \
        -o ${sample_id} \
        -p ${sample_id}

    sed "s/${sample_id}/${sample_id},${params.fastfinder_value}/g" \
        ${sample_id}/${sample_id}.MLST_res.txt \
        | sed 's/\\t/,/g' \
        | sed 's/,,/,Nvx,/g' \
        | sed 's/,0,/,FAILED,/g' \
        | sed 's/,999,/,Nvx,/g' \
        > ${sample_id}.csv

    cp ${sample_id}/${sample_id}.MLST_res.txt ${sample_id}.tsv

    # Remove unnecessary intermediate alignment and BLAST index files
    find ${sample_id} -type f \\( \
        -name "*.sam" -o \
        -name "*.bam" -o \
        -name "*.bai" -o \
        -name "*.nhr" -o \
        -name "*.nin" -o \
        -name "*.nsq" \
    \\) -delete
    """
}

/*
* Merge MLST CSV outputs generated by mompS
* Input   : multiple per-sample MLST CSV files
* Output  : single merged and sorted MLST table (CSV)
* Purpose : aggregate sequence type (ST) and allele calls across all sample_id
*/
process MERGE_MOMPS {
    label 'momps'
    publishDir "${params.result}", mode: 'copy'

    input:
        path(mlst_files)
        path(fastfinder_files)
    
    output:
        path("Fastfinder_mompSResults_${params.suffix}.csv"), emit : fastfinder
        path("MLST_mompSResults_${params.suffix}.tsv"), emit : mlst

    script:
    """
    printf "Sample ID,${params.fastfinder_desc},flaA,pilE,asd,mip,mompS,proA,neuA,ST\n" \
        > Fastfinder_mompSResults_${params.suffix}.csv
    printf "Sample_ID\tMompS_flaA\tMompS_pilE\tMompS_asd\tMompS_mip\tMompS_mompS\tMompS_proA\tMompS_neuA\tMompS_ST\n" \
        > MLST_mompSResults_${params.suffix}.tsv

    if [ ${mlst_files.size()} -ne 0 ]; then
        cat ${fastfinder_files} | sort -k1,1 >> Fastfinder_mompSResults_${params.suffix}.csv
        cat ${mlst_files} | sort -k1,1 >> MLST_mompSResults_${params.suffix}.tsv
    fi
    """
}


// -----------------------------------------------------------------------------
/*
* Identify closest reference genomes using FastANI
* Input   : contigs FASTA file
* Output  : fastANI similarity report
* Purpose : estimate ANI against a reference genome collection for strain identification
*/
process STRAIN_FASTANI {
    label 'fastani'
    publishDir "${params.result}/6_fastANI", mode: 'copy'

    input:
        tuple val(sample_id), path(contigs)

    output:
        tuple val(sample_id), path("${sample_id}.tsv")

    script:
    """
    fastANI \
        -q ${contigs} \
        --rl ${params.fastani_genomes} \
        -t ${task.cpus} \
        -o "${sample_id}.tsv"
    
    genomes_dir=\$(dirname "${params.fastani_genomes}")
    sed -i -e "s|\${genomes_dir}/||g" ${sample_id}.tsv
    sed -i -e "s|\\.fasta||g" ${sample_id}.tsv
    sed -i -e "s|^|${params.suffix}\\t${sample_id}\\t|g" ${sample_id}.tsv
    """
}

/*
* Merge FastANI strain identification results
* Input   : list of FastANI TSV files (one per sample)
* Output  : consolidated summary table (Sample_ID, FastANI_strain, FastANI_value)
* Purpose : extract best ANI hit per sample
*/
process MERGE_STRAIN_FASTANI {
    label 'fastani'
    publishDir "${params.result}/dev/Summary", mode: 'copy'

    input:
        path(fastani_files)

    output:
        path("FastANI_summary.tsv")

    script:
    def threshold = params.fastani_min.toFloat()
    """
    echo -e "Sample_ID\tFastANI_strain\tFastANI_value" > FastANI_summary.tsv

    for f in ${fastani_files}; do
        read line < "\$f"

        sample_id=\$(awk '{print \$2}' <<< "\$line")
        strain_full=\$(awk '{print \$4}' <<< "\$line")
        ani=\$(awk '{print \$5}' <<< "\$line")

        strain=\$(echo "\$strain_full" | cut -d'/' -f2)

        awk -v ani="\$ani" -v thr=${threshold} -v strain="\$strain" -v id="\$sample_id" '
        BEGIN {
            if (ani < thr) {
                print id "\\tPotential new spp.\\t" ani
            } else {
                print id "\\t" strain "\\t" ani
            }
        }' >> FastANI_summary.tsv
    done
    """
}

// -----------------------------------------------------------------------------
/*
* Evaluation of core genome MultiLocus Sequence Typing
* Input   : strain + contigs files + 1 params.X_set of genes, ptf, nb
* Output  : allele TSV renamed per strain/nb
* Purpose : run chewBBACA AlleleCall per strain group
*/
process MLST_CHEWBBACA {
    label 'chewbbaca'

    publishDir "${params.result}/7_chewBBACA", mode: 'copy',
        pattern: "*genes.tsv"
    publishDir "${params.result}/dev/7_chewBBACA", mode: 'copy',
        pattern: "*results_*"

    input:
        tuple val(strain),
            path(contigs_files),
            path(gene),
            path(ptf),
            val(nb),
            val(uid)

    output:
        tuple val(strain),
            val(nb),
            path("${strain}_alleles_${nb}genes.tsv"),
            emit: mlst
        path("results_${strain}-${nb}genes/*")

    script:
    """
    mkdir -p results contigs

    # Rebuild a clean input directory
    cp ${contigs_files.join(' ')} contigs/

    chewBBACA.py AlleleCall \
        -i contigs \
        -g ${gene} \
        --ptf ${ptf} \
        -o results

    # Ensure single deterministic output file
    result_file=\$(find results -type f -path '*/results_alleles.tsv' | head -n 1)

    cp \$result_file ${strain}_alleles_${nb}genes.tsv

    sed -i -e "s/.contigs500b//g" ${strain}_alleles_${nb}genes.tsv
    sed -i -e "s/FILE/Sample_ID/g" ${strain}_alleles_${nb}genes.tsv
    sed -i -e "s/INF-//g" ${strain}_alleles_${nb}genes.tsv

    # Rename output folder
    result_folder=\$(ls results | head -n 1)
    mv results/\$result_folder results_${strain}-${nb}genes
    """
}

/*
* Extract specific genes from chewBBACA cgMLST results
* Input   : strain + nb genes + allele TSV
* Output  : specific allele TSV
* Purpose : extract params.alleles_set
*/
process EXTRACT_ALLELES {
    label 'chewbbaca'
    publishDir "${params.result}/7_chewBBACA", mode: 'copy'
    
    input:
        tuple val(strain),
            val(nb),
            path(allele_tsv)
    
    output:
        tuple val(strain),
            val(nb),
            path("${strain}_alleles_${nb}genes.filt.tsv")

    script:
    def allele_args = params.alleles_set.collect {
        "${it.id}:${it.name}"
    }.join(',')
    """
    assembly_mlst_extract_alleles.sh \
        "${allele_args}" \
        ${allele_tsv} \
        ${strain}_alleles_${nb}genes.filt.tsv
    """
}

/*
* Merge chewBBACA allele extraction results across strains
* Input   : list of allele TSV files (heterogeneous gene sets)
* Output  : unified allele matrix with all loci merged
* Purpose : build a complete presence/absence allele matrix
*/
process MERGE_EXTRACT_ALLELES {
    label 'python'
    publishDir "${params.result}/dev/Summary", mode: 'copy'
    
    input:
        path(allele_files)

    output:
        path("Alleles_extracted.tsv")

    script:
    """
    if [ ${allele_files.size()} -eq 0 ]; then
        printf "Sample_ID" \
            > Alleles_extracted.tsv
    else
        assembly_mlst_merge_extract_alleles.py \
            -o Alleles_extracted.tsv \
            ${allele_files.join(' ')}
    fi
    """
}

/*
* Generate a visualisation from chewBBACA results
* Input   : strain + nb genes + allele TSV
* Output  : Newick tree file
* Purpose : QC chewBBACA alleles results
*/
process CHEWBBACA_GRAPETREE {
    label 'grapetree'
    publishDir "${params.result}/8_GrapeTree", mode: 'copy'

    input:
        tuple val(strain),
            val(nb),
            path(allele_tsv)

    output:
        path("${strain}_alleles_${nb}genes*.nwk")

    script:
    """
    suffix=""
    case "${allele_tsv}" in
        *.filt.tsv) suffix=".filt" ;;
    esac
    
    n_profiles=\$(
        tail -n +2 ${allele_tsv} |
        cut -f2- |
        sort -u |
        wc -l
    )

    if [ "\$n_profiles" -lt 2 ]; then
        echo "Only one unique allelic profile found" \
            > ${strain}_alleles_${nb}genes\${suffix}.nwk
        exit 0
    fi

    grapetree \
        -p ${allele_tsv} \
        -m MSTreeV2 \
        > ${strain}_alleles_${nb}genes\${suffix}.nwk
    """
}

/*
* Generate a visualisation from El Gato or mompS results
* Input   : software_name + results CSV
* Output  : Newick tree file + metadata TSV
* Purpose : QC results
*/
process LP_GRAPETREE {
    label 'grapetree'
    publishDir "${params.result}/8_GrapeTree", mode: 'copy'

    input:
        val(software)
        path(results_tsv)

    output:
        path("Lp_${software}.nwk"), emit: tree
        path("Lp_${software}.metadata.tsv"), emit: metadata

    script:
    """
    assembly_mlst_grapetree_tsv.sh \
        ${results_tsv} \
        tmp.tsv \
        Lp_${software}.metadata.tsv
    
    n_profiles=\$(
        tail -n +2 tmp.tsv |
        cut -f2- |
        sort -u |
        wc -l
    )

    if [ "\$n_profiles" -lt 2 ]; then
        echo "Only one unique allelic profile found" \
            > Lp_${software}.nwk
        exit 0
    fi

    grapetree \
        -p tmp.tsv \
        -m MSTreeV2 \
        > Lp_${software}.nwk
    """
}

/*
* Merge user and generated ReporTree MLST allele tables.
* Input   : user MLST allele TSV + generated MLST allele TSV
* Output  : merged MLST allele TSV + warning report
* Purpose : retain only shared columns and combine MLST allele for ReporTree
*/
process MERGE_REPORTREE_TSV {
    label 'python'

    publishDir "${params.result}/dev/Rsync", mode: 'copy',
        pattern: "*genes_MLSTchewbbaca.tsv"
    publishDir "${params.result}/dev/9_ReporTree", mode: 'copy',
        pattern: "*"

    input:
        tuple val(strain),
            val(nb),
            path(allele_tsv)

    output:
        tuple val(strain),
            val(nb),
            path("${strain}_${nb}genes_MLST.merged.tsv"),
            emit: mlst
        path("${strain}_${nb}genes_MLST.warning.tsv")

    script:
    """
    previous_MLST="${params.rep_partition}/${strain}_${nb}genes_MLSTchewbbaca.tsv"

    if [ -f "\${previous_MLST}" ]; then
        assembly_mlst_merge_mlst_tsv.py \
            --tsv-user \${previous_MLST} \
            --tsv ${allele_tsv} \
            --output-name "${strain}_${nb}genes_MLST.merged.tsv" \
            --warning "${strain}_${nb}genes_MLST.warning.tsv"
    else
        cp ${allele_tsv} "${strain}_${nb}genes_MLST.merged.tsv"
        echo "No ChewBBACA MLST TSV file given." > "${strain}_${nb}genes_MLST.warning.tsv"
    fi

    cp "${strain}_${nb}genes_MLST.merged.tsv" "${strain}_${nb}genes_MLSTchewbbaca.tsv"
    """
}

/*
* Prepare inputs for ReporTree visualisation
* Input   : strain + nb genes + allele TSV
* Output  : strain + nb genes + allele TSV
* Purpose : MLST newick tree and strain clustering
*/
// WARNING : col1 of ${allele_tsv} must be in the same order as col1 of cgMLST0.tsv 
process CHEWBBACA_REPORTREE {
    label 'chewbbaca'
    publishDir "${params.result}/dev/9_ReporTree", mode: 'copy'

    input:
        tuple val(strain),
            val(nb),
            path(allele_tsv)

    output:
        tuple val(strain),
            val(nb),
            path("${strain}_${nb}genes/cgMLST0.tsv"),
            emit : mlst
        path("${strain}_${nb}genes/*")

    script:
    """
    chewBBACA.py ExtractCgMLST \
        -i ${allele_tsv} \
        -o ${strain}_${nb}genes \
        --t 0

    awk -F'\\t' -v OFS='\\t' '
    NR==FNR { if(FNR>1) c1[FNR]=\$1; next }
    FNR==1 { print; next }
    { \$1=c1[FNR]; print }
    ' "${allele_tsv}" ${strain}_${nb}genes/cgMLST0.tsv \
    > ${strain}_${nb}genes/cgMLST0.tmp

    mv ${strain}_${nb}genes/cgMLST0.tmp ${strain}_${nb}genes/cgMLST0.tsv
    """
}

/*
* Merge Elgato and user metadata
* Input   : Elgato metadata TSV + user metadata TSV
* Output  : Merged metadata TSV
* Purpose : Keep ID and ST from Elgato, rename Sample_ID to ID, and append matching user metadata
*/
process LP_MERGE_METADATA {
    label 'python'
    publishDir "${params.result}/dev/9_ReporTree", mode: 'copy'

    input:
        path(metadata_elgato)

    output:
        path("Metadata_${params.suffix}.merged.tsv")

    script:
    """
    python3 << 'EOF'
    import pandas as pd

    elgato = pd.read_csv("${metadata_elgato}", sep="\\t", dtype=str).fillna("NA")
    user = pd.read_csv("${params.rep_metadata}", sep="\\t", dtype=str).fillna("NA")

    elgato = elgato[["Sample_ID", "ST"]].rename(columns={"Sample_ID": "ID"})

    merged = elgato.merge(user, on="ID", how="left").fillna("NA")

    merged.to_csv(
        "Metadata_${params.suffix}.merged.tsv",
        sep="\\t",
        index=False
    )
    EOF
    """
}

/*
* Generate ReporTree visualisation
* Input   : strain + nb genes + cgMLST allele tsv + metadata tsv
* Output  : ReporTree outputs
* Purpose : MLST newick tree and strain clustering
*/
//TODO: modifier ligne commande avec Christophe
process VISU_REPORTREE {
    label 'reportree'

    publishDir "${params.result}/dev/Rsync", mode: 'copy',
        pattern: "*_*genes_partitions.tsv"
    publishDir "${params.result}/9_ReporTree", mode: 'copy',
        pattern: "*genes/**"

    input:
        tuple val(strain),
            val(nb),
            path(allele_tsv),
            path(cgmlst_ref)
        path(metadata_tsv)

    output:
        path("${strain}_${nb}genes/*"), emit: folder
        path("${strain}_${nb}genes_partitions.tsv"), emit: partition

    script:
    """
    mkdir -p ${strain}_${nb}genes
    n_profiles=\$(tail -n +2 ${allele_tsv} | cut -f2- | sort -u | wc -l)

    if [ "\$n_profiles" -lt 2 ]; then
        echo "${strain}_alleles_${nb}genes : Only one unique allelic profile found" \
            > ${strain}_${nb}genes/Warning_ReporTree.txt
        exit 0
    fi

    zoom_opts=""
    case "${params.rep_zoom}" in
        all)
            zoom_opts="--zoom-all --site-inclusion ${params.rep_site_inclusion}"
            ;;
        none)
            zoom_opts=""
            ;;
        analyse)
            sample_list=\$(tail -n +2 ${metadata_tsv} | cut -f1 | paste -sd "," -)
            if [ -n "\${sample_list}" ]; then
                zoom_opts="--sample_of_interest \${sample_list} --zoom-cluster-of-interest ${params.rep_interest} --site-inclusion ${params.rep_site_inclusion}"
            else
                zoom_opts=""
            fi
            ;;
        *)
            zoom_opts="--sample_of_interest ${params.rep_zoom} --zoom-cluster-of-interest ${params.rep_interest} --site-inclusion ${params.rep_site_inclusion}"
            ;;
    esac

    nomenclature_opts=""
    if [ -f "${params.rep_partition}/${strain}_${nb}genes_partitions.tsv" ]; then
        nomenclature_opts="--nomenclature-file ${params.rep_partition}/${strain}_${nb}genes_partitions.tsv"
    fi

    threshold_opts=""
    if [ ${params.rep_min_allele} != "none" ]; then
        threshold_opts="--threshold ${params.rep_min_allele}-${params.rep_max_allele}"
    fi

    reportree.py \
        --n_proc ${task.cpus} \
        -m ${metadata_tsv} \
        -a ${allele_tsv} \
        -out ${strain}${nb}genes \
        -l ${cgmlst_ref} \
        --loci-called ${params.rep_loci_called} \
        --analysis grapetree \
        --method MSTreeV2 \
        --columns_summary_report "ST,Year,Origin,Linked_to" \
        --partitions2report 'all' \
        --metadata2report ${params.rep_col_metadata} \
        \$nomenclature_opts \
        \$zoom_opts \
        \$threshold_opts

    mv ${strain}${nb}genes* ${strain}_${nb}genes/.
    cp ${strain}_${nb}genes/${strain}${nb}genes_partitions.tsv ${strain}_${nb}genes_partitions.tsv
    """
}


// -----------------------------------------------------------------------------
/*
* Information about summary generation
* Input   : sample_id + optional files
* Output  : per-sample summary TSV file (one line, no header)
* Purpose : aggregate heterogeneous results, 1 row per sample, using NA for missing data
*/
process ASSEMBLY_MLST_SUMMARY_TABLE  {
    label 'python'
    publishDir "${params.result}", mode: 'copy'

    input:
        path(summary_files)

    output:
        path("Summary_table.tsv")

    script:
    """
    assembly_mlst_summary_table.py \
        -o Summary_table.tsv \
        ${summary_files.join(' ')}
    """
}

/*
* Information about HTML report generation
* Input   : merged summary TSV (output of ASSEMBLY_MLST_SUMMARY_TABLE)
* Output  : HTML report grouped by FastANI_strain + AMR summary section
* Purpose : generate human-readable interactive report for QC / MLST / AMR overview
*/
process ASSEMBLY_MLST_SUMMARY_HTML  {
    label 'python'
    publishDir "${params.result}", mode: 'copy'

    input:
        path(summary_tsv)
        path(software_tsv)

    output:
        path("Summary_table.html")

    script:
    """
    assembly_mlst_summary_html.py \
        -i ${summary_tsv} \
        -o "Summary_table.html" \
        --software ${software_tsv} \
        --analyse_id ${params.analyse_ID} \
        --sequencing_id ${params.suffix}
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
    input:
        val(suffix)
        val(input_dir)
        val(result)

        val(paired_end)
        val(adapters)
        val(decontamination)
        val(downsampling)
        val(momps)
        val(snpeff_other)
        
        val(bbtools_downsampled)

        val(min_quality)
        val(min_length)

        val(bbwrap_ref)
        val(bbwrap_path)
        val(bbwrap_min_id)
        val(bbwrap_max_indel)
        val(bbwrap_bwr)
        val(bbwrap_bw)
        val(bbwrap_min_hits)
        val(bbwrap_qtrim)
        val(bbwrap_trimq)
        val(bbwrap_qin)

        val(kraken2_db)
        val(format_mpa)

        val(spp_target)
        val(min_ratio_target)
        val(min_ratio_legio)
        val(min_ratio_legia)
        val(elgato_depth)

        val(minimap_ref)
        val(minimap_frag)
        val(minimap_optF)
        val(minimap_optk)
        val(minimap_optw)
        val(minimap_optA)
        val(minimap_optB)
        val(minimap_optO)
        val(minimap_optE)
        val(minimap_optr)
        val(minimap_optp)
        val(minimap_optN)
        val(minimap_optf)
        val(minimap_optn)
        val(minimap_optm)
        val(minimap_opts)
        val(minimap_optg)
        val(minimap_optheap)
        val(minimap_optsec)

        val(freeb_targets)
        val(freeb_theta)
        val(freeb_ploidy)
        val(freeb_best_n)
        val(freeb_haplo_len)
        val(freeb_max_it)
        val(freeb_max_dep)
        val(freeb_min_mapqual)
        val(freeb_min_basequal)
        val(freeb_min_var)
        val(freeb_min_dep)

        val(bcf_min_freq)
        val(bcf_qa)

        val(snpeff_amr_config)
        val(snpeff_amr_scheme)
        val(snpeff_other_config)
        val(snpeff_other_scheme)

        val(min_length_contig)
        
        val(fastani_genomes)
        val(fastani_min)

        val(rep_metadata)
        val(rep_partition)
        val(rep_interest)
        val(rep_zoom)
        val(rep_site_inclusion)
        val(rep_min_allele)
        val(rep_max_allele)
        val(rep_loci_called)
        val(rep_col_metadata)

        val(lb_set_json)
        val(lp_set_json)
        val(alleles_set_json)

    output:
        path("pipeline_${suffix}.txt")

    script:
    """
    assembly_mlst_create_info.sh \
        "${suffix}" \
        "${input_dir}" \
        "${result}" \
        "${paired_end}" \
        "${adapters}" \
        "${decontamination}" \
        "${downsampling}" \
        "${momps}" \
        "${snpeff_other}" \
        "${bbtools_downsampled}" \
        "${min_quality}" \
        "${min_length}" \
        "${bbwrap_ref}" \
        "${bbwrap_path}" \
        "${bbwrap_min_id}" \
        "${bbwrap_max_indel}" \
        "${bbwrap_bwr}" \
        "${bbwrap_bw}" \
        "${bbwrap_min_hits}" \
        "${bbwrap_qtrim}" \
        "${bbwrap_trimq}" \
        "${bbwrap_qin}" \
        "${kraken2_db}" \
        "${format_mpa}" \
        "${spp_target}" \
        "${min_ratio_target}" \
        "${min_ratio_legio}" \
        "${min_ratio_legia}" \
        "${elgato_depth}" \
        "${minimap_ref}" \
        "${minimap_frag}" \
        "${minimap_optF}" \
        "${minimap_optk}" \
        "${minimap_optw}" \
        "${minimap_optA}" \
        "${minimap_optB}" \
        "${minimap_optO}" \
        "${minimap_optE}" \
        "${minimap_optr}" \
        "${minimap_optp}" \
        "${minimap_optN}" \
        "${minimap_optf}" \
        "${minimap_optn}" \
        "${minimap_optm}" \
        "${minimap_opts}" \
        "${minimap_optg}" \
        "${minimap_optheap}" \
        "${minimap_optsec}" \
        "${freeb_targets}" \
        "${freeb_theta}" \
        "${freeb_ploidy}" \
        "${freeb_best_n}" \
        "${freeb_haplo_len}" \
        "${freeb_max_it}" \
        "${freeb_max_dep}" \
        "${freeb_min_mapqual}" \
        "${freeb_min_basequal}" \
        "${freeb_min_var}" \
        "${freeb_min_dep}" \
        "${bcf_min_freq}" \
        "${bcf_qa}" \
        '${snpeff_amr_config}' \
        '${snpeff_amr_scheme}' \
        '${snpeff_other_config}' \
        '${snpeff_other_scheme}' \
        '${min_length_contig}' \
        '${fastani_genomes}' \
        '${fastani_min}' \
        '${rep_metadata}' \
        '${rep_partition}' \
        '${rep_interest}' \
        '${rep_zoom}' \
        '${rep_site_inclusion}' \
        '${rep_min_allele}' \
        '${rep_max_allele}' \
        '${rep_loci_called}' \
        '${rep_col_metadata}' \
        '${lb_set_json}' \
        '${lp_set_json}' \
        '${alleles_set_json}'
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

process BBTOOLS_INFO {
    label 'bbtools'

    input:
        path(file) 

    output: 
        path("bbtools_${params.suffix}.txt")

    script:
    """
    software_track_file="bbtools_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "BBTOOLS VERSION" >> \$software_track_file
    bbmap.sh --version >> \$software_track_file 2>&1 || true
    """
}

process SEQKIT_INFO {
    label 'seqkit'

    input:
        path(file) 

    output: 
        path("seqkit_${params.suffix}.txt")

    script:
    """
    software_track_file="seqkit_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "SEQKIT VERSION" >> \$software_track_file
    seqkit version >> \$software_track_file || true
    """
}

process KRAKEN2_INFO {
    label 'kraken2'

    input:
        path(file) 

    output: 
        path("kraken2_${params.suffix}.txt")

    script:
    """
    software_track_file="kraken2_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "KRAKEN2 VERSION" >> \$software_track_file
    kraken2 --version >> \$software_track_file || true
    """
}

process PYTHON_INFO {
    label 'python'

    input:
        path(file) 

    output: 
        path("python_${params.suffix}.txt")

    script:
    """
    software_track_file="python_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "PYTHON PACKAGES VERSION" >> \$software_track_file
    python3 --version >> \$software_track_file || true
    python3 -c "import numpy, pandas, matplotlib; \
        print('numpy=='+numpy.__version__); \
        print('pandas=='+pandas.__version__); \
        print('matplotlib=='+matplotlib.__version__)" \
        >> \$software_track_file || true
    """
}

process KRONA_INFO {
    label 'krona'

    input:
        path(file) 

    output: 
        path("krona_${params.suffix}.txt")

    script:
    """
    software_track_file="krona_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "KRONA VERSION" >> \$software_track_file
    ktImportText 2>&1 | grep -oE 'KronaTools [0-9.]+' | awk '{print \$2}' >> \$software_track_file
    """
}

process ELGATO_INFO {
    label 'elgato'

    input:
        path(file) 

    output: 
        path("elgato_${params.suffix}.txt")

    script:
    """
    software_track_file="elgato_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "EL GATO VERSION" >> \$software_track_file
    el_gato.py --version >> \$software_track_file || true
    """
}

process MINIMAP_INFO {
    label 'minimap2'

    input:
        path(file) 

    output: 
        path("minimap_${params.suffix}.txt")

    script:
    """
    software_track_file="minimap_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "MINIMAP2 VERSION" >> \$software_track_file
    minimap2 --version >> \$software_track_file || true
    """
}

process FREEBAYES_INFO {
    label 'freebayes'

    input:
        path(file) 

    output: 
        path("freebayes_${params.suffix}.txt")

    script:
    """
    software_track_file="freebayes_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "FREEBAYES VERSION" >> \$software_track_file
    freebayes --version >> \$software_track_file || true
    """
}

process BCFTOOLS_INFO {
    label 'bcftools'

    input:
        path(file) 

    output: 
        path("bcftools_${params.suffix}.txt")

    script:
    """
    software_track_file="bcftools_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "BCFTOOLS VERSION" >> \$software_track_file
    bcftools --version >> \$software_track_file || true
    """
}

process SNPEFF_INFO {
    label 'snpeff'

    input:
        path(file) 

    output: 
        path("snpeff_${params.suffix}.txt")

    script:
    """
    software_track_file="snpeff_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "SNPEFF VERSION" >> \$software_track_file
    snpEff -version >> \$software_track_file || true
    """
}

process FASTANI_INFO {
    label 'fastani'

    input:
        path(file) 

    output: 
        path("fastani_${params.suffix}.txt")

    script:
    """
    software_track_file="fastani_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "FASTANI VERSION" >> \$software_track_file
    fastANI --version &>> \$software_track_file || true
    """
}

process SPADES_INFO {
    label 'spades'

    input:
        path(file) 

    output: 
        path("spades_${params.suffix}.txt")

    script:
    """
    software_track_file="spades_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "SPADES VERSION" >> \$software_track_file
    spades.py --version >> \$software_track_file || true
    """
}

process SAMTOOLS_INFO {
    label 'samtools'

    input:
        path(file) 

    output: 
        path("samtools_${params.suffix}.txt")

    script:
    """
    software_track_file="samtools_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "SAMTOOLS VERSION" >> \$software_track_file
    samtools --version >> \$software_track_file || true
    """
}

process QUAST_INFO {
    label 'quast'

    input:
        path(file) 

    output: 
        path("quast_${params.suffix}.txt")

    script:
    """
    software_track_file="quast_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "QUAST VERSION" >> \$software_track_file
    quast.py --version >> \$software_track_file || true
    """
}

process MOMPS_INFO {
    label 'momps'

    input:
        path(file) 

    output: 
        path("momps_${params.suffix}.txt")

    script:
    """
    software_track_file="momps_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "MOMPS VERSION" >> \$software_track_file
    echo ${task.ext.version} >> \$software_track_file || true
    """
}

process CHEWBBACA_INFO {
    label 'chewbbaca'

    input:
        path(file) 

    output: 
        path("chewbbaca_${params.suffix}.txt")

    script:
    """
    software_track_file="chewbbaca_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "CHEWBBACA VERSION" >> \$software_track_file
    chewBBACA.py --version >> \$software_track_file || true
    """
}

process GRAPETREE_INFO {
    label 'grapetree'

    input:
        path(file) 

    output: 
        path("grapetree_${params.suffix}.txt")

    script:
    """
    software_track_file="grapetree_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "GRAPETREE VERSION" >> \$software_track_file
    echo ${task.ext.version} >> \$software_track_file || true
    """
}

process REPORTREE_INFO {
    label 'reportree'

    input:
        path(file) 

    output: 
        path("reportree_${params.suffix}.txt")

    script:
    """
    software_track_file="reportree_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "REPORTREE VERSION" >> \$software_track_file
    reportree.py --version >> \$software_track_file || true
    """
}

process PUBLISH_INFO {
    label 'fastqc'
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