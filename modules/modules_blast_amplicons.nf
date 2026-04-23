#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Processes for workflow_blast_amplicons.nf


// -----------------------------------------------------------------------------
/*
* Reads quality control
* Input   : paired-end FASTQ files (R1, R2)
* Output  : FastQC reports (HTML + ZIP)
* Purpose : assess sequencing quality before or after downstream processing
*/
process QC_FASTQC {
    label 'fastqc'
    publishDir "${params.result}/0_QC/${read_type}", mode: 'copy'

    input:
        val(read_type)
        tuple val(sample_id), val(r1), val(r2)

    output:
        tuple val(sample_id), path("*.zip"), emit: zip_files
        tuple val(sample_id), path("*.html"), emit: html_files

    script:
    """
    fastqc \
        "${r1}" "${r2}" \
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
    publishDir "${params.result}/0_QC/${read_type}", mode: 'copy'

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
* Output  : unpaired reads, quality-filtered paired FASTQ files + adaptor trimming
* Purpose : remove low-quality reads and ensure high-confidence amplicon pairs
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
    def adapter_option = params.adapter ? "" : "--disable_adapter_trimming"

    """
    fastp \
        -i "${r1}" \
        -I "${r2}" \
        -o "${sample_id}_trimR1.fastq.gz" \
        -O "${sample_id}_trimR2.fastq.gz" \
        --qualified_quality_phred ${params.min_quality} \
        --length_required ${params.min_length} \
        --thread ${task.cpus} \
        ${adapter_option}
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
* Format Kraken2 output for visualization
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
    blast_amplicons_mpa_modified.py ${mpa} ${sample_id}_mpaModif.tsv
    """
}

/*
* Krona visualization from Kraken2 MPA-style report formatted
* Input   : Kraken2 results formatted (.tsv)
* Output  : interactive Krona HTML visualization
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
        tuple val(sample_id), path("${sample_id}_totalreads.txt")

    script:
    """
    seqkit stats -T ${r1} | \
        cut -f4 | \
        grep -v "num_seqs" | \
        sed 's/,//' \
        > "${sample_id}_totalreads.txt"
    """
}

/*
* Family-level abundance barplot from Kraken2 MPA report
* Input   : MPA-style taxonomic profile
* Output  : TSV table + barplot (top 15 families)
* Purpose : clean visualization of dominant bacterial families
*/
process MPA_FAMILY_BARPLOT {
    label 'python'
    // PublishDir in config file

    input:
        tuple val(sample_id), path(mpa_modif), path(total_reads)

    output:
        tuple val(sample_id),
            path("${sample_id}_familyBarplot.tsv"),
            path("${sample_id}_familyBarplot.png")

    script:
    """
    total=\$(cat ${total_reads} | head -n 1)

    blast_amplicons_mpa_family_barplot.py \
        ${mpa_modif} \
        \$total \
        ${sample_id}_familyBarplot.tsv \
        ${sample_id}_familyBarplot.png
    """
}

// -----------------------------------------------------------------------------
/*
* Paired-end reads merging 
* Input   : R1 and R2 FASTQ files
* Output  : merged sequences in FASTQ format
* Purpose : reconstruct full amplicon sequences and simplify downstream analysis
*/
process MERGE_FASTQ {
    label 'flash'
    publishDir "${params.result}/dev/2_Blast", mode: 'copy'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id), path("${sample_id}_flashOut/${sample_id}.extendedFrags.fastq.gz")

    script:
    def dovetail = params.dovetail_overlap ? 
        "--allow-outies" : 
        ""
    """
    flash ${r1} ${r2} \
        -d ${sample_id}_flashOut \
        -o ${sample_id} \
        -m ${params.min_overlap} \
        -M ${params.max_overlap}  \
        --compress \
        -t ${task.cpus} \
        ${dovetail}
    """
}

/*
* FASTQ to FASTA conversion
* Input   : merged FASTQ file
* Output  : merged sequences in FASTA format
* Purpose : reconstruct full amplicon sequences and simplify downstream analysis
*/
process FASTQ_TO_FASTA {
    label 'bbtools'
    publishDir "${params.result}/dev/2_Blast", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)

    output:
        tuple val(sample_id), path("${sample_id}_merged.fasta")

    script:
    """
    reformat.sh \
        in="${fastq}" \
        out="${sample_id}_merged.fasta" \
        -Xmx${task.memory.toGiga()}g
    """
}

/*
* Sequence dereplication with abundance tracking
* Input   : FASTA sequences
* Output  : dereplicated FASTA with size annotation
* Purpose : reduce redundancy while preserving sequence counts
*/
process DEREPLICATE_FASTA {
    label 'vsearch'
    publishDir "${params.result}/dev/2_Blast", mode: 'copy'

    input:
        tuple val(sample_id), path(fasta)

    output:
        tuple val(sample_id), path("${sample_id}_derep.fasta"), emit : fasta
        tuple path("${sample_id}_derep.uc"), path("${sample_id}_derep.log")

    script:
    """
    vsearch \
        --derep_fulllength "${fasta}" \
        --output "${sample_id}_derep.fasta" \
        --sizeout \
        --relabel Seq \
        --uc "${sample_id}_derep.uc" \
        --log "${sample_id}_derep.log" \
        --threads ${task.cpus}
    """
}

/*
* Total sequences count recovery from dereplicated FASTA
* Input   : dereplicated FASTA file with size annotations (size=)
* Output  : total number of sequences (sum of all size values)
* Purpose : reconstruct the original seq count prior to dereplication
*/
process COUNT_DEREP_FASTA {
    label 'python'
    publishDir "${params.result}/dev/2_Blast", mode: 'copy'

    input:
        tuple val(sample_id), path(fasta)

    output:
        tuple val(sample_id), path("${sample_id}_totalseq.txt")

    script:
    """
    grep -o "size=[0-9]*" ${fasta} | \
        sed 's/size=//' | \
        awk '{s+=\$1} END {print s}' \
        > "${sample_id}_totalseq.txt"
    """
}

// -----------------------------------------------------------------------------
/*
* Sequence alignment against reference database
* Input   : FASTA sequences
* Output  : BLAST tabular results (TSV)
* Purpose : assign taxonomy by sequence similarity search
*/
process BLASTN_FASTA {
    label 'blast'
    publishDir "${params.result}/dev/2_Blast", mode: 'copy'

    input:
        tuple val(sample_id), path(fasta)

    output:
        tuple val(sample_id), path("${sample_id}_blastresults.tsv")

    script:
    """
    blastn \
        -task megablast \
        -db "${params.blast_db}" \
        -query "${fasta}" \
        -out "${sample_id}_blastresults.tsv" \
        -outfmt "6 qseqid sseqid pident length qlen slen qcovhsp bitscore" \
        -num_threads ${task.cpus}
    """
}

// -----------------------------------------------------------------------------
/*
* BLAST filtering and taxonomic assignment with abundance tracking
* Input   : BLAST TSV results (dereplicated sequences with size in headers)
* Output  : table TSV (read_id, taxon, size)
* Purpose : assign one taxon per sequence and preserve abundance information
*/
process FILTER_BLASTN {
    label 'python'
    publishDir "${params.result}/2_Blast", mode: 'copy'

    input:
        tuple val(sample_id), path(blast)

    output:
        tuple val(sample_id), path("${sample_id}_strictfiltered.tsv"), emit:strict
        tuple val(sample_id), path("${sample_id}_loosefiltered.tsv"), emit:loose

    script:
    """
    blast_amplicons_filter_blast.py \
        "${blast}" \
        ${params.perc_id} \
        ${params.query_cov} \
        ${params.min_qlen} \
        ${params.delta_bitscore} \
        ${params.loose_id} \
        ${params.loose_cov} \
        ${params.loose_qlen} \
        "${sample_id}_strictfiltered.tsv" \
        "${sample_id}_loosefiltered.tsv"
    """
}

// -----------------------------------------------------------------------------
/*
* Sequence counting and visualization
* Input   : results table TSV (read_id, taxon, size)
* Output  : graphical summary
* Purpose : quantify taxa and generate interpretable plots
*/
process PLOT_BLASTFILT {
    label 'python'
    publishDir "${params.result}/2_Blast", mode: 'copy'

    input:
        val(filt_type)
        tuple val(sample_id), path(filtered), path(total_seq)

    output:
        tuple val(sample_id), path("${sample_id}_${filt_type}Barplot.png")

    script:
    """
    total=\$(cat ${total_seq})

    blast_amplicons_plot_blast.py \
        "${filtered}" \
        \${total} \
        "${sample_id}_${filt_type}Barplot.png"
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
        val(adapter)
        val(decontamination)
        val(downsampling)
        val(kraken2_assign)

        val(min_quality)
        val(min_length)

        val(bbwrap_ref)
        val(bbwrap_min_id)
        val(bbwrap_max_indel)
        val(bbwrap_bwr)
        val(bbwrap_bw)
        val(bbwrap_min_hits)
        val(bbwrap_qtrim)
        val(bbwrap_trimq)
        val(bbwrap_qin)
        val(bbwrap_path)

        val(bbtools_downsampled)

        val(kraken2_db)

        val(min_overlap)
        val(max_overlap)
        val(dovetail_overlap)

        val(blast_db)
        val(perc_id)
        val(loose_id)
        val(query_cov)
        val(loose_cov)
        val(min_qlen)
        val(loose_qlen)
        val(delta_bitscore)

    output:
        path("pipeline_${suffix}.txt")

    script:
    """
    blast_amplicons_create_info.sh \
        "${suffix}" \
        "${input_dir}" \
        "${result}" \
        "${paired_end}" \
        "${adapter}" \
        "${decontamination}" \
        "${downsampling}" \
        "${kraken2_assign}" \
        "${min_quality}" \
        "${min_length}" \
        "${bbwrap_ref}" \
        "${bbwrap_min_id}" \
        "${bbwrap_max_indel}" \
        "${bbwrap_bwr}" \
        "${bbwrap_bw}" \
        "${bbwrap_min_hits}" \
        "${bbwrap_qtrim}" \
        "${bbwrap_trimq}" \
        "${bbwrap_qin}" \
        "${bbwrap_path}" \
        "${bbtools_downsampled}" \
        "${kraken2_db}" \
        "${min_overlap}" \
        "${max_overlap}" \
        "${dovetail_overlap}" \
        "${blast_db}" \
        "${perc_id}" \
        "${loose_id}" \
        "${query_cov}" \
        "${loose_cov}" \
        "${min_qlen}" \
        "${loose_qlen}" \
        "${delta_bitscore}"
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

process FLASH_INFO {
    label 'flash'

    input:
        path(file) 

    output: 
        path("flash_${params.suffix}.txt")

    script:
    """
    software_track_file="flash_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "FLASH VERSION" >> \$software_track_file
    flash -v >> \$software_track_file || true
    """
}

process VSEARCH_INFO {
    label 'vsearch'

    input:
        path(file) 

    output: 
        path("vsearch_${params.suffix}.txt")

    script:
    """
    software_track_file="vsearch_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "VSEARCH VERSION" >> \$software_track_file
    vsearch --version >> \$software_track_file || true
    """
}

process BLAST_INFO {
    label 'blast'

    input:
        path(file) 

    output: 
        path("blast_${params.suffix}.txt")

    script:
    """
    software_track_file="blast_${params.suffix}.txt"
    cat $file > \$software_track_file

    echo "" >> \$software_track_file

    echo "BLAST+ VERSION" >> \$software_track_file
    blastn -version >> \$software_track_file || true
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
