#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Processes for workflow_elgato_sbt.nf


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
    publishDir "${params.result}/0_QC", mode: 'copy'

    input:
        val(read_type)
        path(fastqc_zip)

    output:
        path("${read_type}_multiQC_report.html")

    script:
    """
    multiqc ${fastqc_zip} \
        --filename "${read_type}_multiQC_report.html" \
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
    def adapter_option = params.adapters ? 
        "--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" : 
        "--disable_adapter_trimming"

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
    elgato_sbt_mpa_modified.py ${mpa} ${sample_id}_mpaModif.tsv
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

    elgato_sbt_mpa_family_barplot.py \
        ${mpa_modif} \
        \$total \
        ${sample_id}_familyBarplot.tsv \
        ${sample_id}_familyBarplot.png
    """
}


// -----------------------------------------------------------------------------
/*
* Derive Legionella pneumophila MLST profile using el_gato
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
* Merge MLST TSV outputs generated by el_gato
* Input   : multiple per-sample MLST TSV files
* Output  : single merged and sorted MLST table (TSV)
* Purpose : aggregate sequence type (ST) and allele calls across all sample_id
*/
process MERGE_ELGATO {
    label 'elgato'

    publishDir "${params.result}", mode: 'copy',
        pattern: "*.csv"
    publishDir "${params.result}/2_ElGato", mode: 'copy',
        pattern: "*.tsv"

    input:
        path(mlst_files)
        path(fastfinder_files)
    
    output:
        path("MLST_ElGatoResults_${params.suffix}.tsv"), emit: mlst
        path("Fastfinder_ElGatoResults_${params.suffix}.csv"), emit: fastfinder

    script:
    """
    printf "Sample_ID\tST\tflaA\tpilE\tasd\tmip\tmompS\tproA\tneuA\n" \
        > "MLST_ElGatoResults_${params.suffix}.tsv"
    cat ${mlst_files}  \
        | sort -t\$'\\t' -k1,1 \
        >> "MLST_ElGatoResults_${params.suffix}.tsv"

    
    printf "Sample ID,${params.fastfinder_desc},ST,flaA,pilE,asd,mip,mompS,proA,neuA\n" \
        > "Fastfinder_ElGatoResults_${params.suffix}.csv"
    cat ${fastfinder_files} \
        | sort -t',' -k1,1 \
        >> "Fastfinder_ElGatoResults_${params.suffix}.csv"
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

        val(elgato_depth)

    output:
        path("pipeline_${suffix}.txt")

    script:
    """
    elgato_sbt_create_info.sh \
        "${suffix}" \
        "${input_dir}" \
        "${result}" \
        "${paired_end}" \
        "${adapters}" \
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
        "${elgato_depth}"
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

process ELGATO_INFO{
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