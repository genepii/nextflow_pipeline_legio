# Configuration file

This configuration file defines:

- the execution environment (Singularity, SLURM profile, resources),
- the computational resources allocated to each pipeline process,
- the analysis parameters,
- the locations of input, output and reference files.

---

> [!IMPORTANT]
> All parameter values defined in this configuration file are used directly during pipeline execution when the workflow is launched through SLURM.
>
> The values specified in this file determine the behavior and resources allocated to the pipeline. Please check them carefully before execution.
>
> Parameters should only be modified by changing their values. Do not rename, remove, or modify parameter names, as they are required by the pipeline configuration.
>
> Before launching the pipeline, ensure that:
> - parameter names remain unchanged;
> - only parameter values are adapted according to the required execution settings;
> - resource allocation and execution options are compatible with the SLURM environment.
>
> This configuration file is the single source of truth for all parameters used during pipeline execution.

---

# Execution environment

## Singularity configuration

| Parameter | Description |
|-----------|-------------|
| `enabled` | Enables execution inside Singularity containers. |
| `autoMounts` | Automatically mounts host directories into the container. |
| `runOptions` | Additional options passed to Singularity at runtime. |
| `envWhitelist` | Environment variables allowed to be inherited by the container. |

---

# Execution profile

The `local` profile runs the workflow directly on the local machine.

| Parameter | Value | Description |
|-----------|------:|-------------|
| `process.executor` | `local` | Local execution. |
| `process.cpus` | `2` | Default number of CPUs per process. |
| `process.memory` | `200 GB` | Maximum memory available. |
| `process.maxForks` | `1` | Maximum number of simultaneously running processes. |

The `slurm` profile runs the workflow on the local machine via SLURM.

| Parameter | Value | Description |
|-----------|------:|-------------|
| `process.executor` | `slurm` | SLURM execution. |
| `process.queue` | `diag_iai` | SLURM partition. |
| `process.time` | `24h` | Maximum running time. |
| `process.cpus` | `4` | Default number of CPUs per process. |
| `process.memory` | `200 GB` | Maximum memory available. |
| `process.maxForks` | `1` | Maximum number of simultaneously running processes. |

---

# Global process settings

| Parameter | Description |
|-----------|-------------|
| `scratch = true` | Executes each process inside a temporary scratch directory managed by Nextflow. |

---

# Process resources

Each process label defines the container image and computational resources used.

| Process | Container | CPUs | Max concurrent tasks | Initial memory | Max retries |
|---------|-----------|-----:|---------------------:|---------------:|------------:|
| FastQC | `fastqc_v0.12.1.sif` | 1 | 5 | 1 GB | 5 |
| MultiQC | `multiqc_v1.33.sif` | 1 | 1 | 1 GB | 5 |
| Fastp | `fastp_v1.3.2.sif` | 2 | 4 | 2 GB | 4 |
| BBTools | `bbtools_v39.81.sif` | 8 | 2 | 90 GB | 2 |
| SeqKit | `seqkit_v2.13.0.sif` | 1 | 2 | 1 GB | 5 |
| Kraken2 | `kraken2_v2.17.1.sif` | 6 | 2 | 90 GB | 2 |
| Python | `python_plot_parsing_3.11.sif` | 4 | 4 | 1 GB | 5 |
| Krona | `kronatools_2.8.1.sif` | 4 | 4 | 1 GB | 5 |
| ElGato | `elgato_1.22.0.sif` | 4 | 2 | 50 GB | 3 |
| Minimap2 | `minimap2_2.31.sif` | 8 | 2 | 90 GB | 2 |
| Samtools | `samtools_1.23.1.sif` | 4 | 2 | 50 GB | 3 |
| FreeBayes | `freebayes_1.3.10.sif` | 4 | 2 | 50 GB | 3 |
| BCFtools | `bcftools_1.23.1.sif` | 4 | 2 | 50 GB | 3 |
| SnpEff | `snpeff_5.4c.sif` | 4 | 2 | 50 GB | 3 |
| MOMP-S | `momps_2026_05_22.sif` | 4 | 2 | 50 GB | 3 |
| SPAdes | `spades_4.2.0.sif` | 8 | 2 | 90 GB | 2 |
| QUAST | `quast_5.3.0.sif` | 1 | 5 | 1 GB | 5 |
| FastANI | `fastani_1.34.sif` | 4 | 2 | 50 GB | 3 |
| chewBBACA | `chewbbaca_3.3.10.sif` | 4 | 2 | 50 GB | 3 |
| Grapetree | `grapetree_2.1.sif` | 1 | 5 | 1 GB | 5 |
| ReporTree | `reportree_2.6.1.sif` | 8 | 2 | 90 GB | 2 |

All processes use the following retry strategy:

- Retry only when the exit status is `137` (typically memory exhaustion).
- Memory is doubled after each retry.
- Other failures terminate the process.

---

# Pipeline parameters

## General parameters

| Parameter | Description |
|-----------|-------------|
| `analyse_id` | Analysis identifier automatically generated from the execution date (`YYYYMMDD`) (STR). Used to create the analysis output directory name. |
| `suffix` | Sequencing run identifier. Used as a subdirectory name to locate input FASTQ files and organize results (STR). |
| `input_dir_prefix` | Base directory containing sequencing run folders. The final input directory is automatically built as `${input_dir_prefix}/Legionella-Amplicons-${suffix}` (PATH). |
| `output_dir_prefix` | Base directory where final analysis outputs are stored. The final output directory is automatically built as `${output_dir_prefix}/${suffix}/${analyse_id}_Assembly-MLST` (PATH). |
| `save_dir_prefix` | Base directory used to store raw FASTQ files from the sequencing run. The final save directory is automatically built as `${save_dir_prefix}/${suffix}` (PATH). |
| `tmp_dir_prefix` | Base directory used for temporary files generated during execution. The final temporary directory is automatically built as `${tmp_dir_prefix}/${suffix}` (PATH). |
| `work_dir_prefix` | Base directory used for Nextflow work files and intermediate results. The final work directory is automatically built as `${work_dir_prefix}/${suffix}/${analyse_id}_Assembly-MLST` (PATH). |
| `paired_end` | Enables paired-end analysis (`true`) or single-end analysis (`false`) (true/false). |
| `adapters` | Enables adapter trimming using Fastp (true/false). |
| `decontamination` | Enables host-read removal using BBWrap (true/false). |
| `downsampling` | Enables read downsampling before downstream analyses (true/false). |
| `momps` | Enables MOMP-S analysis (true/false). |
| `snpeff_other` | Enables annotation using an additional SnpEff database (true/false). |
| `rep_zoom` | ReporTree zoom mode: analyze samples from metadata, all samples, no zoom, or selected sample IDs (`analyse`, `all`, `none`, or one/multiple sample IDs such as `SampleA` or `SampleA,SampleB`). |

---

# Reformat / downsampling

| Parameter | Description |
|-----------|-------------|
| `bbtools_downsampled` | Fraction of reads retained after downsampling, e.g. `0.5` keeps 50% of reads (FLOAT). |

---

# FASTP parameters

| Parameter | Description |
|-----------|-------------|
| `min_quality` | Minimum Phred quality score required to retain bases (INT). |
| `min_length` | Minimum read length after trimming (INT). |

---

# BBWrap decontamination

| Parameter | Description |
|-----------|-------------|
| `bbwrap_ref` | Reference FASTA file. If empty, the indexed database specified by `bbwrap_path` is used (PATH). |
| `bbwrap_path` | Path to the indexed BBMap reference database (PATH). |
| `bbwrap_min_id` | Minimum alignment identity (FLOAT). |
| `bbwrap_max_indel` | Maximum allowed indel length (INT). |
| `bbwrap_bwr` | Bandwidth ratio used during alignment (FLOAT). |
| `bbwrap_bw` | Alignment bandwidth (INT). |
| `bbwrap_min_hits` | Minimum number of matching seeds required (INT). |
| `bbwrap_qtrim` | Quality trimming mode (`r`, `l` or `rl`). |
| `bbwrap_trimq` | Quality threshold used for trimming (INT). |
| `bbwrap_qin` | Input FASTQ quality encoding, typically `33` (INT). |

---

# Kraken2 classification

| Parameter | Description |
|-----------|-------------|
| `kraken2_db` | Path to the Kraken2 database (PATH). |
| `format_mpa` | Produces MetaPhlAn-compatible output format. If disabled, MPA-related processes should also be removed from the workflow (true/false). |

---

# ElGato species identification

| Parameter | Description |
|-----------|-------------|
| `spp_target` | Kraken2 species name used as target organism (STR). |
| `min_ratio_target` | Minimum proportion of target species reads (FLOAT). |
| `min_ratio_legio` | Minimum proportion of reads assigned to the genus *Legionella* (FLOAT). |
| `min_ratio_legia` | Minimum proportion of reads assigned to the group *Legionellaceae* (FLOAT). |
| `elgato_depth` | Minimum sequencing depth required for allele calling (INT). |

---

# Minimap2 antibiotic resistance mapping

| Parameter | Description |
|-----------|-------------|
| `minimap_ref` | Reference FASTA containing resistance genes (PATH). |
| `minimap_opt*` | Minimap2 alignment options controlling sensitivity, scoring, filtering and output behaviour (see [manual](https://lh3.github.io/minimap2/minimap2.html)). |

---

# FreeBayes variant calling

| Parameter | Description |
|-----------|-------------|
| `freeb_targets` | BED file containing targeted resistance regions (PATH). |
| `freeb_theta` | Expected heterozygosity parameter (FLOAT). |
| `freeb_ploidy` | Assumed ploidy level (INT). |
| `freeb_best_n` | Maximum number of alternate alleles reported (INT). |
| `freeb_haplo_len` | Minimum haplotype length (INT). |
| `freeb_max_it` | Maximum iterations (INT). |
| `freeb_max_dep` | Maximum coverage depth considered (INT). |
| `freeb_min_mapqual` | Minimum mapping quality (INT). |
| `freeb_min_basequal` | Minimum base quality (INT). |
| `freeb_min_var` | Minimum number of supporting reads for variants (INT). |
| `freeb_min_dep` | Minimum depth required for variant calling (INT). |

---

# BCFtools filtering

| Parameter | Description |
|-----------|-------------|
| `bcf_min_freq` | Minimum allele frequency (FLOAT). |
| `bcf_saf` | Minimum number of supporting forward reads (INT). |
| `bcf_sar` | Minimum number of supporting reverse reads (INT). |
| `bcf_qa` | Minimum alignment quality (INT). |

---

# SnpEff annotation

| Parameter | Description |
|-----------|-------------|
| `snpeff_amr_config` | SnpEff configuration file for antimicrobial resistance annotation (PATH). |
| `snpeff_amr_scheme` | Annotation database name for resistance genes (STR). |
| `snpeff_other_config` | Additional SnpEff configuration file (PATH). |
| `snpeff_other_scheme` | Additional annotation database name (STR). |

---

# Assembly parameters

| Parameter | Description |
|-----------|-------------|
| `min_length_contig` | Minimum contig length retained after assembly (INT). |

---

# FastFinder metadata

| Parameter | Description |
|-----------|-------------|
| `fastfinder_desc` | Metadata field names for FastFinder, comma-separated values (STR). |
| `fastfinder_value` | Expected values associated with the selected metadata fields, comma-separated values (STR). |

---

# FastANI

| Parameter | Description |
|-----------|-------------|
| `fastani_genomes` | Reference genome list used for ANI comparison (PATH). |
| `fastani_min` | Minimum ANI threshold (INT). |

---

# ReporTree

| Parameter | Description |
|-----------|-------------|
| `rep_metadata` | Metadata file used by ReporTree (PATH). |
| `rep_partition` | ReporTree partitions directory containing partition files named according to the `${strain}_${nb}genes_partitions.tsv` format, e.g. `Lp_4707genes_partitions.tsv` (PATH). |
| `rep_interest` | Cluster threshold used to define the group of samples included in the zoom analysis (INT). |
| `rep_site_inclusion` | Minimum proportion of sites required for the zoom analysis (FLOAT). |
| `rep_min_allele` | Minimum allele threshold for clustering (INT or `none` if no threshold wanted). |
| `rep_max_allele` | Maximum allele threshold for clustering (INT). |
| `rep_loci_called` | Minimum proportion of called loci required (FLOAT). |
| `rep_col_metadata` | Metadata column used for grouping (STR). |

---

# chewBBACA and Grapetree

| Parameter | Description |
|-----------|-------------|
| `lb_set` | Configuration list for *Legionella longbeachae* cgMLST/wgMLST schemas. |
| `lp_set` | Configuration list for *Legionella pneumophila* cgMLST/wgMLST schemas. |
| `genes` | chewBBACA schema directory (PATH). |
| `ptf` | Protein translation file (PATH). |
| `lit` | Locus information table (PATH). |
| `nb` | Number of loci in the scheme (INT). |
| `grapetree` | Enables Grapetree analysis (true/false). |
| `reportree` | Enables ReporTree analysis (true/false). |

---

# Specific alleles

| Parameter | Description |
|-----------|-------------|
| `alleles_set` | List of specific alleles to extract and report. Allele identifiers must correspond to alleles available in the `{genes}` schema files. |
| `id` | Allele identifier (STR). |
| `name` | Display name associated with the allele (STR). |

---

# Working directory

| Parameter | Description |
|-----------|-------------|
| `workDir` | Nextflow working directory containing intermediate files (PATH). |
