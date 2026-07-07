# Configuration file

This configuration file defines:

- the execution environment (Singularity, local profile, resources),
- the computational resources allocated to each pipeline process,
- the default analysis parameters,
- the locations of input, output and reference files.

---

> [!IMPORTANT]
> Some parameter values shown in this configuration file are **default values only** and are **overridden by the corresponding `start_*.sh` launcher script**.
>
> These parameters should preferably be modified through the launcher script (or its command-line options when available), as the launcher defines the values actually used during pipeline execution.
>
> In the configuration file, this section begins with:
>
> ```groovy
> // -----------------------------------------------------------------------------
> // Notes : Settings below are set by the `start_XXX.sh` script;
> //     values shown are indicative only and apply if the script is run directly.
> // -----------------------------------------------------------------------------
> ```
>
> and ends with:
>
> ```groovy
> // -------------------------------------------------------------------------
> // Notes : Settings below are set by the `XXX.config` file.
> // -------------------------------------------------------------------------
> ```
>
> Parameters located between these two comment blocks are managed by the launcher script. All parameters below the second comment block are read directly from the configuration file.

---

# Execution environment

## User and temporary directory

| Parameter | Description |
|-----------|-------------|
| `user` | Current Linux username (`$USER`). If unavailable, defaults to `blast-user`. |
| `baseScratch` | Scratch directory used for temporary files and as the container home directory. |

---

## Singularity configuration

| Parameter | Description |
|-----------|-------------|
| `enabled` | Enables execution inside Singularity containers. |
| `autoMounts` | Automatically mounts host directories into the container. |
| `runOptions` | Additional options passed to Singularity at runtime. |

### `runOptions`

| Option | Description |
|--------|-------------|
| `--bind /srv/scratch/iai` | Makes the scratch filesystem available inside the container. |
| `--bind /srv/net/cluqumngs` | Makes the sequencing storage available inside the container. |
| `--env TMPDIR=${baseScratch}` | Uses the scratch directory for temporary files. |
| `--env HOME=${baseScratch}` | Sets the container home directory to the scratch directory. |

---

# Execution profile

The `local` profile runs the workflow directly on the local machine.

| Parameter | Value | Description |
|-----------|------:|-------------|
| `process.executor` | `local` | Local execution. |
| `process.cpus` | `2` | Default number of CPUs per process. |
| `process.memory` | `200 GB` | Maximum memory available. |
| `process.maxForks` | `1` | Maximum number of simultaneously running processes. |

---

# Global process settings

Before each process execution, the pipeline ensures that the temporary directory exists and exports it as `TMPDIR`.

| Parameter | Description |
|-----------|-------------|
| `beforeScript` | Creates `${baseScratch}` if necessary and exports it as the temporary directory used by all processes. |

---

# Process resources

Each process label defines the container image and computational resources used.

| Process | Container | CPUs | Max concurrent tasks | Initial memory | Memory increase | Retry policy | Maximum retries |
|---------|-----------|-----:|---------------------:|---------------:|-----------------|----------------|----------------:|
| FastQC | `fastqc_v0.12.1.sif` | 1 | 10 | 5 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 2 |
| MultiQC | `multiqc_v1.33.sif` | 1 | 1 | 15 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 2 |
| Fastp | `fastp_v1.3.2.sif` | 2 | 4 | 1 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 5 |
| BBTools | `bbtools_v39.81.sif` | 8 | 2 | 90 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 2 |
| SeqKit | `seqkit_v2.13.0.sif` | 1 | 2 | 1 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 5 |
| Kraken2 | `kraken2_v2.17.1.sif` | 8 | 2 | 90 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 2 |
| Python | `python_plot_3.11.sif` | 4 | 4 | 1 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 5 |
| KronaTools | `kronatools_2.8.1.sif` | 4 | 4 | 1 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 5 |
| FLASH | `flash_v1.2.11.sif` | 8 | 2 | 90 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 2 |
| BLAST+ | `blast_v2.17.0.sif` | 8 | 2 | 90 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 2 |
| VSEARCH | `vsearch_v2.30.6.sif` | 4 | 4 | 90 GB | Doubled after each retry | Retry on exit code `137` (out-of-memory) | 2 |

All processes use the following retry strategy:

- Retry only when the exit status is `137` (typically memory exhaustion).
- Memory is doubled after each retry.
- Other failures terminate the process.

---

# Pipeline parameters

The following parameters define the analysis to be performed.

## General parameters

| Parameter | Description |
|-----------|-------------|
| `analyse_ID` | Analysis identifier automatically generated from the current date (STR). |
| `suffix` | Sequencing run identifier. Usually corresponds to the sequencing date (`YYYYMMDD`) (STR). |
| `input_dir` | Directory containing the input FASTQ files (PATH). |
| `result` | Output directory where all analysis results are written (PATH). |
| `paired_end` | Enables paired-end analysis (`true`) or single-end analysis (`false`) (true/false). |
| `adapters` | Enables adapter trimming using Fastp (true/false). |
| `decontamination` | Enables host-read removal using BBWrap (true/false). |
| `downsampling` | Enables read downsampling before downstream analyses (true/false). |
| `kraken2_assign` | Enables taxonomic classification using Kraken2 (true/false). |

---

# Downsampling

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

# Kraken2

| Parameter | Description |
|-----------|-------------|
| `kraken2_db` | Path to the Kraken2 database (PATH). |
| `format_mpa` | Produces MetaPhlAn-compatible output format. If disabled, MPA-related processes should also be removed from the workflow (true/false). |

---

# Merged reads (FLASH)

| Parameter | Description |
|-----------|-------------|
| `min_overlap` | Minimum overlap length required to merge paired-end reads (INT). |
| `max_overlap` | Maximum overlap length allowed for merging (INT). |
| `dovetail_overlap` | Allows merging of reads with dovetail overlaps (true/false). |

---

# BLASTN amplicon classification

| Parameter | Description |
|-----------|-------------|
| `blast_db` | Path to the BLAST nucleotide database (PATH). |
| `perc_id` | Minimum sequence identity (%) required for a strict match (INT). |
| `loose_id` | Minimum sequence identity (%) required for a loose match (INT). |
| `query_cov` | Minimum query coverage (%) per HSP for a strict match (INT). |
| `loose_cov` | Minimum query coverage (%) per HSP for a loose match (INT). |
| `min_qlen` | Minimum query sequence length (bp) for a strict match (INT). |
| `loose_qlen` | Minimum query sequence length (bp) for a loose match (INT). |
| `delta_bitscore` | Maximum allowed bitscore difference between the best hit and alternative hits before reporting an ambiguous assignment (INT). |

---

# Working directory

| Parameter | Description |
|-----------|-------------|
| `workDir` | Nextflow working directory where intermediate files are stored. By default: `${params.result}/work` (PATH). |