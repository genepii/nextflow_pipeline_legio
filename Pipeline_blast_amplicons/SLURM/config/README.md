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

# Process resources

Each process label defines the container image.
There are also labels to specify CPU usage, memory usage and the number of forks permitted. 

---

# Pipeline parameters

The following parameters define the analysis to be performed.

## General parameters

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
| `workDir` | Nextflow working directory where intermediate files are stored (PATH). |