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
| `scratch = true` | Executes each process inside a temporary scratch directory. |
| `QIIME2_DISABLE_PROVENANCE = 1` | Disables QIIME2 provenance generation to reduce output size and execution time. |

---

# Process resources

Each process label defines the container image.
There are also labels to specify CPU usage, memory usage and the number of forks permitted. 

---

# Pipeline parameters

The following parameters define the analysis to be performed.

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
| `all_in_one` | If `true`, processes all samples jointly; otherwise processes samples independently (true/false). |
| `adapters` | Enables adapter trimming using Fastp (true/false). |
| `classifier` | Taxonomic classification method (`sklearn`, `blast`, or `vsearch`). |

---

# FASTP parameters

| Parameter | Description |
|-----------|-------------|
| `min_quality` | Minimum Phred quality score required to retain bases (INT). |
| `min_length` | Minimum read length after trimming (INT). |

---

# DADA2 denoising

| Parameter | Description |
|-----------|-------------|
| `trim_left_f` | Number of bases removed from the beginning of forward reads (INT). |
| `trim_left_r` | Number of bases removed from the beginning of reverse reads (INT). |
| `trunc_len_f` | Forward read truncation length. Used alone for single-end analyses (INT). |
| `trunc_len_r` | Reverse read truncation length (INT). |
| `reads_learn` | Number of reads used to learn the DADA2 error model (INT). |
| `fold_parents` | Minimum abundance fold difference required for parent sequence identification (INT). |
| `min_reads` | Minimum number of reads required for a sample to be retained. Set to `0` when `all_in_one = true`, otherwise `5000`. |

---

# Classifier training

| Parameter | Description |
|-----------|-------------|
| `db` | Name of the reference database (STR). |
| `reads` | Reference FASTA sequences used to train the classifier (PATH). |
| `taxa` | Taxonomic annotation file corresponding to the reference sequences (PATH). |

---

# Taxonomic classification

## sklearn

| Parameter | Description |
|-----------|-------------|
| `sklearn_confidence` | Minimum confidence score required for taxonomic assignment (FLOAT). |

### BLAST

| Parameter | Description |
|-----------|-------------|
| `blast_identity` | Minimum sequence identity (FLOAT). |
| `blast_maxaccepts` | Maximum number of accepted hits (INT). |
| `blast_query_cov` | Minimum query coverage (FLOAT). |

### VSEARCH

| Parameter | Description |
|-----------|-------------|
| `vsearch_identity` | Minimum sequence identity (FLOAT). |
| `vsearch_maxaccepts` | Maximum number of accepted hits (INT). |
| `vsearch_query_cov` | Minimum query coverage (FLOAT). |

---

# Working directory

| Parameter | Description |
|-----------|-------------|
| `workDir` | Nextflow working directory where intermediate files are stored (PATH). |
