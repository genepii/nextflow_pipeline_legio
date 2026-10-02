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
| `user` | Current Linux username (`$USER`). If unavailable, defaults to `qiime2-user`. |
| `baseScratch` | Directory used for temporary files and as the container home directory. |
| `matplotlibConfig` | Directory used for Matplotlib’s configuration and cache files. |

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
| `executor.cpus` | `16` | Maximum number of CPUs available to the executor. |
| `executor.memory` | `300 GB` | Maximum amount of memory available to the executor. |

---

# Global process settings

| Parameter | Description |
|-----------|-------------|
| `cpus` | `1` | Default number of CPUs allocated to each process. |
| `maxForks` | `1` | Maximum number of tasks from the same process that can run simultaneously. |
| `memory` | `1GB` | Default amount of memory allocated to each process. |
| `errorStrategy` | `137` | Defines how Nextflow handles a process when an error occurs. |
| `maxRetries` | `2` | Maximum number of times a failed process can be retried. |
| `scratch = true` | Executes each process inside a temporary scratch directory. |
| `QIIME2_DISABLE_PROVENANCE = 1` | Disables QIIME2 provenance generation to reduce output size and execution time. |
| `TMPDIR = baseScratch` | Sets the location for temporary files. |
| `MPLCONFIGDIR = matplotlibConfig` | Sets the location for Matplotlib’s configuration and cache files. |

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
| `analyse_ID` | Analysis identifier automatically generated from the current date (STR). |
| `suffix` | Sequencing run identifier. Usually corresponds to the sequencing date (`YYYYMMDD`) (STR). |
| `input_dir` | Directory containing the input FASTQ files (PATH). |
| `result` | Output directory where all analysis results are written (PATH). |
| `paired_end` | Enables paired-end analysis (`true`) or single-end analysis (`false`) (true/false). |
| `all_in_one` | If `true`, processes all samples jointly; otherwise processes samples independently (true/false). |
| `adapters` | Enables adapter trimming (true/false). |
| `denoising` | If `true`, reads are denoised using DADA2 prior to classification; otherwise reads are merged into sequences and then deduplicated (true/false). |
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

# VSEARCH deduplicating

| Parameter | Description |
|-----------|-------------|
| `trunc_qual` | Quality threshold used for read truncation (INT). |
| `min_overlap` | Minimum required overlap between reads (INT). |
| `max_diffs` | Maximum number of differences allowed between reads (INT). |
| `min_mergelenght` | Minimum length required for merged reads (INT). |

---

# Classifier training

| Parameter | Description |
|-----------|-------------|
| `db` | Name of the reference database (STR). |
| `path_db` | Path to the reference database (PATH). |
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
| `workDir` | Nextflow working directory where intermediate files are stored. By default: `${params.result}/work` (PATH). |