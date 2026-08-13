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

| Process   | Container                   | CPUs | Max concurrent tasks | Initial memory | Max retries |
| --------- | --------------------------- | ---: | -------------------: | -------------: | ----------: |
| Fastp     | `fastp_v1.3.2.sif`          |    5 |                    3 |           1 GB |           4 |
| Python    | `python_plot_tree_3.11.sif` |    2 |                    8 |           1 GB |           5 |
| Snippy    | `snippy_v4.6.0.sif`         |    4 |                    4 |           1 GB |           5 |
| IQ-TREE   | `iqtree_v3.1.3.sif`         |    6 |                    2 |           1 GB |           5 |
| Gubbins   | `gubbins_v3.4.3.sif`        |    2 |                    2 |           1 GB |           5 |
| Grapetree | `grapetree_2.1.sif`         |    2 |                    3 |           1 GB |           5 |
| ReporTree | `reportree_2.6.1.sif`       |    8 |                    2 |          40 GB |           2 |

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
| `suffix` | Comparison run identifier. Used as a subdirectory name to locate input Metadata file and organize results (STR). |
| `output_dir_prefix` | Base directory where final analysis outputs are stored. The final output directory is automatically built as `${output_dir_prefix}/${suffix}/${analyse_id}_Strain-Comparison` (PATH). |
| `tmp_dir_prefix` | Base directory used for temporary files generated during execution. The final temporary directory is automatically built as `${tmp_dir_prefix}/${suffix}` (PATH). |
| `work_dir_prefix` | Base directory used for Nextflow work files and intermediate results. The final work directory is automatically built as `${work_dir_prefix}/${suffix}/${analyse_id}_Strain-Comparison` (PATH). |
| `metadata_user_prefix` | Base directory used to locate the metadata file defining the strains to be compared (PATH). The final path is automatically built as `${metadata_user_prefix}${suffix}` (PATH). |
| `adapters` | Enables adapter trimming using Fastp (true/false). |
| `zoom` | ReporTree zoom mode: analyze samples from metadata, all samples, no zoom, or selected sample IDs (`analyse`, `none`, or one/multiple sample IDs such as `SampleA` or `SampleA,SampleB`). |

---

# FASTA research parameters

| Parameter      | Description                                                                              |
| -------------- | ---------------------------------------------------------------------------------------- |
| `path_fasta`   | List of directories searched for FASTA files (LIST).                                     |
| `fasta_prefix` | Prefix pattern used to filter FASTA files. `false` disables prefix filtering (BOOL/STR). |
| `fasta_suffix` | Suffix pattern used to filter FASTA files. `false` disables suffix filtering (BOOL/STR). |
| `fasta_ext`    | List of accepted FASTA file extensions (LIST).                                           |

---

# FASTQ research parameters

| Parameter      | Description                                                                                                          |
| -------------- | -------------------------------------------------------------------------------------------------------------------- |
| `path_fastq`   | List of directories searched for FASTQ files (LIST).                                                                 |
| `fastq_prefix` | Prefix pattern used to filter FASTQ files. `false` disables prefix filtering (BOOL/STR).                             |
| `fastq_suffix` | List of patterns used to identify paired-end FASTQ files, including `_R1/_R2` and `_1/_2` naming conventions (LIST). |
| `fastq_ext`    | List of accepted FASTQ file extensions (LIST).                                                                       |

---

# FastP parameters

| Parameter | Description |
|-----------|-------------|
| `min_quality` | Minimum Phred quality score required to retain bases (INT). |
| `min_length` | Minimum read length after trimming (INT). |

---

# Snippy parameters

| Parameter     | Description                                                                                                                                 |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `snippy_dir` | Root directory containing Snippy databases and reference genomes (STR).              |
| `snippy_ref`  | List of sequence-type-specific Snippy configurations, including the ST identifier, previous results directory, and reference genome (LIST). |
| `snippy_opts` | Additional command-line options passed to Snippy (STR).                                                                                     |

---

# IQ-TREE paramaters 

| Parameter          | Description                                                                                |
| ------------------ | ------------------------------------------------------------------------------------------ |
| `iqtree_bootstrap` | Number of standard nonparametric bootstrap replicates used to assess branch support (INT). |
| `iqtree_alrt`      | Number of SH-aLRT replicates used to assess branch support (INT).                          |
| `iqtree_model`     | Evolutionary substitution model used for phylogenetic inference (STR).                     |
| `iqtree_color_by`  | Metadata column used to assign colors to samples in the phylogenetic tree (STR).           |

---

# Gubbins parameters

| Parameter         | Description                                                            |
| ----------------- | ---------------------------------------------------------------------- |
| `gubbins_percent` | Percentage threshold used by Gubbins for recombination analysis (INT). |

---

# chewBBACA parameters

| Parameter       | Description                                                                                                                                                               |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `chewbbaca_set` | List of chewBBACA analysis configurations, including the number of loci, gene schema, previous results, nomenclature, and Grapetree/ReporTree activation settings (LIST). |
| `st`            | Sequence type associated with the chewBBACA configuration (INT).                                                                                                          |
| `nb`            | Number of loci in the cgMLST schema (INT).                                                                                                                                |
| `genes`         | Path to the chewBBACA gene schema (STR).                                                                                                                                  |
| `previous`      | Path to previous chewBBACA or ReporTree results used as input/reference data (STR).                                                                                       |
| `nomenclature`  | Path to the nomenclature/partition file, or `none` when no nomenclature file is used (STR).                                                                               |
| `grapetree`     | Enables or disables Grapetree analysis for the corresponding chewBBACA dataset (BOOL).                                                                                    |
| `reportree`     | Enables or disables ReporTree analysis for the corresponding chewBBACA dataset (BOOL).                                                                                    |

---

# Grapetree parameters

| Parameter     | Description                                    |
| ------------- | ---------------------------------------------- |
| `grape_model` | Clustering/tree model used by Grapetree (STR). |


---

# ReporTree parameters 

| Parameter            | Description                                                                                                        |
| -------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `rep_interest`       | Threshold or cluster identifier(s) of interest for ReporTree analysis. Can be a single value or a list (INT/LIST). |
| `rep_site_inclusion` | Minimum proportion of included sites required for the analysis (FLOAT).                                            |
| `rep_min_allele`     | Minimum number of allelic differences used to define the cluster level. `none` considers all clusters (INT/STR).   |
| `rep_max_allele`     | Maximum number of allelic differences considered. `none` uses the cluster defined by `rep_min_allele` (INT/STR).   |
| `rep_loci_called`    | Minimum proportion of loci that must be called for a sample to be included (FLOAT).                                |
| `rep_col_metadata`   | Metadata column used by ReporTree for sample annotation/grouping (STR).                                            |
| `rep_model`          | Clustering model used by ReporTree (STR).                                                                          |
| `rep_analysis`       | Analysis method used by ReporTree (STR).                                                                           |

---

# Working directory

| Parameter | Description |
|-----------|-------------|
| `workDir` | Nextflow working directory containing intermediate files (PATH). |
