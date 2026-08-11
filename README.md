# nextflow_pipeline_legio

A collection of Nextflow pipelines dedicated to the analysis of Legionellosis-related data (CNR des Légionelles).

## Information

To ensure simplicity and independence between tools and workflows, each pipeline has its own dedicated files, all of which include the same tag in their filenames for identification purposes. This is a decision made by the laboratory, and it is therefore possible that certain parts or processes may be similar across multiple pipelines.

- `qiime2_amplicons` for QIIME2 Amplicons Pipeline
- `blast_amplicons` for BLASTN Amplicons Pipeline
- `elgato_sbt` for EL GATO Nested SBT Pipeline
- `assembly_mlst` for MLST and Assembly Pipeline

*/!\ The configuration files and values provided here are for information purposes only if the analyses are not run on machines belonging to HCL. The paths and certain options have been configured for those machines, and only the underlying logic and software can be applied in other environments.*


```
project/
├── docs/               changelogs and documents
├── img/                images used for ReadMe
├── Pipeline_*/         Nextflow pipelines 
└── nextflow_25.10.4    Nextflow executable
```

And inside each Pipeline_* folders : 

```
Pipeline_*/
├── assets/             static data (data/db/etc.) for Nextflow pipelines
├── bin/                scripts for modules files
├── config/             config files for Nextflow pipelines
├── modules/            modules files for Nextflow pipelines
├── SLURM/              sbatch scripts and config for launching Nextflow pipelines in SLURM mode
├── start_*             bash script for launching Nextflow pipelines in LOCAL mode
└── workflow_*          Nextflow main workflows
```

---

## User instructions for HCL pipeline execution (SLURM / local)

This document clarifies the essential rules for running the pipelines in the HCL context. The order below follows operational priority and execution safety.

---

### 1. SLURM execution requirements (NGS-WEB or direct submission)

When launching in SLURM mode (via NGS-WEB or manual submission), the following rules must be strictly respected:

* The run serial number **must end with 8 digits in the format `YYYYMMDD`**.
* Input data must be stored in the correct project:

  * `Legionella-Amplicons` for Blast or Qiime2 pipelines
  * `Legionella-Nested-SBT` for El Gato pipeline
  * `Legionella` for MLST pipeline
* When launched via NGS-WEB:

  * Only the version of the script stored in the INRIA GitLab repository is executed.
  * The local copy in this repository is provided only for grouping, traceability, and backup purposes.

---

### 2. Configuration-driven execution (no script modification required via NGS-WEB)

All input/output paths and software parameters are fully defined by the configuration file:

* No modification of pipeline scripts is required to adapt an analysis.
* Changing analysis settings is done exclusively via the config file.
* When not using NGS-WEB, a different config file than the default can be provided at submission time using `sbatch` options (see `-h` or pipeline README).
    * In this case, certain paths are specified in the configuration file by the launch script; further details can be found in the linked files.

---

### 3. Local and SLURM input data rules

These rules apply in both execution modes:

* **Input files must be `fastq` or `fastq.gz`**.
* **Data must contain `_R1` and `_R2` in filenames**.
* Only the part before the first `_` is used for downstream naming (= sample_id).

---

### 4. Configuration file usage

* Separate configuration files are used for:

  * LOCAL execution
  * SLURM execution
* It is essential to use the correct configuration format depending on the execution mode.
* File hierarchy and internal descriptions clearly indicate ownership and usage context.

---

### 5. Documentation and user support

* A help command is available via `-h` for all pipelines.
* README files are provided for both users and developers.
* These resources should be consulted before modifying or extending pipeline usage.

---

### 6. Data backup and storage policy

* Input and output data are backed up across different servers.
* This mechanism was initially introduced for development needs.
* It may evolve or be simplified in future versions of the pipeline.
