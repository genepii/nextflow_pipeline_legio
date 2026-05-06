# nextflow_pipeline_legio

**[See all the information in the project’s ReadMe; here, only the differences from a standard batch launch are detailed]**

A collection of Nextflow pipelines dedicated to the analysis of Legionellosis-related data (CNR des Légionelles).

*/!\ The configuration files and values provided here are for information purposes only if the analyses are not run on machines belonging to HCL. The paths and certain options have been configured for those machines, and only the underlying logic and software can be applied in other environments.*

```
sbatch_env/
├── config/             config files for Nextflow pipelines
├── launch/             sbatch scripts for launching Nextflow pipelines in SLURM mode
└── dump_params.nf      Nextflow workflow for retrieving write folders from the config file
```

---

- [Nextflow Pipeline Legio](#nextflow_pipeline_legio)
  - [QIIME2 Amplicons Pipeline](#qiime2-amplicons-pipeline-nextflow)
  - [BLASTN Amplicons Pipeline](#blastn-amplicons-pipeline-nextflow)
  - [EL GATO Nested SBT Pipeline](#el-gato-nested-sbt-pipeline-nextflow)


---

## QIIME2 Amplicons Pipeline (Nextflow)

#### Description
This pipeline analyzes Illumina amplicon sequencing data using QIIME2 and Nextflow.  
It supports single-end and paired-end data, optional adapters trimming, and two analysis modes: all samples together or separately; and taxonomic assignment using VSEARCH, BLASTN or a Bayesian classifier (sklearn).

![Pipeline simplified DAG](../img/qiime-amplicons_DAG.png)

---

#### Run

```bash
sbatch start_qiime2_amplicons.sh -s <NGS-WEB_id> [options]
```
---

#### Required

* `-s` : sequencing run ID / NGS-WEB series (ex: for "Legionella-Amplicons-20260331", run_ID = 20260331)

---

<div class="page"/>


## BLASTN Amplicons Pipeline (Nextflow)

#### Description
This pipeline processes Illumina paired-end amplicon sequencing data using Nextflow.
It performs quality control, optional preprocessing steps, and taxonomic identification with a dedicated Legionella workflow.

![Pipeline simplified DAG](../img/blast-amplicons_DAG.png)

---

#### Run

```bash
sbatch start_blast_amplicons.sh -s <NGS-WEB_id> [options]
```

---

#### Required

* `-s` : sequencing run ID / NGS-WEB series (ex: for "Legionella-Amplicons-20260331", run_ID = 20260331)

---

<div class="page"/>

## EL GATO Nested SBT Pipeline (Nextflow)

#### Description
This Nextflow pipeline processes Illumina paired-end Sequencing-Based Typing (SBT) data.
It includes quality control, optional preprocessing steps, taxonomic classification using Kraken2, and MLST profiling (El Gato).

![Pipeline simplified DAG](../img/elgato-nested-sbt_DAG.png)

---

#### Run

```bash
sbatch start_elgato_sbt.sh -s <NGS-WEB_id> [options]
```

---

#### Required

* `-s` : sequencing run ID / NGS-WEB series (ex: for "Legionella-Nested-SBT-20260331", run_ID = 20260331)

---

<div class="page"/>