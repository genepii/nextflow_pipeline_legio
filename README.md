# nextflow_pipeline_legio
A collection of Nextflow pipelines dedicated to the analysis of Legionellosis-related data (CNR des Légionelles).


---

- [Nextflow Pipeline Legionelles](#nextflow_pipeline_legio)
  - [QIIME2 Pipeline](#qiime2-pipeline-nextflow)


---

## QIIME2 Pipeline (Nextflow)

#### Description
This pipeline analyzes Illumina amplicon sequencing data using QIIME2 and Nextflow.  
It supports single-end and paired-end data, optional trimming, and two analysis modes: all samples together or separately.

---

#### Features

* Supports Illumina single-end and paired-end data
* Quality control and filtering
* Optional trimming step
* Automatic QIIME2 data import and setup
* Taxonomic classification
* **Builds a custom classifier during execution**
* Generates QC reports and visual outputs
* Flexible per-sample or global analysis
* Produces a summary file listing all software used and their parameters
* Full execution details available via `--help`

---

#### Run

```bash
./launch_sh/start_qiime2.sh -d <seq_id> [options]
```

---

#### Required

* `-d, --seq_id` : sequencing run ID (ex: for 
"Legionella-Amplicons-20260331", run_ID = 20260331)

#### Options

###### Input / Output

* `-i, --input` : path to the folder containing the raw data (Illumina sequencing)
* `-o, --output` : folder for results files, relevant to user 
* `-w, --work` : permanent backup folder for all output files, relevant to dev
* `-s, --save` : permanent backup folder for input data
* `-m, --tmp` : temporary backup folder for input data, removed when analysis ended

###### Sequencing mode

* `-pe, --paired` : paired-end data (`True/False`)

###### Preprocessing

* `-t, --trim` : enable trimming (`True/False`)

  * removes adapters
  * filters low-quality reads (Phred < 30)

###### Analysis mode

* `-a, --all` : run mode

  * `True` : all samples together
  * `False` : samples processed separately

