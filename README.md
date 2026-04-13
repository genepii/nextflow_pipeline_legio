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

#### Run

```bash
./launch_sh/start_qiime2.sh -d <seq_id> [options]
```

---

#### Required

* `-d, --seq_id` : sequencing run ID

#### Options

###### Input / Output

* `-i, --input` : input folder (raw data)
* `-o, --output` : output folder
* `-s, --save` : save folder
* `-m, --tmp` : temporary folder
* `-w, --work` : Nextflow work directory

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

---

#### Features

* Supports Illumina single-end and paired-end data
* Quality control and filtering
* Optional trimming step
* QIIME2 import and data setup
* Taxonomic classification
* QC reports and visual outputs
* Flexible per-sample or global analysis

