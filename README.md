# nextflow_pipeline_legio

A collection of Nextflow pipelines dedicated to the analysis of Legionellosis-related data (CNR des Légionelles).

```
project/
├── assets/             static data (data/db/etc.) for Nextflow pipelines
├── bin/                scripts for modules files
├── config/             config files for Nextflow pipelines
├── launch_sh/          bash scripts for launching Nextflow pipelines
├── modules/            modules files for Nextflow pipelines
├── nextflow_25.10.4    Nextflow executable 
└── workflows_*         Nextflow main workflows
```

---

# Table of contents

- [Nextflow Pipeline Legionelles](#nextflow_pipeline_legio)
  - [QIIME2 Amplicons Pipeline](#qiime2-amplicons-pipeline-nextflow)
    

---

## QIIME2 Amplicons Pipeline (Nextflow)

#### Description
This pipeline analyzes Illumina amplicon sequencing data using QIIME2 and Nextflow.  
It supports single-end and paired-end data, optional adapters trimming, and two analysis modes: all samples together or separately.

---

#### Features

* Supports Illumina single-end and paired-end data
* Quality control and filtering (QPhred and minimum length)
  * Optional adapters trimming step
* Automatic QIIME2 data import and setup
* Taxonomic classification
* **Builds a custom classifier** if not in `assets/qiime_amplicons/`
* Generates FastQC reports and visual outputs (Krona, Barplot)
* Flexible per-sample or global analysis
* Produces a summary file listing all software used and their parameters (softwaresTrackfile)
* Full execution details available via `--help`

---

#### Run

```bash
./launch_sh/start_qiime2_amplicons.sh -d <seq_id> [options]
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

```
By default : 

/srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/
├── NGS_results/
│   └── 23S-5S/
│       └── <sequencing_id>/
│           └── <YYYYMMDD>_Qiime2-amplicons/
│               └── results files, for User [--output]
│
├── Raw_fastq/
│   └── 23S-5S/
│       └── <sequencing_id>/
│           └── all input files [--save] <== not touched, only for data saving
│
/srv/scratch/iai/bachcl/
├── Raw_fastq/
│   └── Legionella/
│       └── 23S-5S/
│           └── <sequencing_id>/
│               └── all input files [--tmp] <== used during analysis and removed when ended
│
└── result/
    └── Legionella/
        └── 23S-5S/
            └── <sequencing_id>/
                └── <YYYYMMDD>_Qiime2-amplicons/
                    ├── work/ <== removed when analysis ended
                    └── all results and created files, for Dev [--work]
```

###### Sequencing mode

* `-pe, --paired` : paired-end data (`True/False`)

###### Preprocessing

* `-t, --adapters` : enable adaptaters trimming (`True/False`)

  * `True` : add adapters trimming during trimming step 
  * `False` : removes only bad quality reads (QPhred and minimum length)

###### Analysis mode

* `-a, --all` : run mode

  * `True` : all samples together
  * `False` : samples processed separately

###### Help to developpers

* `-c, --config` : path to a new nextflow config file, for developping new parameters