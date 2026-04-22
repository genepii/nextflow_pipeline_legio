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

- [TMP_HCL](#tmp_hcl)
  - [QIIME2 Amplicons Pipeline](#qiime2-amplicons-pipeline-nextflow)
  - [BLASTN Amplicons Pipeline](#blastn-amplicons-pipeline-nextflow)


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


---

<div class="page"/>


## BLASTN Amplicons Pipeline (Nextflow)

#### Description
This pipeline processes Illumina paired-end amplicon sequencing data using Nextflow.
It performs quality control, optional preprocessing steps, and taxonomic identification with a dedicated Legionella workflow.

---

#### Features

* Paired-end Illumina data only (single-end not supported)
* Quality control and filtering (QPhred and minimum length)
  * FastQC reports generated at each processing step
  * Detailed statistics files for traceability
* Optional preprocessing steps:
  * Adapter trimming
  * Decontamination against a human genome database
  * Downsampling to a user-defined fraction of reads
* Optional taxonomic assignment using Kraken2
* Dereplication and merging of FASTQ files into FASTA format
* **Sequence identification using BLASTN against a Legionella-specific database**
  * Two parameter sets are applied:
    * Strict mode: highly stringent thresholds to ensure high-confidence identifications
    * Loose mode: relaxed thresholds to explore borderline matches and assess whether less stringent parameters recover additional relevant signals
* Generation of visual outputs (e.g., barplots or Krona)
* Comprehensive tracking of all tools, versions, and parameters used (softwaresTrackfile)
* Full pipeline usage and parameters accessible via `--help`

---

#### Run

```bash
./launch_sh/start_blast_amplicons.sh -d <seq_id> [options]
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
│           └── <YYYYMMDD>_Blast-amplicons/
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
                └── <YYYYMMDD>_Blast-amplicons/
                    ├── work/ <== removed when analysis ended
                    └── all results and created files, for Dev [--work]
```

###### Preprocessing

* `-a, --adapters` : enable adaptaters trimming (`True/False`)

  * `True` : add adapters trimming during trimming step 
  * `False` : removes only bad quality reads (QPhred and minimum length)

* `-e, --deconta` : enable decontamination of reads (`True/False`)

  * `True` : keep the reads that are not aligned to the human database 
  * `False` : keep all the reads given as input

* `-n, --down` : enable downsampling (`float`)

  * `float` : percentage of reads retained for analysis
  * if the option is not used, all the reads are kept for analysis

###### Optionnal analysis

* `-k, --kraken` : enable classification with Kraken2

  * `True` : perform Kraken2 analysis
  * `False` : skip Kraken2

###### Help to developpers

* `-c, --config` : path to a new nextflow config file, for developping new parameters
