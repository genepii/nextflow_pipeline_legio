# Assembly + MLST Pipeline (Nextflow)

This pipeline assembles Illumina data and then searches for mutations in a set of genes to characterise the strain under study (Multilocus Sequence Typing approach = MLST).
It supports paired-end data and optional adapter trimming. Different treatment approaches are also supported depending on the strain of Legionella identified. 

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

![Pipeline simplified DAG](../img/assembly-mlst_DAG.png)

---

## Features

* Supports Illumina paired-end data
* Quality control and filtering (QPhred and minimum length)
  * Optional adapters trimming step

* Produces a summary file listing all software used and their parameters (softwaresTrackfile)
* Full execution details available via `--help`

---

## Run

```bash
./start_assembly_mlst.sh -d <seq_id> [options]
```

---

## Required

* `-d, --seq_id` : sequencing run ID (ex: for "Legionella-20260331", run_ID = 20260331)

## Options

#### Input / Output

* `-i, --input` : path to the folder containing the raw data (Illumina sequencing)
* `-o, --output` : folder for results files, relevant to user 
* `-w, --work` : permanent backup folder for all output files, relevant to dev
* `-s, --save` : permanent backup folder for input data
* `-m, --tmp` : temporary backup folder for input data, removed when analysis ended

```
By default : 

/srv/autofs/nfs4/cluqumngs/TMP_IAI/04_CNR_Legionella/
├── NGS_results/
│   └── Genomes/
│       └── <sequencing_id>/
│           └── <YYYYMMDD>_Assembly-MLST/
│               └── results files, for User [--output]
│
├── Raw_fastq/
│   └── Genomes/
│       └── <sequencing_id>/
│           └── all input files [--save] <== not touched, only for data saving
│
/srv/scratch/iai/bachcl/
├── Raw_fastq/
│   └── Legionella/
│       └── Genomes/
│           └── <sequencing_id>/
│               └── all input files [--tmp] <== used during analysis and removed when ended
│
└── result/
    └── Legionella/
        └── Genomes/
            └── <sequencing_id>/
                └── <YYYYMMDD>_Assembly-MLST/
                    ├── work/ <== removed when analysis ended
                    └── all results and created files, for Dev [--work]
```

#### Preprocessing

* `-t, --adapters` : enable adaptaters trimming

  * `True` : removes Illumina TruSeq Adapters (`by default`)
  * `False` : removes only bad quality reads - QPhred and minimum length

> Illumina TruSeq Adapter Read 1 :
> AGATCGGAAGAGCACACGTCTGAACTCCAGTCA 
>
> Illumina TruSeq Adapter Read 2 :
> AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT


#### Analysis mode


#### Help to developpers

* `-c, --config` : path to a new nextflow config file, for developping new parameters

---

