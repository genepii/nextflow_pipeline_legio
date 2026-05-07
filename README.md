# nextflow_pipeline_legio

A collection of Nextflow pipelines dedicated to the analysis of Legionellosis-related data (CNR des Légionelles).

## Information

To ensure simplicity and independence between tools and workflows, each pipeline has its own dedicated files, all of which include the same tag in their filenames for identification purposes. This is a decision made by the laboratory, and it is therefore possible that certain parts or processes may be similar across multiple pipelines.

- `qiime2_amplicons` for the QIIME2 Amplicons Pipeline
- `blast_amplicons` for the BLASTN Amplicons Pipeline
- `elgato_sbt` for the EL GATO Nested SBT Pipeline

*/!\ The configuration files and values provided here are for information purposes only if the analyses are not run on machines belonging to HCL. The paths and certain options have been configured for those machines, and only the underlying logic and software can be applied in other environments.*


```
project/
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
