# QIIME2 Amplicons Pipeline (Nextflow + SLURM)

[See all the information in the project’s ReadMe; here, only the differences from a standard batch launch are detailed]

*/!\ The configuration files and values provided here are for information purposes only if the analyses are not run on machines belonging to HCL. The paths and certain options have been configured for those machines, and only the underlying logic and software can be applied in other environments.*

**To run the pipeline via NGS-WEB, the nextflow_pipeline_legio repository must be located in the `/srv/scratch/DEPOT/pipelines` directory**
**and start_\*.sbatch must be placed in the NGS-WEB repository.**

```
SLURM/
├── config/             config files for Nextflow pipeline
├── dump_params.nf      Nextflow workflow for retrieving write folders from the config file
└── start_*             sbatch script for launching Nextflow pipelines in SLURM mode
```

---

This pipeline analyzes Illumina amplicon sequencing data using QIIME2 and Nextflow.  
It supports single-end and paired-end data, optional adapters trimming, and two analysis modes: all samples together or separately; and taxonomic assignment using VSEARCH, BLASTN or a Bayesian classifier (sklearn).

---

#### Run

```bash
sbatch start_qiime2_amplicons.sbatch -s <NGS-WEB_id> [options]
```
---

#### Required

* `-s` : sequencing run ID / NGS-WEB series (ex: for "Legionella-Amplicons-20260331", run_ID = 20260331)

