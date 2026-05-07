# BLASTN Amplicons Pipeline (Nextflow + SLURM)

**[See all the information in the project’s ReadMe; here, only the differences from a standard batch launch are detailed]**

A collection of Nextflow pipelines dedicated to the analysis of Legionellosis-related data (CNR des Légionelles).

*/!\ The configuration files and values provided here are for information purposes only if the analyses are not run on machines belonging to HCL. The paths and certain options have been configured for those machines, and only the underlying logic and software can be applied in other environments.*

```
SLURM/
├── config/             config files for Nextflow pipeline
├── dump_params.nf      Nextflow workflow for retrieving write folders from the config file
└── start_*             sbatch script for launching Nextflow pipelines in SLURM mode
```

---

This pipeline processes Illumina paired-end amplicon sequencing data using Nextflow.
It performs quality control, optional preprocessing steps, and taxonomic identification with a dedicated Legionella workflow.


---

#### Run

```bash
sbatch start_blast_amplicons.sh -s <NGS-WEB_id> [options]
```

---

#### Required

* `-s` : sequencing run ID / NGS-WEB series (ex: for "Legionella-Amplicons-20260331", run_ID = 20260331)
