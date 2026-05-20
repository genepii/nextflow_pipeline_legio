# EL GATO Nested SBT Pipeline (Nextflow + SLURM)

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

This Nextflow pipeline processes Illumina paired-end Sequencing-Based Typing (SBT) data.
It includes quality control, optional preprocessing steps, taxonomic classification using Kraken2, and MLST profiling (El Gato).


---

#### Run

```bash
sbatch start_elgato_sbt.sbatch -s <NGS-WEB_id> [options]
```

---

#### Required

* `-s` : NGS-WEB series (ex: "Legionella-Nested-SBT-20260331")

#### Help to developpers

* `-c` : path to a new nextflow config file, for developping new parameters
