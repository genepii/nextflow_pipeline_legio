# Assembly + MLST Pipeline (Nextflow + SLURM)

[See all the information in the project’s ReadMe; here, only the differences from a standard batch launch are detailed]

*/!\ The configuration files and values provided here are for information purposes only if the analyses are not run on machines belonging to HCL. The paths and certain options have been configured for those machines, and only the underlying logic and software can be applied in other environments.*

**To run the pipeline via NGS-WEB, the nextflow_pipeline_legio repository must be located in the `/srv/scratch/DEPOT/pipelines` directory**
**and {script}.sbatch must be placed in the NGS-WEB repository.**

```
SLURM/
├── config/             config files for Nextflow pipeline
├── dump_params.nf      Nextflow workflow for retrieving write folders from the config file
└── {script}.sbatch     sbatch script for launching Nextflow pipelines in SLURM mode
```

---

This pipeline assembles Illumina data and then searches for mutations in a set of genes to characterise the strain under study (Multilocus Sequence Typing approach = MLST).
It supports paired-end data and optional adapter trimming. Different treatment approaches are also supported depending on the strain of Legionella identified. 

---

#### Run

```bash
sbatch MLST-assembly.sbatch -s <NGS-WEB_id> [options]
```

---

#### Required

* `-s` : NGS-WEB series (ex: "Legionella-20260331")

#### Help to developpers

* `-c` : path to a new nextflow config file, for developping new parameters
