# Strain Comparison Pipeline (Nextflow + SLURM)

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

This pipeline performs comparative genomic analyses between strains defined in a metadata file containing sample IDs, sequence types (ST), collection years and origin. It automatically retrieves the corresponding FASTA and FASTQ files, validates the input list, and performs multiple complementary analyses : SNP comparison, Phylogenetic analysis and/or cgMLST analysis.

---

#### Run

```bash
sbatch Strain_Comparison.sbatch -s <Comparison_id> [options]
```

---

#### Required

* `-s` : comparison ID (ex: for metadata file "Strains_20260730.txt", comparison_ID = 20260730)

#### Help to developpers

* `-c` : path to a new nextflow config file, for developping new parameters
