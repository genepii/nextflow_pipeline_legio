# Changelog

All notable changes to FastANI database will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added 

- Retrieve an old metadata file to merge it with the new one and save it to the local database.
    - All files will contain at least the ID and ST columns.
    - Only one metadata file per strain.
- Addition of labels for memory and CPU resources.

### Changed

- Reorganisation of the modules and the main workflow.
    - Renamed processes.
- Improvements to the management of input files and validation checks.
    - The ‘Linked_to’ column has been made optional.
- Adaptation of the pipeline for execution across multiple machines and users.
    - Adjustment of DB path.
    - Change permissions for folders created by the launch script.
- Improvements to the management and publication of output.
    - Directory structure.
    - Saving MLST results for all schemas, not just the one used in ReporTree.
- El Gato version update from 1.22.0 to 1.23.0.
- Reorganisation of the HTML report.
    - Addition of a button to copy the contents of the tables.
    - Visual changes.

### Fixed

- Fixed the input validation logic.
- Fixed SLURM execution issues relating to time limits and return codes.
- Fixed issues relating to TMP_DIR and WORK_DIR during execution via SLURM.
- Fixed issue relating to the creation of the TMP_DIR directory when launching via bash.

### Known issues

- ReporTree issue if IDs consist solely of numeric characters and the leading zeros are omitted.
    - Awaiting a fix from the developers; issue raised on GitHub.

---

## [Beta] - 2026-07-08

### Added 

- Addition and structuring of the MLST assembly and typing pipeline.
    - Integration of genome assembly from reads.
    - Integration of MLST analysis.
    - Integration of MLST result clustering according to different schemes.
- Searching for and managing input FASTQ files, with support for:
    - *.fastq
    - *.fastq.gz
- Searching for and managing input metadata file.
- Incorporation of different analysis paths depending on the results obtained during the analysis.
- Generation and centralisation of assembly and MLST results.
    - Generation of an HTML report highlighting outliers or AMRs.
        - Addition of sections dedicated to:
            - quality indicators;
            - column descriptions;
            - Assembly & MLST summary;
            - AMR research.
        - Dynamic generation of report sections based on the strains present in the data.
    - Production of summary files suitable for downstream analyses.
    - Publication of results in a dedicated directory structure.
- Saving new MLST/clustering results to the local database.
