# Changelog

All notable changes to the pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added 

- Retrieve an old metadata file to merge it with the new one and save it to the local database.
    - All files will contain at least the ID and ST columns.
    - Only one metadata file per strain (Lp, Lb, etc.).
    - Adding allele results for genes of interest.
- Addition of labels for memory and CPU resources.
- Addition of option to change the kmers value used by the ElGato software.
- Checking that FastQ files are not empty at the start.
    - The `min_raw_seq` parameter (config file) specifies the minimum number of sequences required.
- Added an option to prevent new cgMLSTs from being saved to the ReporTree database. 
    - Also save the database and the NWK trees to a folder specified in the start_*.sh script for manual use, and avoid making any changes to the database used by the pipeline.

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
    - MLST_ELGATO module is ignored if an error occurs, to address a specific case of sample contamination where the new version crashes compared to version 1.22.
- Reorganisation of the HTML report.
    - Addition of a button to copy the contents of the tables.
    - Visual changes.
    - Group all strains of L. pneumophila under the name L. pneumophila.
    - Depth highlighted < 20X, not 30X as previously stated.
    - Highlighting ‘nan’ values when the El Gato method fails for L. pneumophila samples.

### Fixed

- input validation logic.
- SLURM execution issues relating to time limits and return codes.
- TMP_DIR and WORK_DIR during execution via SLURM.
- creation of the TMP_DIR directory when launching via bash.
- the absence of a sample in the metadata causes the script to terminate.

### Known issues

- ReporTree issue if IDs consist solely of numeric characters and the leading zeros are omitted.
    - Awaiting a fix from the developers; issue raised on GitHub.
- Specific error when Neu and NeuAH are detected in the same sample. El Gato stops abruptly, so the module is ignored and ‘nan’ is entered into the final table (shown in red in the HTML). Only in El Gato version > 1.22.

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
