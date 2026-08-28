# Changelog

All notable changes to the pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added 
- Addition of labels for memory and CPU resources.

### Changed

- Adaptation of the pipeline for execution across multiple machines and users.
    - Adjustment of DB path.
    - Change permissions for folders created by the launch script.

### Fixed

- Fixed SLURM execution issues relating to time limits and return codes.
- Fixed issues relating to TMP_DIR and WORK_DIR during execution via SLURM.
- Fixed issue relating to the creation of the TMP_DIR directory when launching via bash.

---

## [Beta] - 2026-05-06

### Added 

- Implementation of a BLAST pipeline dedicated to amplicons.
- Management of amplicon sequences as inputs to the pipeline.
- Automation of similarity searches against reference databases.
- Generation of structured BLAST results for the interpretation of amplicons.
- Organisation of results by sample/amplicon.
