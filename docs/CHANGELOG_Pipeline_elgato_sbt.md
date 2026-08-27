# Changelog

All notable changes to FastANI database will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added 
- Addition of labels for memory and CPU resources.

### Changed

- Adaptation of the pipeline for execution across multiple machines.
    - Adjustment of DB path.
- El Gato version update from 1.22.0 to 1.23.0.

### Fixed

- Fixed SLURM execution issues relating to time limits and return codes.
- Fixed issues relating to TMP_DIR and WORK_DIR during execution via SLURM.
- Fixed issue relating to the creation of the TMP_DIR directory when launching via bash.

---

## [Beta] - 2026-05-06

### Added 

- Implementation of the Elgato SBT analysis pipeline.
- Automation of the preparation of sequences required for typing.
- Integration of the typing and result generation stages.
- Structuring of results to enable their use in downstream reports.
