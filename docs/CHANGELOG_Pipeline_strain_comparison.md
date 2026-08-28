# Changelog

All notable changes to the pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added 

- Generation of an HTML report for each comparison.
    - Addition of a table of contents to the HTML report.
    - Addition of sections dedicated to:
        - Summary of the analysis steps carried out;
        - Data analysed;
        - Detection of SNPs;
        - Phylogenetic trees generated.
    - Improved presentation of HTML headings and anchors.
- Addition of labels for memory and CPU resources.

### Changed

- Adaptation of the pipeline for execution across multiple machines and users.
    - Adjustment of DB path.
    - Change permissions for folders created by the launch script.

### Fixed

- Fixed SLURM execution issues relating to time limits and return codes.
- Fixed issues relating to TMP_DIR and WORK_DIR during execution via SLURM.
- Fixed issue relating to the creation of the TMP_DIR directory when launching via bash.

### Known issues

- ReporTree issue if IDs consist solely of numeric characters and the leading zeros are omitted.
    - Awaiting a fix from the developers; issue raised on GitHub.

---

## [Alpha] - 2026-08-12

### Added 

- Implementation of a strain comparison pipeline.
- Addition of a process to search for FASTA assemblies / FASTQ files associated with samples.
    - Recursive search for available assemblies/reads.
    - Selection of the most recent assembly/reads where multiple versions exist.
- Enrichment of the metadata table with the path to the FASTA assembly / FASTQ files.
