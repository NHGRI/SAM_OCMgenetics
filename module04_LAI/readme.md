# Genomic Data Preparation and Local Ancestry Inference Pipeline

This repository contains a two-step pipeline for preparing genotype data, merging study and reference panels, and performing **local ancestry inference (LAI)** followed by **Tractor association analysis**.

All paths are replaced with symbolic placeholders (e.g. `<<STUDY_BFILE>>`) - Configuration can be handled via environment variables or a local config file

------------------------------------------------------------------------

## Overview of the Pipeline

------------------------------------------------------------------------

## Script 1: `01_dataPrep_qc_merge_subset.sh`

### Purpose

This script performs **data preparation and harmonization** prior to local ancestry inference. It ensures that study and reference datasets are compatible, high-quality, and split into chromosome-level inputs required by downstream tools.

### Key Steps

1.  **Quality control (QC)**
    -   MAF and genotype missingness filtering
    -   Autosomes only
2.  **Variant harmonization**
    -   Convert variant IDs to `CHR:POS`
    -   Remove duplicate variants
3.  **Merge study and reference panels**
    -   Identify overlapping SNPs
    -   Handle strand flips
    -   Exclude unresolved mismatches
4.  **Post-merge QC**
    -   Apply additional MAF / GENO filters
    -   Remove strand-ambiguous A/T and G/C SNPs
5.  **Dataset splitting**
    -   Separate merged data into:
        -   study-only samples
        -   reference-only samples
6.  **Region subsetting and chromosome split**
    -   Extract target genomic regions (e.g. candidate loci)
    -   Generate per-chromosome PLINK files (chr1–chr22)

### Inputs (symbolic)

-   Study PLINK dataset: `<<STUDY_BFILE>>`
-   Reference PLINK dataset: `<<REF_BFILE>>`
-   Optional phenotype and sex update files
-   Sample keep lists for study and reference
-   Genomic region range file

### Outputs

-   QC’d and merged PLINK datasets
-   Study-only and reference-only subsets
-   Per-chromosome PLINK files ready for phasing and LAI

------------------------------------------------------------------------

## Script 2: `02_LAI_rfmix_tractor.sh`

### Purpose

This script performs **local ancestry inference and ancestry-aware association testing** using SHAPEIT, RFMix, and Tractor.

It assumes Script 1 has already been run.

### Key Steps

1.  **Phasing (SHAPEIT)**
    -   Phase study and reference samples separately
    -   Per chromosome (1–22)
    -   Uses genetic maps
2.  **VCF conversion**
    -   Convert phased haplotypes to VCF format
3.  **Local ancestry inference (RFMix)**
    -   Infer ancestry tracts using phased study and reference data
    -   Outputs local ancestry calls per individual and chromosome
4.  **Tractor analysis**
    -   Extract ancestry-specific haplotype dosages
    -   Run ancestry-aware regression (e.g. logistic or linear)
    -   Generate per-chromosome Tractor association results

### Inputs (symbolic)

-   Per-chromosome study PLINK datasets
-   Per-chromosome reference PLINK datasets
-   Genetic maps for phasing
-   RFMix sample-to-population map
-   Tractor phenotype file

### Outputs

-   RFMix local ancestry files
-   Tractor ancestry-specific association result files
-   One output file per chromosome

------------------------------------------------------------------------

## Software Requirements

The pipeline assumes the following tools are available (via modules or PATH):

-   **PLINK 1.9**
-   **PLINK 2.0**
-   **SHAPEIT v2**
-   **RFMix**
-   **Python (≥3.6)** for Tractor
-   **Tractor** (GitHub: <https://github.com/Atkinson-Lab/Tractor>)

Exact versions can be adapted to your computing environment.

------------------------------------------------------------------------

## Configuration and GitHub Safety

All scripts use **symbolic placeholders** 
