# Genomic Analysis Pipeline - GEMMA

## Overview

This repository contains bash scripts for genomic variant association analysis. All specific file paths have been replaced with configurable environment variables for portability and security.

## Setup Instructions

### 1. Configure Paths

Edit `config_paths.sh` to set your environment-specific paths:

``` bash
export BASE_DIR="/path/to/base/directory"
export IMPUTATION_DIR="/path/to/imputation/directory"
export PROJECT_DIR="/path/to/project/directory"
```

### 2. Source Configuration

Before running any analysis script, source the configuration file:

``` bash
source config_paths.sh
```

### 3. Create Directory Structure

Create all necessary directories:

``` bash
source config_paths.sh
create_directories
```

## Path Mapping Reference

| Variable         | Description                               | Original Path |
|-------------------|-------------------------|----------------------------|
| `BASE_DIR`       | Base directory for original genotype data | NA            |
| `IMPUTATION_DIR` | Imputated data directory                  | NA            |
| `PROJECT_DIR`    | Project root directory                    | NA            |
| `DATA_DIR`       | Data directory (relative)                 | NA            |
| `PHENO_FILE`     | Phenotype file                            | NA            |
| `SEX_FILE`       | Sex update file                           | NA            |

## Analysis Pipeline

### MDS Analysis (cmd06)

``` bash
source config_paths.sh
bash cmd06_mds.sh
bash cmd06_2_mds_ocm50k.sh
```

-   Calculates MDS eigenvalues for pooled and population-specific data
-   Generates IBD estimates
-   Performs population stratification analysis

### GEMMA GWAS - Array Data (cmd07)

``` bash
source config_paths.sh
bash cmd07_gemma_array_1_GRM.sh      # Generate GRM
bash cmd07_gemma_array_2_maf_run.sh  # Run GEMMA
bash cmd07_gemma_array_3_plot.sh     # Plot results
```

-   Filters to 711 samples with valid phenotypes
-   Generates genetic relationship matrices
-   Runs GWAS with MAF thresholds (0.01, 0.05)
-   Adjusts for sex and principal components
-   Generates Manhattan and QQ plots

### GEMMA Association - OCM50k Region (cmd09)

``` bash
source config_paths.sh
bash cmd09_gemma_impuOCM50k_1_dataPrep.sh  # Prepare data
bash cmd09_gemma_impuOCM50k_2_run.sh       # Run GEMMA
```

-   Focuses on OCM 50kb region

## R Scripts

### one line summary

-   **r_scp_Gemma_dataPrep_cov_array.R**: Builds GEMMA covariate files (intercept + sex + selected MDS PCs) for the array genotype dataset.
-   **r_scp_Gemma_dataPrep_cov_oriOCM50k.R**: Builds GEMMA covariate files (intercept + sex + selected MDS PCs) for the original OCM50k dataset.
-   **r_scp_plot_GEMMA_gwas_arg.R**: Generates QQ and Manhattan plots from a GEMMA GWAS association output specified via command-line arguments.
-   **r_proc_plot_GEMMA_OCM50k.R**: Produces QQ and Manhattan plots for GEMMA results in the OCM50k region, including an alternative labeled Manhattan plot.

### Configuration (paths) for R script

All cleaned scripts remove hard-coded absolute paths and instead use environment variables:

-   `SAMOCM_PROJECT_DIR`: directory where these scripts live (default: current working directory).
-   `SAMOCM_Analyst_DATA_DIR`: path to `Analyst``_data` (default: `file.path(SAMOCM_PROJECT_DIR, '..', '``Analyst``_data')`).

## Key Output Files

### Genetic Relationship Matrices (GRM)

-   `${GEMMA_ARRAY_DIR}/*standrel.sXX.txt`

### GWAS Results

-   `${GEMMA_ARRAY_DIR}/gemma_out_*.assoc.txt` - Array-based GWAS
-   `${GEMMA_IMPU_DIR}/gemma_out_*.assoc.txt` - Imputation-based GWAS
-   `${GEMMA_OCM50K_DIR}/gemma_out_*.assoc.txt` - OCM50k region GWAS

### MDS/PCA Results

-   `${MDS_DIR}/*_mds.mds` - MDS coordinates

### Sample Lists

-   `${LIST_DIR}/list.id.rm` - Samples to remove (n=19)
-   `${LIST_DIR}/list.id.JAM` - Jamaican samples
-   `${LIST_DIR}/list.id.MAL` - Malawian samples

## Module Requirements

``` bash
module load plink/1.9.0-beta4.4
module load GEMMA/0.96
module load R
```

## Analysis Parameters

### MAF Thresholds

-   0.01 (1% minor allele frequency)
-   0.05 (5% minor allele frequency)

### Covariate Models

-   `SEX` - Sex only
-   `SEX_C1` - Sex + MDS C1
-   `SEX_C1-C2` - Sex + MDS C1-2
-   `SEX_C1-C3` - Sex + MDS C1-3

### LD Pruning

-   Window: 50 SNPs
-   Step: 10 SNPs
-   r² threshold: 0.1

## Sample Sizes

-   Final analysis: 711 samples
    -   Jamaican (JAM): 340 samples
    -   Malawian (MAL): 371 samples

## Citation

If using these scripts, please cite the original study and acknowledge the data sources.

## Contact

For questions about the analysis pipeline, please contact the study investigators.


