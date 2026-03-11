#!/bin/bash
################################################################################
# Path Configuration File for Genomic Analysis Pipeline
# Source this file before running any analysis scripts: source config_paths.sh
################################################################################

# Base directories for original genotype data
export BASE_DIR="/path/to/base/directory"

# Imputated data directory
export IMPUTATION_DIR="/path/to/imputation/directory"

# Project root directory
export PROJECT_DIR="/path/to/project/directory"

# Data directory (relative path from scripts directory)
export DATA_DIR="../intermediate_data"

# Specific subdirectories
export ARRAY_PHENO_DIR="${DATA_DIR}/array_phenosex"
export IMPU_GENO_DIR="${DATA_DIR}/impuGeno_phenosex"
export ORI_OCM50K_DIR="${DATA_DIR}/oriOCM50k_phenosex"
export MDS_DIR="${DATA_DIR}/mds"
export LIST_DIR="${DATA_DIR}/list"
export GEMMA_OUTPUT_DIR="${DATA_DIR}/gemma_output"

# Specific GEMMA results directories
export GEMMA_ARRAY_DIR="${GEMMA_OUTPUT_DIR}/gemma_results_array"
export GEMMA_IMPU_DIR="${GEMMA_OUTPUT_DIR}/gemma_results_impuGeno"
export GEMMA_OCM50K_DIR="${GEMMA_OUTPUT_DIR}/gemma_results_oriOCM50k"

# Input data files
export PLINK_DIR="${BASE_DIR}/file/name"
export PHENO_FILE="${PLINK_DIR}/file/name"
export SEX_FILE="${BASE_DIR}/file/name"

# Imputation merged file
export IMPU_MERGED="${IMPUTATION_DIR}/path/to/the file"

# OCM50k region file
export SAM_OCM50K="${IMPUTATION_DIR}/path/to/the file"

# Original genotype file stems
export ORI_GSTEM="/original/file/name"
export QC_GSTEM="/afterQC/original/file/name"

# Output file naming conventions
export OUTF_ARRAY="sam_ocm_genoArray"
export OUTF_IMPU="sam_ocm_imputGeno"
export OUTF_OCM50K="sam_ocm_oriOCM50k"

# Create necessary directories
create_directories() {
    mkdir -p ${ARRAY_PHENO_DIR}
    mkdir -p ${IMPU_GENO_DIR}
    mkdir -p ${ORI_OCM50K_DIR}
    mkdir -p ${MDS_DIR}
    mkdir -p ${LIST_DIR}
    mkdir -p ${GEMMA_ARRAY_DIR}
    mkdir -p ${GEMMA_IMPU_DIR}
    mkdir -p ${GEMMA_OCM50K_DIR}
    mkdir -p ${ARRAY_PHENO_DIR}/gemma_input
    mkdir -p ${IMPU_GENO_DIR}/gemma_input
    mkdir -p ${ORI_OCM50K_DIR}/gemma_input
    echo "All directories created successfully"
}

# Validate paths exist
validate_paths() {
    local valid=true
    
    if [ ! -d "${BASE_DIR}" ]; then
        echo "ERROR: BASE_DIR does not exist: ${BASE_DIR}"
        valid=false
    fi
    
    if [ ! -d "${IMPUTATION_DIR}" ]; then
        echo "ERROR: IMPUTATION_DIR does not exist: ${IMPUTATION_DIR}"
        valid=false
    fi
    
    if [ ! -f "${PHENO_FILE}" ]; then
        echo "ERROR: PHENO_FILE does not exist: ${PHENO_FILE}"
        valid=false
    fi
    
    if [ ! -f "${SEX_FILE}" ]; then
        echo "ERROR: SEX_FILE does not exist: ${SEX_FILE}"
        valid=false
    fi
    
    if [ "$valid" = true ]; then
        echo "All critical paths validated successfully"
        return 0
    else
        echo "Path validation failed. Please update config_paths.sh"
        return 1
    fi
}

# Print configuration summary
print_config() {
    echo "================================"
    echo "Path Configuration Summary"
    echo "================================"
    echo "BASE_DIR: ${BASE_DIR}"
    echo "IMPUTATION_DIR: ${IMPUTATION_DIR}"
    echo "PROJECT_DIR: ${PROJECT_DIR}"
    echo "DATA_DIR: ${DATA_DIR}"
    echo "PHENO_FILE: ${PHENO_FILE}"
    echo "SEX_FILE: ${SEX_FILE}"
    echo "================================"
}

# Uncomment to auto-validate when sourced
# validate_paths
# print_config