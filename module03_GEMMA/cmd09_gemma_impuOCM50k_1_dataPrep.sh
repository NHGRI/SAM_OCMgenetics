#!/bin/bash

#####################################################
## We will run GEMMA GWAS on OCM 50k Region dataset
## First extract the 711 samples
## Generate cov for the datasets
## data files saved: ${DATA_DIR}/oriOCM50k_phenosex/
## gemma results:    ${DATA_DIR}/gemma_output/gemma_results_oriOCM50k/
## gemma GWAS plot:  ${DATA_DIR}/gemma_output/
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#####################################################

# Source configuration
source config_paths.sh

##$$$$$$$$$$$$$$$$$$$
##$$ filtering imputed datasets
##$$$$$$$$$$$$$$$$$$$
module load plink/1.9.0-beta4.4

f_ocm50k=${SAM_OCM50K}

curdd_dir=${ORI_OCM50K_DIR}/
mkdir -p $curdd_dir

outf=${OUTF_OCM50K}


echo '============================'
echo 'Update phenotype and sex variable'
plink \
--bfile ${f_ocm50k} \
--chr 1-22 \
--pheno ${PHENO_FILE} \
--update-sex ${SEX_FILE} \
--make-bed \
--out ${curdd_dir}genoimp_21Sep2023_mafQC005_OCM50k

echo '============================'
echo 'reduce to 711 samples'
plink --bfile ${curdd_dir}genoimp_21Sep2023_mafQC005_OCM50k --remove ${DATA_DIR}/list/list.id.rm  --make-bed --out ${curdd_dir}genoimp_21Sep2023_mafQC005_OCM50k_711


echo '============================'
echo 'generate cov values'
module load R
mkdir -p ${curdd_dir}/gemma_input
module load R
R --no-save < r_scp_Gemma_dataPrep_cov_oriOCM50k.R > rout_Gemma_dataPrep_cov_oriOCM50k.rout


mkdir -p ${GEMMA_OCM50K_DIR}/

f_out_gemma=${GEMMA_OCM50K_DIR}/

f_GRM=${QC_GSTEM}_711
outf=${OUTF_OCM50K}

echo '============================'
echo '  check file locations for GEMMA input files: '
wc -l ${curdd_dir}genoimp_21Sep2023_mafQC005_OCM50k_711.{bim,fam}

wc -l ${curdd_dir}gemma_input/inCov_SEX_C1.txt 

ls -l ${GEMMA_ARRAY_DIR}/${f_GRM}_pruned_standrel.sXX.txt