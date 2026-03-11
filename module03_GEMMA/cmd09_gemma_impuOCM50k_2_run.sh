#!/bin/bash

#####################################################
## We will run GEMMA GWAS on natasha's OCM_Region dataset
## First extract the 711 samples
## Generate cov for the datasets
## data files saved: ${DATA_DIR}/oriOCM50k_phenosex/
## gemma results:    ${DATA_DIR}/gemma_output/gemma_results_oriOCM50k/
## gemma GWAS plot:  ${DATA_DIR}/gemma_output/
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#####################################################

# Source configuration
source config_paths.sh

##########################
## run GEMMA with SEXs and PCs
##########################
module load GEMMA/0.96
mkdir -p ${GEMMA_OCM50K_DIR}/

f_out_gemma=${GEMMA_OCM50K_DIR}/
curdd_dir=${ORI_OCM50K_DIR}/
f_GRM=${QC_GSTEM}_711
outf=${OUTF_OCM50K}

echo '============================'
echo '  check file locations for GEMMA input files: '
wc -l ${curdd_dir}genoimp_21Sep2023_mafQC005_OCM50k_711.{bim,fam}

wc -l ${curdd_dir}gemma_input/inCov_SEX_C1.txt 

ls -l ${GEMMA_ARRAY_DIR}/${f_GRM}_pruned_standrel.sXX.txt 

echo '============================'
echo 'Run GEMMA '
## Natasha's dataset is based of 005 cutoff
for maf in 005
do
  for cov in SEX SEX_C1 SEX_C1-C2 SEX_C1-C3
  do
  gemma -bfile ${curdd_dir}genoimp_21Sep2023_mafQC005_OCM50k_711 \
     -k ${GEMMA_ARRAY_DIR}/${f_GRM}_pruned_standrel.sXX.txt \
     -lmm 4 \
  	 -c ${curdd_dir}gemma_input/inCov_${cov}.txt \
  	 -o gemma_out_${outf}_mafQC${maf}_cov${cov}
  done
done

# ## move all the results and log to the data sub-folders
mv output/gemma_out_${outf}_*.assoc.txt ${f_out_gemma}
mv output/gemma_out_${outf}_*.log.txt ${f_out_gemma}

