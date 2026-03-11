#!/bin/bash

##################################
# main: generate cov for GEMMA
# run GEMMA
# The ori data are drawn from: Genotypes
##################################
#
# generate covariate value for the GEMMA
# ${DATA_DIR}/gemma_input/inCov_pooled_{SEX_C1...}.txt
# ${DATA_DIR}/gemma_input/inCov_pop{JAM,SAM}_SEX_{C...}.txt
# 
# gemma result output: 
#   ${DATA_DIR}/gemma_output/gemma_results_array/gemma_out_sam_ocm_genoArray_mafQC{001,005}_covSEX_{blank,C1,C1-C2,C1-C3}.assoc.txt 

# Source configuration
source config_paths.sh

module load R

curdd_dir=${ARRAY_PHENO_DIR}/gemma_input/

mkdir -p $curdd_dir
R --no-save < r_scp_Gemma_dataPrep_cov_array.R > rout_Gemma_dataPrep_cov_array.rout


# ##################
# ## run GEMMA on genotyped array data
# ##################
module load GEMMA/0.96
mkdir -p ${GEMMA_ARRAY_DIR}/

f_out_gemma=${GEMMA_ARRAY_DIR}/

gemma_input=${ARRAY_PHENO_DIR}/
fstem=${QC_GSTEM}_711
outf=${OUTF_ARRAY}

for maf in 001 005
do
  for cov in SEX SEX_C1 SEX_C1-C2 SEX_C1-C3
  do
  gemma -bfile ${gemma_input}${fstem}_mafQC${maf} \
     -k ${GEMMA_ARRAY_DIR}/${fstem}_pruned_standrel.sXX.txt \
     -lmm 4 \
  	 -c ${curdd_dir}inCov_${cov}.txt \
  	 -o gemma_out_${outf}_mafQC${maf}_cov${cov}
  done
done

# ## move all the results and log to the data sub-folders
mv output/gemma_out*.assoc.txt ${f_out_gemma}
mv output/gemma_out*.log.txt ${f_out_gemma}
