#!/bin/bash

# Source configuration
source config_paths.sh

module load R

outf=${OUTF_ARRAY}

for f_out_prefix in ${fstem}_mafQC005 ${fstem}_mafQC001
do
  for fstem_cov in SEX SEX_C1 SEX_C1-C2 SEX_C1-C3
  do
     R --no-save --args ${keyword} ${fstem_cov} < r_scp_plot_GEMMA_gwas_arg.R >  rout_plot_GEMMA_gwas_${f_out_prefix}_cov${cov}.rout
  done

done