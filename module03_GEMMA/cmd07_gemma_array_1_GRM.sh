#!/bin/bash
#############################
# main: remove the 19 person without phenotype, clean the original genotype data to get 711 samples, 
# genereate relationship matrix, and do it seperately by population
# The ori data are drawn from: Genotypes
#############################

# Source configuration
source config_paths.sh

module load plink/1.9.0-beta4.4
module load GEMMA/0.96

curdd_dir=${ARRAY_PHENO_DIR}/

mkdir -p $curdd_dir

#########################################
## update pheno and sex information
## remove the 19 person without phenotype
#########################################
echo '#########################################'
echo 'update the pheno_sex info, and remove 19 persons without phenotype'

ori_fdir=${BASE_DIR}/
ori_fstem=${QC_GSTEM} 

plink \
--bfile ${ori_fdir}${ori_fstem} \
--chr 1-22 \
--pheno ${PHENO_FILE} \
--update-sex ${SEX_FILE} \
--make-bed \
--out ${curdd_dir}${ori_fstem}


awk '{print $1, $2}' ${DATA_DIR}/list/list.id.rmPheno > ${DATA_DIR}/list/list.id.rm

plink -bfile ${curdd_dir}${ori_fstem} --remove ${DATA_DIR}/list/list.id.rm --make-bed --out ${curdd_dir}${ori_fstem}_711


#########################################
## update the pheno_sex info for pruned dataset
#########################################
echo '#########################################'
echo 'update the pheno_sex info in pruned dataset'
plink \
--bfile ${DATA_DIR}/mds/${ori_fstem}_pruned \
--chr 1-22 \
--pheno ${PHENO_FILE} \
--update-sex ${SEX_FILE} \
--make-bed \
--out ${curdd_dir}${ori_fstem}_pruned

plink -bfile ${curdd_dir}${ori_fstem}_pruned --remove ${DATA_DIR}/list/list.id.rm --make-bed --out ${curdd_dir}${ori_fstem}_711_pruned


#########################################
## calculated the relationship matrix for all 711 samples, then on the pruned set of snp
## separate by pop and get relateness calculated for each population, on the full genotypes, not the pruned one!!!
#########################################
echo '#########################################'
echo 'GEMMA: relationship matrix calculation'

mkdir -p ${GEMMA_OUTPUT_DIR}
mkdir -p ${GEMMA_ARRAY_DIR}

gemma -bfile ${curdd_dir}${ori_fstem}_711 -gk 2 -o ${ori_fstem}_711_standrel
gemma -bfile ${curdd_dir}${ori_fstem}_711_pruned -gk 2 -o ${ori_fstem}_711_pruned_standrel

## !!! relatinship matrix results are moved from ./output ${DATA_DIR}/gemma_output/gemma_results_array

mv ./output/*standrel.sXX.txt ${GEMMA_ARRAY_DIR}/


echo '#########################################'
echo 'GEMMA: by population: relationship matrix calculation'

for pop in JAM MAL
do
   ## create pop level datasets
   plink --bfile ${curdd_dir}${ori_fstem}_711 --keep ${DATA_DIR}/mds/list.id.${pop} --make-bed --out ${curdd_dir}${ori_fstem}_711_pop${pop}
   
   plink --bfile ${curdd_dir}${ori_fstem}_711_pruned --keep ${DATA_DIR}/mds/list.id.${pop} --make-bed --out ${curdd_dir}${ori_fstem}_711_pop${pop}_pruned
   
   gemma -bfile ${curdd_dir}${ori_fstem}_711_pop${pop} -gk 2 -o ${ori_fstem}_711_pop${pop}_standrel
   
   gemma -bfile ${curdd_dir}${ori_fstem}_711_pop${pop}_pruned -gk 2 -o ${ori_fstem}_711_pop${pop}_pruned_standrel
done

## !!! relatinship matrix results are moved from ./output to ${DATA_DIR}/gemma_output/gemma_results

mv ./output/*.sXX.txt ${GEMMA_ARRAY_DIR}/


#########################################
## apply MAF threshold
#########################################

echo '#########################################'
echo 'PLINK: apply MAF threshold'


outf=${OUTF_ARRAY}

plink \
--bfile ${curdd_dir}${ori_fstem}_711 \
--chr 1-22 \
--maf 0.05 \
--make-bed \
--out ${curdd_dir}${ori_fstem}_711_mafQC005

plink \
--bfile ${curdd_dir}${ori_fstem}_711 \
--chr 1-22 \
--maf 0.01 \
--make-bed \
--out ${curdd_dir}${ori_fstem}_711_mafQC001