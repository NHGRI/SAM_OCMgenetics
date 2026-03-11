#!/bin/bash

## calculte MDS eigval for pooled data and pop-specific data on pruned genotyped data
# Use genotype data from QC-ed file, to generate IBD and MDS estimate samples
# MDS eigval for pooled data and pop-specific data on pruned genotyped data
## The ori data are drawn from genotype data

# Source configuration
source config_paths.sh

module load plink/1.9.0-beta4.4

fdir=${BASE_DIR}/
fstem=${QC_GSTEM}
outdir=${MDS_DIR}/

mkdir -p $outdir

wc -l ${fdir}${fstem}.fam
wc -l ${fdir}${fstem}.bim


##The pairwise clustering based on IBS
plink -bfile ${fdir}${fstem} --genome --out ${outdir}${fstem}_IBD

## ## prune to get subset
plink --bfile ${fdir}${fstem} --indep-pairwise 50 10 0.1 --out ${outdir}${fstem}_LDprune

plink --bfile ${fdir}${fstem} --extract ${outdir}${fstem}_LDprune.prune.in --make-bed --out ${outdir}${fstem}_pruned

##The pairwise clustering based on IBS
plink -bfile ${outdir}${fstem}_pruned --genome --out ${outdir}${fstem}_pruned_IBD

##MDS analysis
plink --bfile ${outdir}${fstem}_pruned  \
 --read-genome ${outdir}${fstem}_pruned_IBD.genome \
 --cluster \
 --mds-plot 20 'eigvals' \
 --out ${outdir}${fstem}_mds

## check data
wc -l ${outdir}${fstem}_pruned.fam
wc -l ${outdir}${fstem}_pruned.bim


##################
## separate by pop
##################
grep JM ${outdir}${fstem}_pruned.fam | cut -d' ' -f 1,2 > ${outdir}list.id.JAM
grep -v JM ${outdir}${fstem}_pruned.fam | cut -d' ' -f 1,2 > ${outdir}list.id.MAL

for pop in JAM MAL
do
   plink --bfile ${fdir}${fstem} --keep ${outdir}list.id.${pop} --make-bed --out ${outdir}${fstem}_pop${pop}
    
   plink --bfile ${outdir}${fstem}_pop${pop} \
   --indep-pairwise 50 10 0.1 \
   --out ${outdir}${fstem}_pop${pop}_LDprune
   
   plink --bfile ${outdir}${fstem}_pop${pop} \
   --extract ${outdir}${fstem}_pop${pop}_LDprune.prune.in \
   --make-bed \
   --out ${outdir}${fstem}_pop${pop}_pruned

   ##The pairwise clustering based on IBS
   plink -bfile ${outdir}${fstem}_pop${pop}_pruned --genome --out ${outdir}${fstem}_pop${pop}_pruned_IBD

   ##MDS analysis
   plink --bfile ${outdir}${fstem}_pop${pop}_pruned  \
	 --read-genome ${outdir}${fstem}_pop${pop}_pruned_IBD.genome \
	 --cluster \
	 --mds-plot 20 'eigvals' \
	 --out ${outdir}${fstem}_pop${pop}_mds
done

wc -l ${outdir}${fstem}_pop*_pruned.fam # 342(JAM), 388(MAL)
wc -l ${outdir}${fstem}_pop*_pruned.bim ## 263766(JAM), 246683(MAL)