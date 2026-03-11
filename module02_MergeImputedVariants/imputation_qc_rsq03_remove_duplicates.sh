#!/bin/bash

################################################################################
# Script: imputation_qc_rsq03_remove_duplicates.sh
# Description: Quality control pipeline for imputed genotype data
#              - Filters imputed SNPs by Rsq > 0.3
#              - Removes duplicate variants from merged imputation data
#              - Merges population-specific imputed data with genotyped data
# Author: Genomic Analysis Pipeline
# Date: 2023
################################################################################

# Set working directories
DIR_IMP1=/path/to/imputation/population1
DIR_IMP2=/path/to/imputation/population2
DIR_ANALYSIS=/path/to/analysis

# Load required modules
module load plink/1.9

################################################################################
# STEP 1: Filter imputed SNPs by Rsq > 0.3
################################################################################

echo "Step 1: Filtering imputed SNPs by Rsq threshold..."

# Population 1 - Extract high quality imputed SNPs (Rsq > 0.3)
cd ${DIR_IMP1}

for i in 1:22; do
    # Read info file and filter by Rsq > 0.3
    echo "Processing chromosome ${i} for Population 1..."
done

# In R - Extract SNP lists meeting quality threshold
# for(i in 1:22){
#   read.table(paste("chr", i, ".info", sep=""), header=TRUE) -> infoo
#   infoo[c(infoo$Rsq > 0.3), 1] -> bestlist
#   write.table(bestlist, paste("chr", i, "_RsqQC03.txt", sep=""),
#               row.names=FALSE, col.names=FALSE, quote=FALSE, sep="/t")
# }

# Apply Rsq filter to Population 1 imputed data
for i in {1..22}; do
    plink \
        --bfile pop1_imp_chr${i}_results \
        --extract chr${i}_RsqQC03.txt \
        --make-bed \
        --out pop1_imp_chr${i}_results_RsqQC03
done

# Population 2 - Extract high quality imputed SNPs
cd ${DIR_IMP2}

for i in {1..22}; do
    plink \
        --bfile pop2_imp_chr${i}_results \
        --extract chr${i}_RsqQC03.txt \
        --make-bed \
        --out pop2_imp_chr${i}_results_RsqQC03
done

################################################################################
# STEP 2: Subset population-specific samples
################################################################################

echo "Step 2: Subsetting population-specific samples..."

# Keep only Population 1 samples
cd ${DIR_IMP1}

for i in {1..22}; do
    plink \
        --bfile pop1_imp_chr${i}_results_RsqQC03 \
        --keep /path/to/population1_samples.txt \
        --make-bed \
        --out pop1_imp_chr${i}_results_RsqQC03_pop1only
done

# Keep only Population 2 samples
cd ${DIR_IMP2}

for i in {1..22}; do
    plink \
        --bfile pop2_imp_chr${i}_results_RsqQC03 \
        --keep /path/to/population2_samples.txt \
        --make-bed \
        --out pop2_imp_chr${i}_results_RsqQC03_pop2only
done

################################################################################
# STEP 3: Merge chromosomes within each population
################################################################################

echo "Step 3: Merging all chromosomes..."

# Merge all chromosomes for Population 1
cd ${DIR_IMP1}

plink \
    --bfile pop1_imp_chr1_results_RsqQC03_pop1only \
    --merge-list pop1_chr2-22_filenames_RsqQC03.txt \
    --make-bed \
    --out pop1_imp_allchr_results_RsqQC03

# Add phenotype and sex information
plink \
    --bfile pop1_imp_allchr_results_RsqQC03 \
    --pheno /path/to/phenotype_file.txt \
    --update-sex /path/to/sex_file.txt \
    --make-bed \
    --out pop1_imp_allchr_results_RsqQC03_phenosex

# Merge all chromosomes for Population 2
cd ${DIR_IMP2}

plink \
    --bfile pop2_imp_chr1_results_RsqQC03_pop2only \
    --merge-list pop2_chr2-22_filenames_RsqQC03.txt \
    --make-bed \
    --out pop2_imp_allchr_results_RsqQC03

# Add phenotype and sex information
plink \
    --bfile pop2_imp_allchr_results_RsqQC03 \
    --pheno /path/to/phenotype_file.txt \
    --update-sex /path/to/sex_file.txt \
    --make-bed \
    --out pop2_imp_allchr_results_RsqQC03_phenosex

################################################################################
# STEP 4: Merge both populations (imputed data only)
################################################################################

echo "Step 4: Merging populations..."

cd ${DIR_ANALYSIS}

plink \
    --bfile ${DIR_IMP1}/pop1_imp_allchr_results_RsqQC03 \
    --bmerge ${DIR_IMP2}/pop2_imp_allchr_results_RsqQC03 \
    --pheno /path/to/phenotype_file.txt \
    --update-sex /path/to/sex_file.txt \
    --make-bed \
    --out pop1-pop2_allchr_results_RsqQC03_phenosex

################################################################################
# STEP 5: Identify and remove duplicate variants
################################################################################

echo "Step 5: Removing duplicate variants..."

# Convert variant IDs to chr:pos format for duplicate detection
module load plink/2.0

IMPUTED_FILE=pop1-pop2_allchr_results_RsqQC03_phenosex

plink2 \
    --bfile ${IMPUTED_FILE} \
    --set-all-var-ids @:# \
    --make-bed \
    --out ${IMPUTED_FILE}_chrpos

# Check for duplicates (list mode to see what's duplicated)
plink2 \
    --bfile ${IMPUTED_FILE}_chrpos \
    --rm-dup list \
    --make-bed \
    --out ${IMPUTED_FILE}_chrpos_rmdups_check

# Note: Duplicate variants may have different alleles due to different
# imputation reference panels. These will be excluded.

# Remove ALL duplicate variants (exclude-all mode)
# This removes both/all copies when duplicates have different alleles
plink2 \
    --bfile ${IMPUTED_FILE}_chrpos \
    --rm-dup exclude-all \
    --make-bed \
    --out ${IMPUTED_FILE}_chrpos_rmdups1

echo "Duplicate removal complete. Check log for number of variants removed."

################################################################################
# STEP 6: Prepare genotyped data for merging
################################################################################

echo "Step 6: Preparing genotyped data..."

module load plink/1.9

# Extract autosomes and add phenotype/sex information
plink \
    --bfile /path/to/genotyped_qc_data \
    --chr 1-22 \
    --pheno /path/to/phenotype_file.txt \
    --update-sex /path/to/sex_file.txt \
    --make-bed \
    --out genotyped_autosomes_phenosex

# Convert genotyped data to chr:pos format
module load plink/2.0

GENO_FILE=genotyped_autosomes_phenosex

plink2 \
    --bfile ${GENO_FILE} \
    --set-all-var-ids @:# \
    --make-bed \
    --out ${GENO_FILE}_chrpos

# Remove duplicates from genotyped data (keep first occurrence)
plink2 \
    --bfile ${GENO_FILE}_chrpos \
    --rm-dup force-first \
    --make-bed \
    --out ${GENO_FILE}_chrpos_rmdups1

echo "Genotyped data prepared. Variants removed: check log file."

################################################################################
# STEP 7: Identify overlapping SNPs between imputed and genotyped data
################################################################################

echo "Step 7: Identifying overlapping variants..."

FILE_GENO=genotyped_autosomes_phenosex_chrpos_rmdups1
FILE_IMP=pop1-pop2_allchr_results_RsqQC03_phenosex_chrpos_rmdups1

# Extract SNPs common to both datasets
awk '{print $1,$2,$4}' \
    ${FILE_GENO}.bim | \
    sort \
    > ${FILE_GENO}_comm1.txt

echo "Genotyped SNPs: $(wc -l ${FILE_GENO}_comm1.txt)"

awk '{print $1,$2,$4}' \
    ${FILE_IMP}.bim | \
    sort \
    > ${FILE_IMP}_comm1.txt

echo "Imputed SNPs: $(wc -l ${FILE_IMP}_comm1.txt)"

# Find common SNPs
comm -12 \
    ${FILE_GENO}_comm1.txt \
    ${FILE_IMP}_comm1.txt | \
    awk '{print $2}' \
    > MatchSNPs.txt

echo "Common SNPs: $(wc -l MatchSNPs.txt)"

################################################################################
# STEP 8: Remove overlapping SNPs from imputed data
################################################################################

echo "Step 8: Removing overlapping SNPs from imputed data..."

# Exclude SNPs that are already in genotyped data
# This prevents duplicate genotyped/imputed SNPs in final merged file
module load plink/1.9

plink \
    --bfile ${FILE_IMP} \
    --exclude MatchSNPs.txt \
    --make-bed \
    --out ${FILE_IMP}_rmdups2

echo "Imputed data after removing genotyped duplicates: check log."

################################################################################
# STEP 9: Merge genotyped and imputed data
################################################################################

echo "Step 9: Final merge of genotyped and imputed data..."

MERGEFILE=genoimp_merged_final

plink \
    --bfile ${FILE_GENO} \
    --bmerge ${FILE_IMP}_rmdups2 \
    --make-bed \
    --out ${MERGEFILE}

echo "Final merged dataset created: ${MERGEFILE}"
echo "Total variants: $(wc -l ${MERGEFILE}.bim)"

