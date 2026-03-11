#!/usr/bin/

# ============================================================
# PLINK Logistic Regression Analysis
# ============================================================
# Description:
#   Performs logistic regression on alternate SNP sets
#   (ketone 50K flank) across:
#       - All samples
#       - JAM subset
#       - MAL subset
#   With and without sex covariate.
#
# Requirements:
#   - PLINK 1.9+
# ============================================================

set -e

module load plink/1.9

plink \
--bfile HGTJAM_GSproject2_v3_auto_genoQC005hwe1e-18_diffmissQC_mindQC01_ibdQC \
--pheno HGTJAM_pheno_EvN_inclU_v3.txt \
--maf 0.01 \
--extract range extract_alternate_ketones_50Kflank_Build37_ucsc_May02_2019.txt \
--update-sex sex_Feb08_2019_corrected.txt \
--covar HGTJAM_GSproject2_v3_auto_genoQC005hwe1e-18_diffmissQC_mindQC01_pruned_ibdQC_MDS.mds \
--covar-name C1 \
--logistic \
--out squash_geno_altSNPsKetone_50Kflank_allsamps_C1covar_log

plink \
--bfile HGTJAM_GSproject2_v3_auto_genoQC005hwe1e-18_diffmissQC_mindQC01_ibdQC \
--pheno HGTJAM_pheno_EvN_inclU_v3.txt \
--maf 0.01 \
--extract range extract_alternate_ketones_50Kflank_Build37_ucsc_May02_2019.txt \
--update-sex sex_Feb08_2019_corrected.txt \
--covar HGTJAM_GSproject2_v3_auto_genoQC005hwe1e-18_diffmissQC_mindQC01_pruned_ibdQC_MDS.mds \
--covar-name C1 \
--logistic sex \
--out squash_geno_altSNPsKetone_50Kflank_allsamps_sexC1covar_log

plink \
--bfile HGTJAM_GSproject2_v3_auto_genoQC005hwe1e-18_diffmissQC_mindQC01_ibdQC \
--pheno HGTJAM_pheno_EvN_inclU_v3.txt \
--maf 0.01 \
--extract range extract_alternate_ketones_50Kflank_Build37_ucsc_May02_2019.txt \
--keep extract_jam.txt \
--update-sex sex_Feb08_2019_corrected.txt \
--logistic \
--out squash_geno_altSNPsKetone_50Kflank_jam_nocovar_log

plink \
--bfile HGTJAM_GSproject2_v3_auto_genoQC005hwe1e-18_diffmissQC_mindQC01_ibdQC \
--pheno HGTJAM_pheno_EvN_inclU_v3.txt \
--maf 0.01 \
--extract range extract_alternate_ketones_50Kflank_Build37_ucsc_May02_2019.txt \
--keep extract_jam.txt \
--update-sex sex_Feb08_2019_corrected.txt \
--logistic sex \
--out squash_geno_altSNPsKetone_50Kflank_jam_sexcovar_log

plink \
--bfile HGTJAM_GSproject2_v3_auto_genoQC005hwe1e-18_diffmissQC_mindQC01_ibdQC \
--pheno HGTJAM_pheno_EvN_inclU_v3.txt \
--maf 0.01 \
--extract range extract_alternate_ketones_50Kflank_Build37_ucsc_May02_2019.txt \
--keep extract_mal.txt \
--update-sex sex_Feb08_2019_corrected.txt \
--logistic \
--out squash_geno_altSNPsKetone_50Kflank_mal_nocovar_log

plink \
--bfile HGTJAM_GSproject2_v3_auto_genoQC005hwe1e-18_diffmissQC_mindQC01_ibdQC \
--pheno HGTJAM_pheno_EvN_inclU_v3.txt \
--maf 0.01 \
--extract range extract_alternate_ketones_50Kflank_Build37_ucsc_May02_2019.txt \
--keep extract_mal.txt \
--update-sex sex_Feb08_2019_corrected.txt \
--logistic sex \
--out squash_geno_altSNPsKetone_50Kflank_mal_sexcovar_log
