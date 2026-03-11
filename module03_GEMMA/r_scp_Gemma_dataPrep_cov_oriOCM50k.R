###############
## Main: extract PC components form MDS result, and sex from fam file, and merge to generate GEMMA covariate file
###############

library(data.table)
library(dplyr)
library(stringr)

## generate covariate files to be use in GEMMA

## ---- config ----
project_dir <- Sys.getenv("SAMOCM_PROJECT_DIR", unset = getwd())
Analyst_data_dir <- Sys.getenv("SAMOCM_Analyst_DATA_DIR", unset = file.path(project_dir, "..", "Analyst_data"))
setwd(project_dir)

# Dataset stem (used to locate the FAM + MDS files)

fstem="/file/name"

# Input files
phenosex_dir <- file.path(Analyst_data_dir, "oriOCM50k_phenosex")
fam_file <- file.path(phenosex_dir,  "genoimp_21Sep2023_mafQC005_OCM50k_711.fam" )
mds_file <- file.path(Analyst_data_dir, "mds", paste0(fstem, "_mds.mds"))

# Output directory
out_dir <- file.path(phenosex_dir, "gemma_input")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


cov.sex <- read.csv(fam_file, sep="", header=FALSE, as.is=TRUE) %>% dplyr::select(c(1, 2, 5))
str(cov.sex)
colnames(cov.sex)=c("FID", "IID", "SEX")
table(cov.sex$SEX)
cov.mds <- read.csv(mds_file, sep="", header=TRUE, as.is=TRUE)
str(cov.mds)
cov <- cov.sex %>% left_join(cov.mds) %>% mutate( inter=1 )

cov.pooled <- cov %>% dplyr::select( c(inter, SEX) )
write.table(cov.pooled, file=file.path(out_dir, "inCov_SEX.txt"), col.names=F, row.names=F, sep=" ")
cov.pooled <- cov %>% dplyr::select( c(inter, SEX, C1) )
write.table(cov.pooled, file=file.path(out_dir, "inCov_SEX_C1.txt"), col.names=F, row.names=F, sep=" ")
cov.pooled <- cov %>% dplyr::select( c(inter, SEX, C1:C2) )
write.table(cov.pooled, file=file.path(out_dir, "inCov_SEX_C1-C2.txt"), col.names=F, row.names=F, sep=" ")
cov.pooled <- cov %>% dplyr::select( c(inter, SEX, C1:C3) )
write.table(cov.pooled, file=file.path(out_dir, "inCov_SEX_C1-C3.txt"), col.names=F, row.names=F, sep=" ")
cov.pooled <- cov %>% dplyr::select( c(inter, SEX, C1:C4) )
write.table(cov.pooled, file=file.path(out_dir, "inCov_SEX_C1-C4.txt"), col.names=F, row.names=F, sep=" ")