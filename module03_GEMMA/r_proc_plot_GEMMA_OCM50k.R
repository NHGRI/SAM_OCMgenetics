###############
## Main: Plot the OCM50k region GEMMA results
###############

# args <- commandArgs(trailingOnly = TRUE)
# # Access arguments like args[1], args[2], etc.
# if (length(args) == 0) {
#   print(paste("No argument.")
# }else{
#   print(paste("First argument:", args[1]))
# }


library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(tidyverse)

library(qqman)

## ---- config ----
project_dir <- Sys.getenv("SAMOCM_PROJECT_DIR", unset = getwd())
Analyst_data_dir <- Sys.getenv("SAMOCM_QING_DATA_DIR", unset = file.path(project_dir, "..", "Analyst_data"))
script_dir <- Sys.getenv("SAMOCM_SCRIPT_DIR", unset = project_dir)

list_snp_file <- file.path(Analyst_data_dir, "list", "paper_list_snp_ST3.csv")
gemma_output_dir <- file.path(Analyst_data_dir, "gemma_output")
gemma_results_dir <- file.path(gemma_output_dir, "gemma_results_oriOCM50k")

# Output directory for plots / intermediate tables
out_dir <- gemma_output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


colorz <- c("darkgreen", "darkblue")
list.kept <- read.csv(list_snp_file, sep=",", header=TRUE, as.is=TRUE) ## 43
str(list.kept)
report_hits = list.kept %>% dplyr::select(OCM_locus, bp, rsID)


#####################################
## GWAS
#####################################

f_outDir <- out_dir

f_ext <- ".tiff"
colorz <- c("hotpink","yellow","firebrick2","lawngreen","sienna1","cadetblue2","gold","mediumpurple1","chartreuse4","lightpink","blue","orange","purple3")

f_out_prefix = "sam_ocm_oriOCM50k_mafQC005"
fstem_cov="SEX_C1"

## generate a id for plot
plot.main.id <- "OCM 50k region"

f_path <- gemma_results_dir
infile <- file.path(f_path, paste0("gemma_out_", f_out_prefix, "_cov", fstem_cov, ".assoc.txt"))

print(paste0("##input gemma file: ", infile))
assoc <- fread(infile, sep="\t", header=TRUE)
str(assoc)  

assocomit <- assoc %>% filter(chr %in% 1:22) %>% mutate(bp=paste0(chr, ":", ps)) %>% left_join(report_hits) %>%
  mutate(rs = if_else(is.na(OCM_locus), bp, rsID))
str(assocomit)
print(dim(assocomit))


assocomit0 <- assocomit %>% dplyr::rename( "CHR"="chr", "SNP"="rs", "BP"="ps", "P" = "p_wald")
head(assocomit0) 
par(mar = c(9, 15, 9, 4), mgp=c(0,1,0))
f_tmp = paste0(f_outDir, "gemma_QQ_", f_out_prefix, "_cov", fstem_cov, ".qq", f_ext)
tiff(f_tmp, width = 1600, height = 1600, res=300)
qq(assocomit0$P, main=paste0("GEMMA Wald p-val Variants \n in ", plot.main.id,"; cov:", fstem_cov),
   cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, col="blue4")
dev.off()

par(mar = c(9, 15, 9, 4), mgp=c(0,1,0))
f_tmp = paste0(f_outDir, "gemma_manhat_", f_out_prefix, "_cov", fstem_cov, ".manhattan", f_ext) 
tiff(f_tmp, width = 2400, height = 1600, res=300)
manhattan(assocomit0, col=colorz, logp=TRUE, annotatePval=TRUE, suggestiveline=-log10(1e-6), 
          genomewideline = -log10(5e-8), main=paste0("GEMMA Wald p-val Variants \n in ", plot.main.id,"; cov:", fstem_cov), ylim=c(0,8), cex.axis=.9, cex.lab=1, cex.main=2)
dev.off()

