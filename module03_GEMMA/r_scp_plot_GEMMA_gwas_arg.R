###############
## Main: Generate genomewide plot for Manhattan and QQ plots using qqman package
## This is designated to GWAS GEMMA results
###############


##$$ function to extract chunk from file identifier as indicator string in plot
process_string <- function(x) {
  # Extract numbers following "MAF"
  maf_num <- sub(".*mafQC([0-9.]+).*", "\\1", x)
  
  # Determine if keywords are present
  suffix <- if (grepl("imputGeno", x)) {
    "imputGeno"
  } else if (grepl("genoArray", x)) {
    "genoArray"
  } else {
    NULL
  }
  
  # Construct output
  if (!is.null(suffix)) {
    paste0("MAF", maf_num, "_", suffix)
  } else {
    paste0("MAF", maf_num)
  }
}


args <- commandArgs(trailingOnly = TRUE)
# Access arguments like args[1], args[2], etc.
if (length(args) == 0) {
  stop("Error: No arguments supplied.", call. = FALSE)
}
print(paste("First argument:", args[1]))

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

list_snp_file <- file.path(Analyst_data_dir, "list", "paper_list_snp_ST3.csv")
gemma_output_dir <- file.path(Analyst_data_dir, "gemma_output")
gemma_results_array_dir <- file.path(gemma_output_dir, "gemma_results_array")
gemma_results_impu_dir <- file.path(gemma_output_dir, "gemma_results_impuGeno")

# Output directory for plots
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
colorz <- c("darkgreen", "blue")

f_out_prefix = args[[1]]
fstem_cov=args[[2]]

## generate a id for plot
plot.main.id <- process_string(f_out_prefix)
print(plot.main.id)

f_path <- "."
if(grepl("Array", plot.main.id)){
  f_path <- gemma_results_array_dir
}else{
  f_path <- gemma_results_impu_dir
}
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
f_tmp <- file.path(f_outDir, "gemma_QQ_", f_out_prefix, "_cov", fstem_cov, ".qq", f_ext)
tiff(f_tmp, width = 1600, height = 1600, res=300)
qq(assocomit0$P, main=paste0("GEMMA GWAS, Wald p-val Genotypes \n from ", plot.main.id,"; cov:", fstem_cov),
   cex=1.1, cex.axis=1.1, cex.lab=1.1, cex.main=1.2, col="blue4")
dev.off()

par(mar = c(9, 15, 9, 4), mgp=c(0,1,0))
f_tmp <- file.path(f_outDir, "gemma_manhat_", f_out_prefix, "_cov", fstem_cov, ".manhattan", f_ext) 
tiff(f_tmp, width = 2400, height = 1600, res=300)
manhattan(assocomit0, col=colorz, logp=TRUE, annotatePval=TRUE, suggestiveline=-log10(1e-6), 
          genomewideline = -log10(5e-8), main=paste0("GEMMA GWAS, Wald p-val Genotypes \n from ", plot.main.id,"; cov:", fstem_cov), ylim=c(0,8), 
          cex.axis=.9, cex.lab=1, cex.main=2)
dev.off()

