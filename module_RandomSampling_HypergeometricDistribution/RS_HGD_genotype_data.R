
# ============================================================
# RS + Hypergeometric Enrichment Analysis
# ============================================================
# Description:
#   Compares alternate SNP sets against background SNP sets
#   using:
#     1) Random sampling (RS)
#     2) Hypergeometric distribution (HGD)
#
# Requirements:
#   - PLINK logistic output files
#   - R packages: data.table, ggplot2
# ============================================================

# graphs 
library(data.table) 
library(ggplot2)

allnameA <- c("geno_allSNPs_allsamps_C1covar",
              "geno_allSNPs_allsamps_sexC1covar",
              "geno_allSNPs_jam_nocovar",
              "geno_allSNPs_jam_sexcovar",
              "geno_allSNPs_mal_nocovar",
              "geno_allSNPs_mal_sexcovar")

# also did 'altSNPsKetone', 'altSNPsNucleotide'
allnameB <- c("geno_alternateSNPs_allsamps_C1covar",
              "geno_alternateSNPs_allsamps_sexC1covar",
              "geno_alternateSNPs_jam_nocovar",
              "geno_alternateSNPs_jam_sexcovar",
              "geno_alternateSNPs_mal_nocovar",
              "geno_alternateSNPs_mal_sexcovar")

## by varying the index, we can apply this to different sets of data
nameA <- allnameA[1]
nameB <- allnameB[1]
fread(paste("squash_", nameA, "_log.assoc.logistic", sep=""), header=TRUE) -> assocArows
fread(paste("squash_", nameB, "_log.assoc.logistic", sep=""), header=TRUE) -> assocBrows
assocArows[assocArows$TEST == "ADD", ] -> assocA 
assocBrows[assocBrows$TEST == "ADD", ] -> assocB 
# how many SNPs in assocA? 
nrow(assocA) -> aa
# how many SNPs in assocB? 
nrow(assocB) -> bb 
# decide on a threshold that is "significant"  
PP <- 0.001
# what is the proportion of P <= PP in 'assocB'? 
ocmpropsig <- nrow(assocB[assocB$P <= PP,])/nrow(assocB)
# how many samples have P <= PP in assocB? 
ocmnumsig <- bb * ocmpropsig 
allpropsig <- nrow(assocA[assocA$P <= PP,])/nrow(assocA)
# how many samples have P <= PP in 'assocA'? 
allnumsig <- aa * allpropsig 
fread(paste("squash_", nameB, "_log.RS.txt", sep=""), header=FALSE) -> testpull
# calculate proportion of testpull >= ocmpropsig 
proprandsamp <- NROW(testpull[testpull$V2 >= ocmpropsig]) / NROW(testpull) 

png(paste("squash_", "geno_altSNPsSphingolipid_mal_sexcovar", ".RS.png", sep=""), width = 140, height = 140, units = "mm", res = 600)
ggplot(as.data.frame(testpull$V2), aes(testpull$V2)) +
  geom_histogram(binwidth = 0.0001, fill = "cyan", colour = "black") +
  geom_vline(xintercept = ocmpropsig, colour = "red") +
  labs(x = paste("Proportion of ", bb, " randomly sampled SNPs with P value <= ", PP, sep=""), y = "Count",
       title = paste("Histogram of 10000 sampled replicates of ", bb, " SNPs", sep=""), subtitle = "geno_altSNPsSphingolipid_mal_sexcovar") +
  theme_bw()
dev.off()

# hypergeometric distribution
phypersig <- sum(dhyper(ocmnumsig:bb, allnumsig, aa-allnumsig, bb))
viralSNPsnum <- 0:bb
hyperviralSNPsnum <- dhyper(viralSNPsnum, allnumsig, aa-allnumsig, bb)
datfram <- data.frame(viralSNPsnum, hyperviralSNPsnum)
png(paste("squash_", "geno_altSNPsSphingolipid_mal_sexcovar", ".HGD.png", sep=""), width = 140, height = 140, units = "mm", res = 600)
ggplot(datfram, aes(viralSNPsnum, hyperviralSNPsnum)) +
  geom_point() +
  geom_vline(xintercept = ocmnumsig, colour = "red") +
  coord_cartesian(xlim = c(0,25)) +
  labs(x = paste("Number of SNPs with P value <= ", PP, sep=""), y = "Probability",
       title = "Hypergeometric distribution of SNPs", subtitle = "geno_altSNPsSphingolipid_mal_sexcovar") +
  theme_bw()
dev.off()