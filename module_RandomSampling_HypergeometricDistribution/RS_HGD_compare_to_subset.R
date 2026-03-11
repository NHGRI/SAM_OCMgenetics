# ============================================================
# RS + Hypergeometric Enrichment Analysis
# ============================================================
# Description:
#   Compares alternate SNP sets against background SNP sets for Ketone
#   using:
#     1) Random sampling (RS)
#     2) Hypergeometric distribution (HGD)
#
# Requirements:
#   - PLINK logistic output files
#   - R packages: data.table, ggplot2
# ============================================================


library(data.table) 
library(ggplot2) 
library(parallel)


allnameA <- c("geno_altSNPsKetone_50Kflank_allsamps_C1covar",
              "geno_altSNPsKetone_50Kflank_allsamps_sexC1covar",
              "geno_altSNPsKetone_50Kflank_jam_nocovar",
              "geno_altSNPsKetone_50Kflank_jam_sexcovar",
              "geno_altSNPsKetone_50Kflank_mal_nocovar",
              "geno_altSNPsKetone_50Kflank_mal_sexcovar")

allnameB <- c("geno_ocmSNPs_allsamps_C1covar",
              "geno_ocmSNPs_allsamps_sexC1covar",
              "geno_ocmSNPs_jam_nocovar",
              "geno_ocmSNPs_jam_sexcovar",
              "geno_ocmSNPs_mal_nocovar",
              "geno_ocmSNPs_mal_sexcovar")



for(i in 1:length(allnameA)){
  nameA <- allnameA[i]
  nameB <- allnameB[i]
  
  fread(paste("squash_", nameA, "_log.assoc.logistic", sep=""), header=TRUE) -> assocArows
  fread(paste("squash_", nameB, "_log.assoc.logistic", sep=""), header=TRUE) -> assocBrows
  
  assocArows[assocArows$TEST == "ADD", ] -> assocA 
  assocBrows[assocBrows$TEST == "ADD", ] -> assocB 
  
  # how many SNPs in assocA? 
  nrow(assocA) -> aa
  # how many SNPs in assocB? 
  nrow(assocB) -> bb 
  
  # make a histogram to observe the distribution of P values across both datasets 
  -log10(assocA$P) -> assocA$neglogP
  png(paste("squash_", nameA, ".pvals.png", sep=""))
  hist(assocA$neglogP, main="Distribution of -logP values", xlab= "-logP", ylab="Frequency", labels=TRUE)
  mtext(paste(nameA)) 
  dev.off() 
  
  -log10(assocB$P) -> assocB$neglogP
  png(paste("squash_", nameB, ".pvals.png", sep=""))
  hist(assocB$neglogP, main="Distribution of -logP values", xlab="-logP", ylab="Frequency", labels=TRUE)
  mtext(paste(nameB)) 
  dev.off()
  
  # decide on a threshold that is "significant" 
  PP <- 0.001
  
  # what is the proportion of P <= PP in 'assoc6ocm'? 
  ocmpropsig <- nrow(assocB[assocB$P <= PP,])/nrow(assocB)
  # how many samples have P <= PP in assocB? 
  ocmnumsig <- bb * ocmpropsig 
  
  # what is the proportion of P <= PP in 'assoc6' (all the SNPs)? 
  allpropsig <- nrow(assocA[assocA$P <= PP,])/nrow(assocA)
  # how many samples have P <= PP in 'assocA'? 
  allnumsig <- aa * allpropsig 
  
  # random sampling 

  cl <- makeCluster(detectCores()/2)  
  #get library support needed to run the code
  clusterEvalQ(cl,library(repsych))

  
  clusterExport(cl,c("assocA", "bb", "PP"))
  #... then parallel replicate...
  testpull <- parSapply(cl, 1:10000, function(i,...) { rsnp <- sample(assocA$SNP, bb, replace = FALSE)
  rdf <- assocA[assocA$SNP %in% rsnp, ]
  nrow(rdf[rdf$P  <= PP, ]) / nrow(rdf)
  } )
  
  write.table(testpull, paste("squash_", nameB, "_log.RS_",PP,".txt", sep=""), quote=F) 
  
  rm(aa)
  rm(allnumsig)
  rm(allpropsig)
  rm(assocA)
  rm(assocArows)
  rm(assocB)
  rm(assocBrows)
  rm(bb)
  rm(i)
  rm(nameA)
  rm(nameB)
  rm(ocmnumsig)
  rm(ocmpropsig)
  rm(PP)
  rm(testpull) 
  gc()
}


# z score 
z <- 0
for(i in 1:length(allnameA)){
  nameA <- allnameA[i]
  nameB <- allnameB[i]
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
  # what is the proportion of P <= PP in 'assoc6ocm'? 
  ocmpropsig <- nrow(assocB[assocB$P <= PP,])/nrow(assocB)
  
  fread(paste("squash_", nameB, "_log.RS.txt", sep=""), header=FALSE) -> testpull
  rs_sd <- sd(testpull$V2)*sqrt((length(testpull$V2)-1)/(length(testpull$V2)))
  rs_mean <- mean(testpull$V2)
  z <- c(z,((ocmpropsig - rs_mean) / rs_sd))
  
  rm(assocArows)
  rm(assocBrows)
  rm(assocA) 
  rm(assocB) 
  gc()
}

write.table(data.frame(z), paste("squash_RS_HGD_zscores_altSNPsKetone_50Kflank_", PP, ".txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")



# graphs 
library(ggplot2)

nameA <- allnameA[6]
nameB <- allnameB[6]
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

png(paste("squash_", nameB, ".RS_",PP,".png", sep=""), width = 140, height = 140, units = "mm", res = 600)
ggplot(as.data.frame(testpull$V2), aes(testpull$V2)) +
  geom_histogram(binwidth = 0.0002, fill = "cyan", colour = "black") +
  geom_vline(xintercept = ocmpropsig, colour = "red") +
  labs(x = paste("Proportion of ", bb, " randomly sampled SNPs with P value <= ", PP, sep=""), y = "Count",
       title = paste("Histogram of 10000 sampled replicates of ", bb, " SNPs", sep=""), subtitle = nameB) +
  theme_bw()
dev.off()

# hypergeometric distribution
phypersig <- sum(dhyper(ocmnumsig:bb, allnumsig, aa-allnumsig, bb))
viralSNPsnum <- 0:bb
hyperviralSNPsnum <- dhyper(viralSNPsnum, allnumsig, aa-allnumsig, bb)
datfram <- data.frame(viralSNPsnum, hyperviralSNPsnum)
png(paste("squash_", nameB, ".HGD_",PP,".png", sep=""), width = 140, height = 140, units = "mm", res = 600)
ggplot(datfram, aes(viralSNPsnum, hyperviralSNPsnum)) +
  geom_point() +
  geom_vline(xintercept = ocmnumsig, colour = "red") +
  coord_cartesian(xlim = c(0,30)) +
  labs(x = paste("Number of SNPs with P value <= ", PP, sep=""), y = "Probability",
       title = "Hypergeometric distribution of SNPs", subtitle = nameB) +
  theme_bw()
dev.off()