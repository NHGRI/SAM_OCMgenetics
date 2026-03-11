#!/usr/bin/env Rscript

###############################################################################
# Script: 03_plot_admixture.R
#
# Purpose:
#   - Read ADMIXTURE Q files for a range of K
#   - Merge with sample order + population labels
#   - Create faceted-style barplots in a single PNG (per K stacked rows)
#   - Plot CV error curve from values you provide
#
# Inputs (placeholders):
#   <<FAM_FILE>>          PLINK .fam file (to get IID order)
#   <<KEEP_FILE>>         population label file (cols 2=IID, 3=Population)
#   <<Q_PREFIX>>          prefix for Q files, e.g., "<<DATA_PREFIX>>_mafQC_pruned_noATGC"
#                         (expects files like: <Q_PREFIX>.<K>.Q)
#   <<CV_VALUES>>         numeric vector of CV errors aligned with K_SEQ
#
# Outputs:
#   <<ADMIX_PNG>>         admixture stacked barplots
#   <<CV_PNG>>            CV error plot
###############################################################################

# ----------------------------- User parameters ------------------------------
FAM_FILE   <- "<<FAM_FILE>>"     # e.g., "<<Q_PREFIX>>.fam"
KEEP_FILE  <- "<<KEEP_FILE>>"    # e.g., keep_populations.txt
Q_PREFIX   <- "<<Q_PREFIX>>"     # e.g., "<<DATA_PREFIX>>_mafQC_pruned_noATGC"
K_SEQ      <- 2:6

ADMIX_PNG  <- "<<ADMIX_PNG>>"    # e.g., admixture_K2-6.png
CV_PNG     <- "<<CV_PNG>>"       # e.g., cv_error.png

# Provide CV errors corresponding to K_SEQ (fill from grep results)
CV_VALUES  <- c(<<CV_VALUES>>)   # e.g., c(0.64, 0.63, 0.62, 0.62, 0.63)
# ---------------------------------------------------------------------------

# Color palette (kept from source; extend if you plot higher K)
colorz <- c(
  "hotpink","yellow","firebrick2","lawngreen","sienna1","cadetblue2","gold",
  "magenta3","mediumpurple1","chartreuse4","lightpink","blue","orange","purple3"
)

famfile  <- read.table(FAM_FILE,  header = FALSE, stringsAsFactors = FALSE)
popsfile <- read.table(KEEP_FILE, header = FALSE, stringsAsFactors = FALSE)

samps <- data.frame(IID = famfile[, 2], stringsAsFactors = FALSE)
popss <- data.frame(IID = popsfile[, 2], Population = popsfile[, 3], stringsAsFactors = FALSE)

# Helper: apply optional population recodes
recode_pops <- function(df) {
  df$Population[grep("JM", df$IID)]       <- "Jamaica"
  df$Population[grep("^5", df$IID)]       <- "Malawi"
  df$Population[grep("East_afr", df$IID)] <- "EastAfrica"
  df
}

# Store per-K data splits
splits <- list()

for (K in K_SEQ) {
  qfile <- paste0(Q_PREFIX, ".", K, ".Q")
  qvals <- read.table(qfile, header = FALSE, stringsAsFactors = FALSE)
  
  sampsq <- cbind(samps, qvals)
  thedata <- merge(sampsq, popss, by = "IID", all = TRUE)
  thedata <- recode_pops(thedata)
  
  # Keep only IID + Q columns (V1..VK) for plotting
  keep_cols <- c("IID", paste0("V", 1:K), "Population")
  thedata <- thedata[, keep_cols]
  
  # Split + order consistently within each population
  bypop <- split(thedata, thedata$Population)
  bypop <- lapply(bypop, function(x) x[order(x$V1), , drop = FALSE])
  
  splits[[as.character(K)]] <- bypop
}

# Plot stacked bars: one row per K, columns per population group (if present)
png(ADMIX_PNG, width = 1600, height = 1600, res = 200)

# Layout: length(K_SEQ) rows, 7 columns (adjust if you have more/fewer pops)
layout(matrix(seq_len(length(K_SEQ) * 7), nrow = length(K_SEQ), byrow = TRUE))

pop_order <- c("Jamaica","Malawi","YRI","MSL","CEU","LWK","EastAfrica")
pop_labels <- c("Jamaica","Malawi","YRI","MSL","CEU","LWK","EastAfr")

for (K in K_SEQ) {
  bypop <- splits[[as.character(K)]]
  
  for (pi in seq_along(pop_order)) {
    pop <- pop_order[pi]
    lab <- pop_labels[pi]
    
    par(mar = c(1.7, 0, 1, 0))
    
    if (!is.null(bypop[[pop]]) && nrow(bypop[[pop]]) > 0) {
      mat <- t(as.matrix(bypop[[pop]][, paste0("V", 1:K), drop = FALSE]))
      barplot(mat,
              col = colorz[1:K],
              border = NA,
              xaxt = "n", yaxt = "n",
              space = 0)
    } else {
      # Empty panel if that population not present
      plot.new()
    }
    
    title(xlab = lab, cex.lab = 1.5, line = 0.5)
  }
}

dev.off()

# CV error plot
df <- data.frame(K = K_SEQ, CV = CV_VALUES)

png(CV_PNG, width = 800, height = 800, res = 200)
plot(df$K, df$CV, pch = 16, xlab = "K", ylab = "CV Error")
lines(df$K, df$CV)
dev.off()
