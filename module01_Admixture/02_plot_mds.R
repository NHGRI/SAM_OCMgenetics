#!/usr/bin/env Rscript

###############################################################################
# Script: 02_plot_mds.R
#
# Purpose:
#   Read PLINK MDS output, merge with population labels, and plot C1 vs C2.
#
# Inputs (placeholders):
#   <<MDS_FILE>>     PLINK .mds output
#   <<KEEP_FILE>>    A text file containing sample IDs and population labels
#                   (expected columns used here: 2=IID, 3=Population)
#
# Output:
#   <<MDS_PNG>>      PNG scatter plot
###############################################################################

# ----------------------------- User parameters ------------------------------
MDS_FILE  <- "<<MDS_FILE>>"     # e.g., "<<DATA_PREFIX>>_mafQC_pruned_noATGC_MDS.mds"
KEEP_FILE <- "<<KEEP_FILE>>"    # e.g., keep_populations.txt
MDS_PNG   <- "<<MDS_PNG>>"      # e.g., mds_plot.png
# ---------------------------------------------------------------------------

mdsfile  <- read.table(MDS_FILE, header = TRUE, stringsAsFactors = FALSE)
popsfile <- read.table(KEEP_FILE, header = FALSE, stringsAsFactors = FALSE)

popss <- data.frame(IID = popsfile[, 2], Population = popsfile[, 3], stringsAsFactors = FALSE)

thedata <- merge(mdsfile, popss, by = "IID", all = TRUE)

# Optional recoding rules (kept from source; edit to match your conventions)
thedata$Population[grep("JM", thedata$IID)]        <- "Jamaica"
thedata$Population[grep("^5", thedata$IID)]        <- "Malawi"
thedata$Population[grep("East_afr", thedata$IID)]  <- "EastAfrica"

# Split for coloring
jam  <- thedata[thedata$Population == "Jamaica", ]
mal  <- thedata[thedata$Population == "Malawi", ]
yri  <- thedata[thedata$Population == "YRI", ]
msl  <- thedata[thedata$Population == "MSL", ]
ceu  <- thedata[thedata$Population == "CEU", ]
lwk  <- thedata[thedata$Population == "LWK", ]
eafr <- thedata[thedata$Population == "EastAfrica", ]

colorz <- c("red", "orange", "yellow", "green", "skyblue", "blue", "purple")

png(MDS_PNG, width = 1600, height = 1600, res = 200)

plot(jam$C1, jam$C2,
     pch = 20, col = colorz[1],
     xlim = range(thedata$C1, na.rm = TRUE),
     ylim = range(thedata$C2, na.rm = TRUE),
     xlab = "C1", ylab = "C2")

points(mal$C1,  mal$C2,  pch = 20, col = colorz[2])
points(yri$C1,  yri$C2,  pch = 20, col = colorz[3])
points(msl$C1,  msl$C2,  pch = 20, col = colorz[4])
points(ceu$C1,  ceu$C2,  pch = 20, col = colorz[5])
points(lwk$C1,  lwk$C2,  pch = 20, col = colorz[6])
points(eafr$C1, eafr$C2, pch = 20, col = colorz[7])

legend("topright",
       legend = c("Jamaica", "Malawi", "YRI", "MSL", "CEU", "LWK", "East Africa"),
       pch = 20, col = colorz[1:7])

dev.off()
