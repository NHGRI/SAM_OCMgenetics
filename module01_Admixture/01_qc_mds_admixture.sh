#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Script: 01_qc_mds_admixture.sh
#
# Description:
#   - Apply MAF filter
#   - LD prune
#   - Remove ambiguous A/T and G/C SNPs
#   - Compute MDS (6 dimensions)
#   - Run ADMIXTURE for a range of K values with CV
#
# Inputs (placeholders; replace with your real paths locally):
#   <<DATA_PREFIX>>            PLINK bed/bim/fam prefix (no extension)
#
# Outputs:
#   A series of PLINK datasets + MDS files + ADMIXTURE logs/Q files
###############################################################################

# ----------------------------- User parameters ------------------------------
DATA_PREFIX="<<DATA_PREFIX>>"     # e.g., /path/to/merged_dataset_round2
MAF_THRESHOLD="0.10"
PRUNE_WINDOW="50"
PRUNE_STEP="10"
PRUNE_R2="0.1"
MDS_DIMS="6"
ADMIX_THREADS="4"
K_MIN="2"
K_MAX="6"
# ---------------------------------------------------------------------------

echo "[INFO] Starting QC pipeline..."
echo "[INFO] DATA_PREFIX=${DATA_PREFIX}"

###############################################################################
# Step 1) Remove SNPs with low MAF
###############################################################################
echo "[INFO] Step 1: MAF filter (MAF >= ${MAF_THRESHOLD})"
module load plink/1.9
plink --bfile "${DATA_PREFIX}" \
  --maf "${MAF_THRESHOLD}" \
  --make-bed \
  --out "${DATA_PREFIX}_mafQC"

###############################################################################
# Step 2) LD prune
###############################################################################
echo "[INFO] Step 2: LD pruning (window=${PRUNE_WINDOW}, step=${PRUNE_STEP}, r2=${PRUNE_R2})"
plink --bfile "${DATA_PREFIX}_mafQC" \
  --indep-pairwise 50 10 0.1 \
  --out "${DATA_PREFIX}_mafQC_pruning"

plink --bfile "${DATA_PREFIX}_mafQC" \
  --extract "${DATA_PREFIX}_mafQC_pruning.prune.in" \
  --make-bed \
  --out "${DATA_PREFIX}_mafQC_pruned"

###############################################################################
# Step 3) Exclude ambiguous A/T and G/C SNPs (strand ambiguous)
###############################################################################
echo "[INFO] Step 3: Excluding A/T and G/C SNPs"
awk '{
  # BIM columns: 1=CHR 2=SNP 3=CM 4=BP 5=A1 6=A2
  if (
    ($5=="A" && $6=="T") || ($5=="T" && $6=="A") ||
    ($5=="C" && $6=="G") || ($5=="G" && $6=="C")
  ) print $2
}' "${DATA_PREFIX}_mafQC_pruned.bim" > "${DATA_PREFIX}_mafQC_pruned_ATGC_snps.txt"

plink --bfile "${DATA_PREFIX}_mafQC_pruned" \
  --exclude "${DATA_PREFIX}_mafQC_pruned_ATGC_snps.txt" \
  --make-bed \
  --out "${DATA_PREFIX}_mafQC_pruned_noATGC"

###############################################################################
# Step 4) MDS (via --genome + --cluster --mds-plot)
###############################################################################
echo "[INFO] Step 4: MDS computation (${MDS_DIMS} dimensions)"
plink --bfile "${DATA_PREFIX}_mafQC_pruned_noATGC" \
  --genome \
  --out "${DATA_PREFIX}_mafQC_pruned_noATGC_IBD"

plink --bfile "${DATA_PREFIX}_mafQC_pruned_noATGC" \
  --read-genome "${DATA_PREFIX}_mafQC_pruned_noATGC_IBD.genome" \
  --cluster \
  --mds-plot "${MDS_DIMS}" \
  --out "${DATA_PREFIX}_mafQC_pruned_noATGC_MDS"

###############################################################################
# Step 5) ADMIXTURE K sweep
###############################################################################
echo "[INFO] Step 5: ADMIXTURE K sweep (K=${K_MIN}..${K_MAX})"
module load admixture
for K in $(seq "${K_MIN}" "${K_MAX}"); do
  echo "[INFO] Running ADMIXTURE K=${K}"
  admixture -j"${ADMIX_THREADS}" --cv \
    "${DATA_PREFIX}_mafQC_pruned_noATGC.bed" "${K}" \
    | tee "${DATA_PREFIX}_mafQC_pruned_noATGC_admix_logK${K}.out"
done

echo "[INFO] Done."
