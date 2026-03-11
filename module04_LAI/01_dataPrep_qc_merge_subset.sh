#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 01_dataPrep_qc_merge_subset.sh
#
# Purpose:
#   1) QC reference (e.g., 1000G subset) and study data
#   2) Normalize variant IDs to CHR:POS to improve matching
#   3) Remove duplicate variants
#   4) Merge study + reference (with flip/exclude logic)
#   5) Post-merge QC + remove ambiguous A/T and G/C SNPs
#   6) Split merged data into study-only and ref-only
#   7) Subset to target regions (e.g., OCM) and split by chromosome
#
# Notes:
#   - All filenames/paths are placeholders. Replace locally or via env vars.
###############################################################################

# ----------------------------- User parameters ------------------------------
# Input datasets (PLINK bed/bim/fam prefix, no extension)
REF_BFILE="<<REF_BFILE>>"          # e.g., 1000G_mslceulwk_only
STUDY_BFILE="<<STUDY_BFILE>>"      # e.g., study_QC_autosomes_phenosex

# Optional phenotype/sex update (set to empty string to skip)
PHENO_FILE="<<PHENO_FILE>>"        # e.g., pheno.txt
SEX_UPDATE_FILE="<<SEX_UPDATE_FILE>>"  # e.g., sex_update.txt

# Keep files for splitting merged dataset
KEEP_STUDY_FAM="<<KEEP_STUDY_FAM>>"     # e.g., original study .fam or keep list
KEEP_REF_LIST="<<KEEP_REF_LIST>>"       # e.g., keep_mslceulwk.txt

# Region extract file (plink --extract range)
REGION_RANGE_FILE="<<REGION_RANGE_FILE>>"   # e.g., extract_OCM_regions.txt

# QC thresholds
#REF_MAF="0.05"
#REF_GENO="0.05"
#STUDY_MAF="0.05"
#POSTMERGE_MAF="0.05"
#POSTMERGE_GENO="0.05"

# Output prefix root
OUTROOT="<<OUTROOT>>"              # e.g., results/tractor_pipeline
# ---------------------------------------------------------------------------

mkdir -p "${OUTROOT}"
cd "${OUTROOT}"

###############################################################################
# 0) Load modules (optional; keep as comments for portability)
###############################################################################
# module load plink/1.9
# module load plink/2.0
# module load samtools
# module load R

###############################################################################
# 1) QC reference dataset (1000G subset)
###############################################################################
echo "[INFO] QC reference dataset..."
REF_QC="${REF_BFILE}_mafQC005genoQC005"

module load plink/1.9

plink --bfile "${REF_BFILE}" \
  --maf 0.05 \
  --geno 0.05 \
  --make-bed \
  --out "${REF_QC}"

###############################################################################
# 2) QC study dataset + optionally add phenotype and update sex
###############################################################################
echo "[INFO] QC study dataset..."
STUDY_QC="${STUDY_BFILE}_mafQC005"

plink --bfile "${STUDY_BFILE}" \
  --chr 1-22 \
  --pheno "${PHENO_FILE}" \
  --update-sex "${SEX_UPDATE_FILE}" \
  --maf 0.05 \
  --make-bed \
  --out "${STUDY_QC}"


###############################################################################
# 3) Normalize variant IDs to CHR:POS using PLINK2
###############################################################################
echo "[INFO] Normalizing variant IDs to CHR:POS..."
module load plink/2.0
plink2 --bfile "${REF_QC}"   --set-all-var-ids @:# --make-bed --out "${REF_QC}_chrpos"
plink2 --bfile "${STUDY_QC}" --set-all-var-ids @:# --make-bed --out "${STUDY_QC}_chrpos"

###############################################################################
# 4) Remove duplicate variants (force-first)
###############################################################################
echo "[INFO] Removing duplicate variants..."
plink2 --bfile "${STUDY_QC}_chrpos" --rm-dup force-first --make-bed --out "${STUDY_QC}_chrpos_rmdups"
plink2 --bfile "${REF_QC}_chrpos"   --rm-dup force-first --make-bed --out "${REF_QC}_chrpos_rmdups"

###############################################################################
# 5) Merge study and reference
#    - Find common CHR SNP BP triplets
#    - Extract matches
#    - Merge with flip/exclude logic
###############################################################################
echo "[INFO] Preparing for merge..."

STUDY_MERGE_IN="${STUDY_QC}_chrpos_rmdups"
REF_MERGE_IN="${REF_QC}_chrpos_rmdups"
MERGE_PREFIX="<<MERGE_PREFIX>>"  # e.g., merged_study_ref

# Build comm files: "CHR SNP BP" then comm -12
awk '{print $1,$2,$4}' "${STUDY_MERGE_IN}.bim" | sort > study_comm.txt
awk '{print $1,$2,$4}' "${REF_MERGE_IN}.bim"   | sort > ref_comm.txt

comm -12 study_comm.txt ref_comm.txt | awk '{print $2}' > MatchSNPs.txt

echo "[INFO] Extracting common SNPs..."
plink --bfile "${STUDY_MERGE_IN}" --extract MatchSNPs.txt --make-bed --out study_MatchSNPs
plink --bfile "${REF_MERGE_IN}"   --extract MatchSNPs.txt --make-bed --out ref_MatchSNPs

echo "[INFO] Merge round 1..."
plink --bfile study_MatchSNPs \
  --bmerge ref_MatchSNPs.bed ref_MatchSNPs.bim ref_MatchSNPs.fam \
  --make-bed \
  --out "${MERGE_PREFIX}_round1" || true

# When merge creates a missnp file, attempt flipping those SNPs in study
if [[ -f "${MERGE_PREFIX}_round1-merge.missnp" ]]; then
  echo "[INFO] Flipping SNPs listed in round1 missnp and retry merge..."
  plink --bfile study_MatchSNPs \
    --flip "${MERGE_PREFIX}_round1-merge.missnp" \
    --make-bed \
    --out study_MatchSNPs_flipped

  echo "[INFO] Merge round 2..."
  plink --bfile study_MatchSNPs_flipped \
    --bmerge ref_MatchSNPs.bed ref_MatchSNPs.bim ref_MatchSNPs.fam \
    --make-bed \
    --out "${MERGE_PREFIX}_round2" || true

  # If still missnp, exclude from both and merge again
  if [[ -f "${MERGE_PREFIX}_round2-merge.missnp" ]]; then
    echo "[INFO] Excluding remaining mismatched SNPs and merging round 3..."
    plink --bfile study_MatchSNPs_flipped \
      --exclude "${MERGE_PREFIX}_round2-merge.missnp" \
      --make-bed \
      --out study_MatchSNPs_flipped_excl

    plink --bfile ref_MatchSNPs \
      --exclude "${MERGE_PREFIX}_round2-merge.missnp" \
      --make-bed \
      --out ref_MatchSNPs_excl

    plink --bfile study_MatchSNPs_flipped_excl \
      --bmerge ref_MatchSNPs_excl.bed ref_MatchSNPs_excl.bim ref_MatchSNPs_excl.fam \
      --make-bed \
      --out "${MERGE_PREFIX}_final"
  else
    mv "${MERGE_PREFIX}_round2.bed" "${MERGE_PREFIX}_final.bed"
    mv "${MERGE_PREFIX}_round2.bim" "${MERGE_PREFIX}_final.bim"
    mv "${MERGE_PREFIX}_round2.fam" "${MERGE_PREFIX}_final.fam"
  fi
else
  mv "${MERGE_PREFIX}_round1.bed" "${MERGE_PREFIX}_final.bed"
  mv "${MERGE_PREFIX}_round1.bim" "${MERGE_PREFIX}_final.bim"
  mv "${MERGE_PREFIX}_round1.fam" "${MERGE_PREFIX}_final.fam"
fi

###############################################################################
# 6) Post-merge QC + remove ambiguous A/T and G/C SNPs
###############################################################################
echo "[INFO] Post-merge QC..."
POSTQC="${MERGE_PREFIX}_final_genomafQC"
plink --bfile "${MERGE_PREFIX}_final" \
  --geno "${POSTMERGE_GENO}" \
  --maf "${POSTMERGE_MAF}" \
  --make-bed \
  --out "${POSTQC}"

echo "[INFO] Removing A/T and G/C SNPs..."
awk '{
  if (
    ($5=="A" && $6=="T") || ($5=="T" && $6=="A") ||
    ($5=="C" && $6=="G") || ($5=="G" && $6=="C")
  ) print $2
}' "${POSTQC}.bim" > "${POSTQC}_ATGC_snps.txt"

plink --bfile "${POSTQC}" \
  --exclude "${POSTQC}_ATGC_snps.txt" \
  --make-bed \
  --out "${POSTQC}_noATGC"

###############################################################################
# 7) Split merged into study-only and reference-only
###############################################################################
echo "[INFO] Splitting merged dataset..."
MERGE_SPLIT_PREFIX="<<MERGE_SPLIT_PREFIX>>"  # e.g., merge_split

plink --bfile "${POSTQC}_noATGC" \
  --keep "${KEEP_STUDY_FAM}" \
  --make-bed \
  --out "${MERGE_SPLIT_PREFIX}_studyONLY"

plink --bfile "${POSTQC}_noATGC" \
  --keep "${KEEP_REF_LIST}" \
  --make-bed \
  --out "${MERGE_SPLIT_PREFIX}_refONLY"

###############################################################################
# 8) Subset to regions (e.g., OCM) and split by chromosome
###############################################################################
echo "[INFO] Subsetting to target regions and splitting by chromosome..."
for SOURCE in studyONLY refONLY; do
  INP="${MERGE_SPLIT_PREFIX}_${SOURCE}"
  SUB="${INP}_subsetREGIONS"

  plink --bfile "${INP}" \
    --extract range "${REGION_RANGE_FILE}" \
    --make-bed \
    --out "${SUB}"

  for chr in {1..22}; do
    plink --bfile "${SUB}" \
      --chr "${chr}" \
      --make-bed \
      --out "${SUB}_chr${chr}"
  done
done

echo "[INFO] Done: data prep + QC + merge + subset + per-chr splits."
