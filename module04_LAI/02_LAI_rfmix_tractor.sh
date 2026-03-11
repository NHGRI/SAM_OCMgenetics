#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 02_LAI_rfmix_tractor.sh
#
# Purpose:
#   1) Phase study and reference per chromosome with SHAPEIT
#   2) Convert phased output to VCF
#   3) Run RFMix for Local Ancestry Inference
#   4) Run Tractor (ExtractTracts + RunTractor)
#
# Assumptions:
#   - You already ran 01_dataPrep_qc_merge_subset.sh and have:
#       <<STUDY_SUB_PREFIX>>_chr{1..22}.bed/bim/fam
#       <<REF_SUB_PREFIX>>_chr{1..22}.bed/bim/fam
#   - Chromosome-specific genetic maps exist for SHAPEIT
###############################################################################

# ----------------------------- User parameters ------------------------------
# Per-chr PLINK prefixes (from script 01)
STUDY_SUB_PREFIX="<<STUDY_SUB_PREFIX>>"   # e.g., merge_split_studyONLY_subsetREGIONS
REF_SUB_PREFIX="<<REF_SUB_PREFIX>>"       # e.g., merge_split_refONLY_subsetREGIONS

# Tools (use module system or absolute paths)
SHAPEIT_BIN="<<SHAPEIT_BIN>>"             # e.g., /path/to/shapeit
RFMIX_BIN="rfmix"                         # if in PATH
TRACTOR_DIR="<<TRACTOR_DIR>>"             # e.g., /path/to/Tractor

# Maps
SHAPEIT_MAP_DIR="<<SHAPEIT_MAP_DIR>>"     # contains genetic_map_chr${chr}_combined_b37.txt
RFMIX_MAPFILE="<<RFMIX_MAPFILE>>"         # e.g., mapfile-for-rfmix.txt (sample/pop map)
RFMIX_GENMAP="<<RFMIX_GENMAP>>"           # e.g., 1000GP3_geneticmap_forRFMIX_ALLchr.txt

# Tractor phenotype file
TRACTOR_PHE="<<TRACTOR_PHE>>"             # e.g., pheno_for_tractor.txt

# Computation
THREADS="8"
NUM_ANCS="3"
TRACTOR_METHOD="logistic"

# Output directory
OUTDIR="<<OUTDIR>>"                       # e.g., results/lai
# ---------------------------------------------------------------------------

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

###############################################################################
# 0) Modules (optional; keep as comments)
###############################################################################
module load shapeit/2.r904
module load rfmix
module load python
module load plink/1.9

###############################################################################
# 1) Phasing per chromosome (study and reference)
###############################################################################
echo "[INFO] Phasing study + reference with SHAPEIT..."
for chr in {1..22}; do
  MAP="${SHAPEIT_MAP_DIR}/genetic_map_chr${chr}_combined_b37.txt"

  "${SHAPEIT_BIN}" \
    --input-bed "${STUDY_SUB_PREFIX}_chr${chr}.bed" \
               "${STUDY_SUB_PREFIX}_chr${chr}.bim" \
               "${STUDY_SUB_PREFIX}_chr${chr}.fam" \
    --input-map "${MAP}" \
    --output-max "${STUDY_SUB_PREFIX}_chr${chr}_phased.haps" \
                 "${STUDY_SUB_PREFIX}_chr${chr}_phased.sample" \
    --thread "${THREADS}" \
    --force

  "${SHAPEIT_BIN}" \
    --input-bed "${REF_SUB_PREFIX}_chr${chr}.bed" \
               "${REF_SUB_PREFIX}_chr${chr}.bim" \
               "${REF_SUB_PREFIX}_chr${chr}.fam" \
    --input-map "${MAP}" \
    --output-max "${REF_SUB_PREFIX}_chr${chr}_phased.haps" \
                 "${REF_SUB_PREFIX}_chr${chr}_phased.sample" \
    --thread "${THREADS}" \
    --force
done

###############################################################################
# 2) Convert phased output to VCF
###############################################################################
echo "[INFO] Converting phased haps/sample to VCF..."
for chr in {1..22}; do
  "${SHAPEIT_BIN}" -convert \
    --input-haps "${STUDY_SUB_PREFIX}_chr${chr}_phased" \
    --output-vcf "${STUDY_SUB_PREFIX}_chr${chr}_phased.vcf"

  "${SHAPEIT_BIN}" -convert \
    --input-haps "${REF_SUB_PREFIX}_chr${chr}_phased" \
    --output-vcf "${REF_SUB_PREFIX}_chr${chr}_phased.vcf"
done

###############################################################################
# 3) RFMix (local ancestry inference)
###############################################################################
echo "[INFO] Running RFMix..."
for chr in {1..22}; do
  "${RFMIX_BIN}" \
    -f "${STUDY_SUB_PREFIX}_chr${chr}_phased.vcf" \
    -r "${REF_SUB_PREFIX}_chr${chr}_phased.vcf" \
    -m "${RFMIX_MAPFILE}" \
    -g "${RFMIX_GENMAP}" \
    -o "rfmix_${STUDY_SUB_PREFIX}_chr${chr}" \
    --chromosome="${chr}"
done

###############################################################################
# 4) Tractor
###############################################################################
echo "[INFO] Running Tractor ExtractTracts + RunTractor..."
for chr in {1..22}; do
  python "${TRACTOR_DIR}/ExtractTracts.py" \
    --msp "rfmix_${STUDY_SUB_PREFIX}_chr${chr}" \
    --vcf "${STUDY_SUB_PREFIX}_chr${chr}_phased" \
    --num-ancs "${NUM_ANCS}"

  python "${TRACTOR_DIR}/RunTractor.py" \
    --hapdose "${STUDY_SUB_PREFIX}_chr${chr}_phased" \
    --phe "${TRACTOR_PHE}" \
    --method "${TRACTOR_METHOD}" \
    --out "tractor_${STUDY_SUB_PREFIX}_chr${chr}.txt"
done

echo "[INFO] Done: LAI + Tractor per-chr outputs written to ${OUTDIR}"
