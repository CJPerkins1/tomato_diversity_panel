#!/usr/bin/env bash
#SBATCH --job-name=ld_stratified
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cjperkins1@arizona.edu

# Load PLINK
module load plink/1.9

# Directories
PLINK_DIR="/xdisk/yadegari/cjperkins1/LD/results"
OUTPUT_DIR="/xdisk/yadegari/cjperkins1/LD/stratified_results"
SAMPLE_LISTS="/xdisk/yadegari/cjperkins1/LD/sample_lists"

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

echo "========================================="
echo "Stratified LD Analysis by Species"
echo "========================================="
echo ""
echo "DEBUG: Checking directories..."
echo "PLINK_DIR: ${PLINK_DIR}"
echo "OUTPUT_DIR: ${OUTPUT_DIR}"
echo "SAMPLE_LISTS: ${SAMPLE_LISTS}"
echo ""
echo "Files in SAMPLE_LISTS directory:"
ls -lh ${SAMPLE_LISTS}/
echo ""
echo "File check results:"
test -f "${SAMPLE_LISTS}/lycopersicum_samples.txt" && echo "  ✓ lycopersicum_samples.txt exists" || echo "  ✗ lycopersicum_samples.txt NOT FOUND"
test -f "${SAMPLE_LISTS}/pimpinellifolium_samples.txt" && echo "  ✓ pimpinellifolium_samples.txt exists" || echo "  ✗ pimpinellifolium_samples.txt NOT FOUND"
test -f "${SAMPLE_LISTS}/cultivated_samples.txt" && echo "  ✓ cultivated_samples.txt exists" || echo "  ✗ cultivated_samples.txt NOT FOUND"
echo ""

# Get chromosomes
CHROMOSOMES=$(awk '{print $1}' ${PLINK_DIR}/concordant_popgen.bim | sort -u)
echo "Chromosomes to process:"
echo "${CHROMOSOMES}"
echo ""

# Process lycopersicum
GROUP="lycopersicum"
echo "========================================="
echo "Processing group: ${GROUP}"
echo "========================================="

SAMPLE_FILE="${SAMPLE_LISTS}/${GROUP}_samples.txt"
echo "Looking for: ${SAMPLE_FILE}"

if [ -f "${SAMPLE_FILE}" ]; then
  N_SAMPLES=$(wc -l < ${SAMPLE_FILE})
  echo "Found! Number of samples: ${N_SAMPLES}"
  echo "First 3 samples:"
  head -3 ${SAMPLE_FILE}
  echo ""
  
  for CHR in ${CHROMOSOMES}; do
    echo "  Processing ${GROUP} - ${CHR}..."
    plink --bfile ${PLINK_DIR}/concordant_popgen \
      --allow-extra-chr \
      --keep ${SAMPLE_FILE} \
      --chr ${CHR} \
      --r2 gz dprime \
      --ld-window-kb 1000 \
      --ld-window 99999 \
      --ld-window-r2 0 \
      --thin 0.1 \
      --out ld_${GROUP}_${CHR}
  done
else
  echo "ERROR: File not found at: ${SAMPLE_FILE}"
fi

echo ""

# Process pimpinellifolium
GROUP="pimpinellifolium"
echo "========================================="
echo "Processing group: ${GROUP}"
echo "========================================="

SAMPLE_FILE="${SAMPLE_LISTS}/${GROUP}_samples.txt"
echo "Looking for: ${SAMPLE_FILE}"

if [ -f "${SAMPLE_FILE}" ]; then
  N_SAMPLES=$(wc -l < ${SAMPLE_FILE})
  echo "Found! Number of samples: ${N_SAMPLES}"
  echo "First 3 samples:"
  head -3 ${SAMPLE_FILE}
  echo ""
  
  for CHR in ${CHROMOSOMES}; do
    echo "  Processing ${GROUP} - ${CHR}..."
    plink --bfile ${PLINK_DIR}/concordant_popgen \
      --allow-extra-chr \
      --keep ${SAMPLE_FILE} \
      --chr ${CHR} \
      --r2 gz dprime \
      --ld-window-kb 1000 \
      --ld-window 99999 \
      --ld-window-r2 0 \
      --thin 0.1 \
      --out ld_${GROUP}_${CHR}
  done
else
  echo "ERROR: File not found at: ${SAMPLE_FILE}"
fi

echo ""

# Process cultivated
GROUP="cultivated"
echo "========================================="
echo "Processing group: ${GROUP}"
echo "========================================="

SAMPLE_FILE="${SAMPLE_LISTS}/${GROUP}_samples.txt"
echo "Looking for: ${SAMPLE_FILE}"

if [ -f "${SAMPLE_FILE}" ]; then
  N_SAMPLES=$(wc -l < ${SAMPLE_FILE})
  echo "Found! Number of samples: ${N_SAMPLES}"
  echo "First 3 samples:"
  head -3 ${SAMPLE_FILE}
  echo ""
  
  for CHR in ${CHROMOSOMES}; do
    echo "  Processing ${GROUP} - ${CHR}..."
    plink --bfile ${PLINK_DIR}/concordant_popgen \
      --allow-extra-chr \
      --keep ${SAMPLE_FILE} \
      --chr ${CHR} \
      --r2 gz dprime \
      --ld-window-kb 1000 \
      --ld-window 99999 \
      --ld-window-r2 0 \
      --thin 0.1 \
      --out ld_${GROUP}_${CHR}
  done
else
  echo "ERROR: File not found at: ${SAMPLE_FILE}"
fi

echo ""
echo "========================================="
echo "All groups processed!"
echo "========================================="
