#!/bin/bash
#SBATCH --job-name=ld_plink
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=ld_plink_%j.out
#SBATCH --error=ld_plink_%j.err

# Load PLINK
module load plink/1.9

VCF_DIR="/xdisk/yadegari/cjperkins1/dp_vcfs/concordant_snps/popgen_filtered"
VCF_FILE="concordant_popgen.filtered.vcf.gz"
OUTPUT_DIR="/xdisk/yadegari/cjperkins1/LD/results"

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

echo "Starting PLINK LD analysis..."
echo "Step 1: Converting VCF to PLINK format..."
plink --vcf ${VCF_DIR}/${VCF_FILE} \
  --double-id \
  --allow-extra-chr \
  --maf 0.05 \
  --geno 0.1 \
  --make-bed \
  --out concordant_popgen

echo ""
echo "Step 2: Calculating LD per chromosome..."
echo "This prevents memory issues by processing one chromosome at a time"
echo ""

# Get list of chromosomes
CHROMOSOMES=$(awk '{print $1}' concordant_popgen.bim | sort -u)

for chr in ${CHROMOSOMES}; do
  echo "========================================"
  echo "Processing chromosome: ${chr}"
  echo "========================================"
  
  plink --bfile concordant_popgen \
    --allow-extra-chr \
    --chr ${chr} \
    --r2 gz dprime \
    --ld-window-kb 1000 \
    --ld-window 99999 \
    --ld-window-r2 0 \
    --thin 0.1 \
    --out ld_${chr}
  
  echo "Chromosome ${chr} complete: ld_${chr}.ld.gz"
  echo ""
done

echo "========================================"
echo "All chromosomes processed successfully!"
echo "========================================"
echo "Output files located in: ${OUTPUT_DIR}"
ls -lh ${OUTPUT_DIR}/ld_*.ld.gz

echo ""
echo "To plot LD decay for a chromosome, run:"
echo "  module load R"
echo "  Rscript ../plot_ld_decay.R ld_chr1.ld.gz"
