#!/usr/bin/env bash
#SBATCH --job-name=gemma_gwas
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cjperkins1@arizona.edu

# Load GEMMA
source /home/u28/cjperkins1/.bashrc
conda activate gemma

# Directories
PLINK_DIR="/xdisk/yadegari/cjperkins1/snp_gwas/data"
WORK_DIR="/xdisk/yadegari/cjperkins1/snp_gwas"
COVAR_DIR="${WORK_DIR}/covariates"
OUTPUT_DIR="${WORK_DIR}/results"

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

echo "========================================="
echo "GEMMA GWAS with PCA Covariates"
echo "========================================="
echo "Start time: $(date)"
echo ""
echo "Using files:"
echo "  PLINK: ${PLINK_DIR}/concordant_gwas.{bed,bim,fam}"
echo "  Phenotypes: ${COVAR_DIR}/gemma_phenotypes.txt"
echo "  Covariates: ${COVAR_DIR}/gemma_covariates_3pcs.txt"
echo ""

# ============================================================================
# Step 1: Calculate Relatedness Matrix (needed for LMM)
# ============================================================================

echo "Step 1: Calculating centered relatedness matrix..."
echo "This may take 30-60 minutes for 3M SNPs..."

gemma -bfile ${PLINK_DIR}/concordant_gwas_no_ref \
  -p ${COVAR_DIR}/gemma_phenotypes.txt \
  -n 1 \
  -gk 1 \
  -o relatedness_matrix

echo "Relatedness matrix saved: output/relatedness_matrix.cXX.txt"
echo ""

# ============================================================================
# Step 2: Run LMM GWAS with covariates for all phenotypes
# ============================================================================

# Count phenotypes (columns - 2 for FID/IID)
N_PHENO=$(head -1 ${COVAR_DIR}/gemma_phenotypes.txt | awk '{print NF-2}')
echo "Number of phenotypes detected: ${N_PHENO}"
echo ""

# Run GWAS for each phenotype
for i in $(seq 1 ${N_PHENO}); do
  echo "========================================="
  echo "Running GWAS for phenotype ${i}..."
  echo "Time: $(date)"
  echo "========================================="
  
  gemma -bfile ${PLINK_DIR}/concordant_gwas_no_ref \
    -k output/relatedness_matrix.cXX.txt \
    -lmm 4 \
    -p ${COVAR_DIR}/gemma_phenotypes.txt \
    -n ${i} \
    -o gwas_pheno${i}_lmm
  
  echo "Results saved: output/gwas_pheno${i}_lmm.assoc.txt"
  echo ""
done

echo "========================================="
echo "GWAS Complete!"
echo "End time: $(date)"
echo "========================================="
echo ""
echo "Output files in: ${OUTPUT_DIR}/output/"
echo ""
ls -lh ${OUTPUT_DIR}/output/*.assoc.txt
echo ""
echo "Next steps:"
echo "1. Download results to local machine"
echo "2. Run visualization and analysis in R"
echo "3. Identify significant SNPs (p < 1e-6)"
