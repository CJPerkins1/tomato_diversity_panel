#!/usr/bin/env bash
#SBATCH --job-name=ld_prune
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

#=========================================================================
# LD Pruning Pipeline
# Input: Filtered concordant SNPs (popgen_filtered)
# Output: LD-pruned datasets for population genomics
#=========================================================================

set -eo pipefail

# Load required modules
module load plink/2.0

# Input/Output paths
WORKDIR="/xdisk/yadegari/cjperkins1/dp_vcfs/concordant_snps"
INPUT_VCF="${WORKDIR}/popgen_filtered/concordant_popgen.filtered.vcf.gz"
OUTDIR="${WORKDIR}/ld_pruned"
PREFIX="concordant_popgen"

# Create output directory
mkdir -p ${OUTDIR}

# Number of threads
THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "========================================================================="
echo "LD Pruning Pipeline - PLINK 2.0"
echo "========================================================================="
echo ""
echo "Job ID: ${SLURM_JOB_ID:-'LOCAL'}"
echo "Running on: $(hostname)"
echo "Cores: ${THREADS}"
echo ""
echo "Input VCF: ${INPUT_VCF}"
echo "Output directory: ${OUTDIR}"
echo ""
echo "========================================================================="

# Step 1: Convert VCF to PLINK binary format
echo ""
echo "[Step 1/5] Converting VCF to PLINK binary format..."
echo "Start time: $(date)"

plink2 \
  --vcf ${INPUT_VCF} \
  --make-bed \
  --out ${OUTDIR}/${PREFIX} \
  --threads ${THREADS} \
  --allow-extra-chr

echo "  ✓ PLINK binary files created"
echo "End time: $(date)"

# Get initial SNP count
INITIAL_SNPS=$(wc -l < ${OUTDIR}/${PREFIX}.bim)
echo "  ✓ Input SNPs: $(printf "%'d" ${INITIAL_SNPS})"

# Step 2: Standard LD pruning (for PCA, ADMIXTURE)
echo ""
echo "[Step 2/5] LD pruning - Standard (r² < 0.2)..."
echo "Start time: $(date)"
echo ""
echo "Parameters:"
echo "  - Window: 50 SNPs"
echo "  - Step: 5 SNPs"
echo "  - r² threshold: 0.2"
echo ""

plink2 \
  --bfile ${OUTDIR}/${PREFIX} \
  --indep-pairwise 50 5 0.2 \
  --out ${OUTDIR}/${PREFIX}_ld02 \
  --threads ${THREADS} \
  --allow-extra-chr

# Extract pruned SNPs
plink2 \
  --bfile ${OUTDIR}/${PREFIX} \
  --extract ${OUTDIR}/${PREFIX}_ld02.prune.in \
  --make-bed \
  --out ${OUTDIR}/${PREFIX}_ld02_pruned \
  --threads ${THREADS} \
  --allow-extra-chr

# Also export as VCF
plink2 \
  --bfile ${OUTDIR}/${PREFIX}_ld02_pruned \
  --export vcf bgz \
  --out ${OUTDIR}/${PREFIX}_ld02_pruned \
  --threads ${THREADS} \
  --allow-extra-chr

PRUNED_SNPS=$(wc -l < ${OUTDIR}/${PREFIX}_ld02_pruned.bim)
echo ""
echo "  ✓ LD-pruned SNPs (r² < 0.2): $(printf "%'d" ${PRUNED_SNPS})"
echo "  ✓ Removed: $(printf "%'d" $((INITIAL_SNPS - PRUNED_SNPS))) SNPs"
echo "  ✓ Retention: $(echo "scale=2; ${PRUNED_SNPS} * 100 / ${INITIAL_SNPS}" | bc)%"
echo "End time: $(date)"

# Step 3: Strict LD pruning (for maximum independence)
echo ""
echo "[Step 3/5] LD pruning - Strict (r² < 0.1)..."
echo "Start time: $(date)"
echo ""
echo "Parameters:"
echo "  - Window: 100 SNPs"
echo "  - Step: 10 SNPs"
echo "  - r² threshold: 0.1"
echo ""

plink2 \
  --bfile ${OUTDIR}/${PREFIX} \
  --indep-pairwise 100 10 0.1 \
  --out ${OUTDIR}/${PREFIX}_ld01 \
  --threads ${THREADS} \
  --allow-extra-chr

# Extract pruned SNPs
plink2 \
  --bfile ${OUTDIR}/${PREFIX} \
  --extract ${OUTDIR}/${PREFIX}_ld01.prune.in \
  --make-bed \
  --out ${OUTDIR}/${PREFIX}_ld01_pruned \
  --threads ${THREADS} \
  --allow-extra-chr

# Also export as VCF
plink2 \
  --bfile ${OUTDIR}/${PREFIX}_ld01_pruned \
  --export vcf bgz \
  --out ${OUTDIR}/${PREFIX}_ld01_pruned \
  --threads ${THREADS} \
  --allow-extra-chr

STRICT_SNPS=$(wc -l < ${OUTDIR}/${PREFIX}_ld01_pruned.bim)
echo ""
echo "  ✓ Strict LD-pruned SNPs (r² < 0.1): $(printf "%'d" ${STRICT_SNPS})"
echo "  ✓ Removed: $(printf "%'d" $((INITIAL_SNPS - STRICT_SNPS))) SNPs"
echo "  ✓ Retention: $(echo "scale=2; ${STRICT_SNPS} * 100 / ${INITIAL_SNPS}" | bc)%"
echo "End time: $(date)"

# Step 4: Thin for phylogenetics
echo ""
echo "[Step 4/5] Distance-based thinning (1 SNP per 50kb)..."
echo "Start time: $(date)"
echo ""

plink2 \
  --bfile ${OUTDIR}/${PREFIX}_ld02_pruned \
  --bp-space 50000 \
  --make-bed \
  --out ${OUTDIR}/${PREFIX}_thinned_50kb \
  --threads ${THREADS} \
  --allow-extra-chr

# Also export as VCF
plink2 \
  --bfile ${OUTDIR}/${PREFIX}_thinned_50kb \
  --export vcf bgz \
  --out ${OUTDIR}/${PREFIX}_thinned_50kb \
  --threads ${THREADS} \
  --allow-extra-chr

THINNED_SNPS=$(wc -l < ${OUTDIR}/${PREFIX}_thinned_50kb.bim)
echo ""
echo "  ✓ Thinned SNPs (50kb spacing): $(printf "%'d" ${THINNED_SNPS})"
echo "  ✓ Removed: $(printf "%'d" $((PRUNED_SNPS - THINNED_SNPS))) SNPs"
echo "  ✓ Retention from LD-pruned: $(echo "scale=2; ${THINNED_SNPS} * 100 / ${PRUNED_SNPS}" | bc)%"
echo "End time: $(date)"

# Step 5: Generate summary report
echo ""
echo "[Step 5/5] Generating summary report..."
echo "Start time: $(date)"

# Per-chromosome counts for each dataset
echo "  Computing per-chromosome statistics..."

for DATASET in "${PREFIX}_ld02_pruned" "${PREFIX}_ld01_pruned" "${PREFIX}_thinned_50kb"; do
  awk '{print $1}' ${OUTDIR}/${DATASET}.bim | sort | uniq -c | \
    awk -v dataset="${DATASET}" '{printf "    %s - %s: %'"'"'d SNPs\n", dataset, $2, $1}' \
    >> ${OUTDIR}/chr_counts_all.tmp
done

cat > ${OUTDIR}/ld_pruning_report.txt << REPORT
========================================================================
LD PRUNING REPORT
========================================================================

Job ID: ${SLURM_JOB_ID:-'LOCAL'}
Date: $(date)
Input: ${INPUT_VCF}
Output directory: ${OUTDIR}

------------------------------------------------------------------------
PRUNING RESULTS
------------------------------------------------------------------------

Input SNPs (filtered):         $(printf "%'12d" ${INITIAL_SNPS})

Standard LD-pruned (r²<0.2):   $(printf "%'12d" ${PRUNED_SNPS})
  Retention:                    $(echo "scale=2; ${PRUNED_SNPS} * 100 / ${INITIAL_SNPS}" | bc)%
  Use for: PCA, ADMIXTURE, structure analysis

Strict LD-pruned (r²<0.1):     $(printf "%'12d" ${STRICT_SNPS})
  Retention:                    $(echo "scale=2; ${STRICT_SNPS} * 100 / ${INITIAL_SNPS}" | bc)%
  Use for: Maximum independence analyses

Thinned (50kb spacing):        $(printf "%'12d" ${THINNED_SNPS})
  Retention:                    $(echo "scale=2; ${THINNED_SNPS} * 100 / ${INITIAL_SNPS}" | bc)%
  Use for: Phylogenetics (RAxML-ng), TreeMix

------------------------------------------------------------------------
OUTPUT FILES
------------------------------------------------------------------------

Standard LD-pruned (RECOMMENDED for most analyses):
  PLINK binary: ${PREFIX}_ld02_pruned.{bed,bim,fam}
  VCF:          ${PREFIX}_ld02_pruned.vcf.gz
  Pruned SNPs:  ${PREFIX}_ld02.prune.in
  Removed SNPs: ${PREFIX}_ld02.prune.out

Strict LD-pruned:
  PLINK binary: ${PREFIX}_ld01_pruned.{bed,bim,fam}
  VCF:          ${PREFIX}_ld01_pruned.vcf.gz

Thinned for phylogenetics:
  PLINK binary: ${PREFIX}_thinned_50kb.{bed,bim,fam}
  VCF:          ${PREFIX}_thinned_50kb.vcf.gz

Original PLINK files (full dataset):
  PLINK binary: ${PREFIX}.{bed,bim,fam}

------------------------------------------------------------------------
PRUNING PARAMETERS
------------------------------------------------------------------------

Standard LD pruning:
  Window: 50 SNPs
  Step: 5 SNPs  
  r² threshold: 0.2
  (Removes one SNP from pairs with r² > 0.2)

Strict LD pruning:
  Window: 100 SNPs
  Step: 10 SNPs
  r² threshold: 0.1

Distance-based thinning:
  Spacing: 50,000 bp (50kb)
  Applied to: Standard LD-pruned set

------------------------------------------------------------------------
PER-CHROMOSOME COUNTS
------------------------------------------------------------------------
$(cat ${OUTDIR}/chr_counts_all.tmp)

========================================================================
REPORT

rm ${OUTDIR}/chr_counts_all.tmp

echo "  ✓ Report saved"
echo "End time: $(date)"

echo ""
echo "========================================================================="
echo "LD PRUNING COMPLETE!"
echo "========================================================================="
echo ""
echo "Summary:"
echo "  Input SNPs:              $(printf "%'12d" ${INITIAL_SNPS})"
echo "  Standard pruned (r²<0.2): $(printf "%'12d" ${PRUNED_SNPS})"
echo "  Strict pruned (r²<0.1):   $(printf "%'12d" ${STRICT_SNPS})"
echo "  Thinned (50kb):          $(printf "%'12d" ${THINNED_SNPS})"
echo ""
echo "Output files in: ${OUTDIR}/"
echo "Report: ${OUTDIR}/ld_pruning_report.txt"
echo ""
echo "========================================================================="
echo ""
echo "RECOMMENDED NEXT STEPS:"
echo ""
echo "1. PCA Analysis:"
echo "   Use: ${PREFIX}_ld02_pruned.bed"
echo "   Tool: plink2 --bfile ${PREFIX}_ld02_pruned --pca"
echo ""
echo "2. ADMIXTURE Analysis:"
echo "   Use: ${PREFIX}_ld02_pruned.bed"
echo "   Tool: admixture ${PREFIX}_ld02_pruned.bed K"
echo ""
echo "3. Phylogenetic Analysis (RAxML-ng):"
echo "   Use: ${PREFIX}_thinned_50kb.vcf.gz"
echo "   Convert to phylip/fasta format first"
echo ""
echo "========================================================================="
echo ""
echo "Total runtime: $SECONDS seconds"
echo "End time: $(date)"
echo ""

