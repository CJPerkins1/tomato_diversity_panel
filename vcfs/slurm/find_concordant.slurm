#!/usr/bin/env bash
#SBATCH --job-name=concordant_snps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=48:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cjperkins1@arizona.edu

# Load bcftools
module load bcftools

# Paths
ONT_NORM="/xdisk/yadegari/cjperkins1/dp_vcfs/ont_raw_vcfs_2026-04-20/normalized"
ILLUMINA_NORM="/xdisk/yadegari/cjperkins1/dp_vcfs/illumina_raw_vcfs_2026-04-16/normalized"
OUTDIR="/xdisk/yadegari/cjperkins1/dp_vcfs/concordant_snps"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

echo "========================================"
echo "Finding Concordant SNPs"
echo "========================================"
echo "Start time: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo ""

# ============================================================================
# STEP 1: Concatenate chromosomes per dataset
# ============================================================================
echo "[Step 1/4] Concatenating ONT chromosomes..."
ls ${ONT_NORM}/ont_SL4.0ch*.snps.vcf.gz | sort -V > ont_snps_files.txt
bcftools concat -a -f ont_snps_files.txt -Oz -o ont_genome.vcf.gz
bcftools index --threads ${SLURM_CPUS_PER_TASK} ont_genome.vcf.gz
ONT_COUNT=$(bcftools view -H ont_genome.vcf.gz | wc -l)
echo "  ✓ ONT: ${ONT_COUNT} SNPs"
echo ""

echo "[Step 2/4] Concatenating Illumina chromosomes..."
ls ${ILLUMINA_NORM}/illumina_SL4.0ch*.snps.vcf.gz | sort -V > illumina_snps_files.txt
bcftools concat -a -f illumina_snps_files.txt -Oz -o illumina_genome.vcf.gz
bcftools index --threads ${SLURM_CPUS_PER_TASK} illumina_genome.vcf.gz
ILLUMINA_COUNT=$(bcftools view -H illumina_genome.vcf.gz | wc -l)
echo "  ✓ Illumina: ${ILLUMINA_COUNT} SNPs"
echo ""

# ============================================================================
# STEP 2: Find concordant positions (intersection)
# ============================================================================
echo "[Step 3/4] Finding concordant SNPs (present in BOTH technologies)..."
bcftools isec \
    -n=2 \
    -w1 \
    -Oz -o illumina_concordant.vcf.gz \
    illumina_genome.vcf.gz \
    ont_genome.vcf.gz

bcftools isec \
    -n=2 \
    -w2 \
    -Oz -o ont_concordant.vcf.gz \
    illumina_genome.vcf.gz \
    ont_genome.vcf.gz

bcftools index --threads ${SLURM_CPUS_PER_TASK} illumina_concordant.vcf.gz
bcftools index --threads ${SLURM_CPUS_PER_TASK} ont_concordant.vcf.gz

CONCORDANT_COUNT=$(bcftools view -H illumina_concordant.vcf.gz | wc -l)
echo "  ✓ Concordant SNPs: ${CONCORDANT_COUNT}"
echo ""

# ============================================================================
# STEP 3: Merge all samples at concordant sites
# ============================================================================
echo "[Step 4/4] Merging all samples (146 Illumina + 63 ONT = 209 total)..."
bcftools merge \
    --threads ${SLURM_CPUS_PER_TASK} \
    --force-samples \
    --missing-to-ref \
    -Oz -o merged_concordant_snps.vcf.gz \
    illumina_concordant.vcf.gz \
    ont_concordant.vcf.gz

bcftools index --threads ${SLURM_CPUS_PER_TASK} merged_concordant_snps.vcf.gz

SAMPLE_COUNT=$(bcftools query -l merged_concordant_snps.vcf.gz | wc -l)
echo "  ✓ Merged: ${SAMPLE_COUNT} samples at ${CONCORDANT_COUNT} SNP sites"
echo ""

# ============================================================================
# Summary Statistics
# ============================================================================
echo "Generating summary statistics..."
{
    echo "=== Concordant SNP Summary ==="
    echo ""
    echo "Date: $(date)"
    echo ""
    echo "Input datasets:"
    echo "  ONT SNPs: ${ONT_COUNT}"
    echo "  Illumina SNPs: ${ILLUMINA_COUNT}"
    echo ""
    echo "Concordant SNPs (intersection): ${CONCORDANT_COUNT}"
    echo "  Concordance rate (ONT): $(awk "BEGIN {printf \"%.2f%%\", (${CONCORDANT_COUNT}/${ONT_COUNT})*100}")"
    echo "  Concordance rate (Illumina): $(awk "BEGIN {printf \"%.2f%%\", (${CONCORDANT_COUNT}/${ILLUMINA_COUNT})*100}")"
    echo ""
    echo "Sample counts:"
    echo "  ONT samples: $(bcftools query -l ont_genome.vcf.gz | wc -l)"
    echo "  Illumina samples: $(bcftools query -l illumina_genome.vcf.gz | wc -l)"
    echo "  Total merged samples: ${SAMPLE_COUNT}"
} > concordant_summary.txt

cat concordant_summary.txt

echo ""
echo "========================================"
echo "Complete!"
echo "========================================"
echo "End time: $(date)"
echo ""
echo "Output files in: ${OUTDIR}"
echo "  - ont_genome.vcf.gz              (all ONT SNPs)"
echo "  - illumina_genome.vcf.gz         (all Illumina SNPs)"
echo "  - illumina_concordant.vcf.gz     (Illumina SNPs also in ONT)"
echo "  - ont_concordant.vcf.gz          (ONT SNPs also in Illumina)"
echo "  - merged_concordant_snps.vcf.gz  (209 samples at concordant sites)"
echo "  - concordant_summary.txt         (statistics)"
echo ""
