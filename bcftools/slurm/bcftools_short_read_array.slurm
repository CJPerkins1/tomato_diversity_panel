#!/usr/bin/env bash
#SBATCH --job-name=tomato_bcftools_short_reads_array
#SBATCH --array=1-153
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=04:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu

set -euo pipefail

# -----------------------
# LOAD MODULES
# -----------------------
module load bcftools/1.19
module load samtools

# -----------------------
# USER SETTINGS
# -----------------------
INPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/aligned_short_reads"
REF="/xdisk/rpalaniv/cjperkins1/dp_alignment/reference/S_lycopersicum_chromosomes.4.00.fa"
OUT_DIR="/xdisk/yadegari/cjperkins1/dp_vcfs/output_short_reads"

mkdir -p ${OUT_DIR}

# Get array of BAM files
BAM_FILES=(${INPUT_DIR}/*.bam)
BAM=${BAM_FILES[$((SLURM_ARRAY_TASK_ID - 1))]}

# Check if BAM file exists
if [ ! -f "${BAM}" ]; then
    echo "ERROR: BAM file not found: ${BAM}"
    exit 1
fi

SAMPLE=$(basename ${BAM} .bam)
OUT_VCF=${OUT_DIR}/${SAMPLE}.vcf.gz

echo "=========================================="
echo "Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Processing: ${SAMPLE}"
echo "Input: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUT_VCF}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "=========================================="

# -----------------------
# Variant Calling (Illumina Paired-End Optimized)
# -----------------------
bcftools mpileup \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f ${REF} \
    -a AD,DP,SP \
    --min-MQ 30 \
    --min-BQ 20 \
    -Ou ${BAM} | \
bcftools call \
    --threads ${SLURM_CPUS_PER_TASK} \
    -mv \
    -P 1e-4 \
    -Oz \
    -o ${OUT_VCF}

# Index VCF
bcftools index ${OUT_VCF}

# -----------------------
# Generate Stats
# -----------------------
bcftools stats ${OUT_VCF} > ${OUT_DIR}/${SAMPLE}.stats.txt

echo "=========================================="
echo "Finished: ${SAMPLE}"
echo "VCF: ${OUT_VCF}"
echo "Stats: ${OUT_DIR}/${SAMPLE}.stats.txt"
echo "=========================================="
