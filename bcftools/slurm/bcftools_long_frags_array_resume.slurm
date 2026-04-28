#!/usr/bin/env bash
#SBATCH --job-name=tomato_bcftools_ont_frags_array
#SBATCH --array=1-63
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=36:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu

set -euo pipefail

module load bcftools/1.19
module load samtools

INPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/aligned_long_reads_150bp_bwamem2"
REF="/xdisk/rpalaniv/cjperkins1/dp_alignment/reference/S_lycopersicum_chromosomes.4.00.fa"
OUT_DIR="/xdisk/yadegari/cjperkins1/dp_vcfs/output_ont_frags_2026-04-14"

mkdir -p ${OUT_DIR}

# Robust BAM file listing
mapfile -t BAM_FILES < <(find ${INPUT_DIR} -maxdepth 1 -name "*.bam" | sort)

if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "ERROR: No BAM files found in ${INPUT_DIR}"
    exit 1
fi

echo "Total BAM files found: ${#BAM_FILES[@]}"

BAM=${BAM_FILES[$((SLURM_ARRAY_TASK_ID - 1))]}

if [ ! -f "${BAM}" ]; then
    echo "ERROR: BAM file not found: ${BAM}"
    exit 1
fi

SAMPLE=$(basename ${BAM} .bam)
OUT_VCF=${OUT_DIR}/${SAMPLE}.vcf.gz
OUT_INDEX=${OUT_VCF}.csi
OUT_STATS=${OUT_DIR}/${SAMPLE}.stats.txt

echo "=========================================="
echo "Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Processing: ${SAMPLE}"
echo "Input: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUT_VCF}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "=========================================="

# --- Resume logic ---
# A run is complete only if all three outputs exist AND the VCF is valid
IS_COMPLETE=false

if [ -f "${OUT_VCF}" ] && [ -f "${OUT_INDEX}" ] && [ -f "${OUT_STATS}" ]; then
    # Validate the VCF is not truncated/corrupt
    if bcftools view --no-version -h "${OUT_VCF}" > /dev/null 2>&1; then
        IS_COMPLETE=true
    else
        echo "WARNING: ${OUT_VCF} exists but failed validation — will rerun."
    fi
fi

if [ "${IS_COMPLETE}" = true ]; then
    echo "SKIPPING: ${SAMPLE} already complete (VCF + index + stats all present and valid)."
    echo "=========================================="
    exit 0
fi

# Clean up any partial outputs before rerunning
echo "Removing any partial outputs for ${SAMPLE}..."
rm -f "${OUT_VCF}" "${OUT_INDEX}" "${OUT_STATS}"

# --- Main variant calling ---
bcftools mpileup \
    --threads ${SLURM_CPUS_PER_TASK} \
    -f ${REF} \
    -a AD,DP,SP,ADF,ADR \
    --min-MQ 30 \
    --min-BQ 10 \
    -Ou ${BAM} | \
bcftools call \
    --threads ${SLURM_CPUS_PER_TASK} \
    -mv \
    -P 1.1e-3 \
    -Oz \
    -o ${OUT_VCF}

bcftools index ${OUT_VCF}
bcftools stats ${OUT_VCF} > ${OUT_STATS}

echo "=========================================="
echo "Finished: ${SAMPLE}"
echo "VCF: ${OUT_VCF}"
echo "Stats: ${OUT_STATS}"
echo "=========================================="
