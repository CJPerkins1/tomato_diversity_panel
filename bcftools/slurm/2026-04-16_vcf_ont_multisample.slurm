#!/usr/bin/env bash
#SBATCH --job-name=vcf_ont_multisample
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=96:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cjperkins1@arizona.edu

set -euo pipefail

module load bcftools/1.19
module load samtools/1.19

# ============================================
# PATHS
# ============================================
INPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/aligned_long_reads_150bp_bwamem2"
REF="/xdisk/rpalaniv/cjperkins1/dp_alignment/reference/S_lycopersicum_chromosomes.4.00.fa"
OUT_DIR="/xdisk/yadegari/cjperkins1/dp_vcfs/ont_raw_vcfs_$(date +%Y-%m-%d)"

mkdir -p ${OUT_DIR}

# ============================================
# CHROMOSOME ARRAY
# ============================================
CHR=$(printf "SL4.0ch%02d" ${SLURM_ARRAY_TASK_ID})

echo "=========================================="
echo "Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Chromosome: ${CHR}"
echo "Reference: ${REF}"
echo "Output Dir: ${OUT_DIR}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "=========================================="

# ============================================
# CREATE BAM LIST (all samples in directory)
# ============================================
BAM_LIST="${OUT_DIR}/ont_bam_list.txt"

if [ ! -f "${BAM_LIST}" ]; then
    echo "Creating BAM list..."
    find ${INPUT_DIR} -maxdepth 1 -name "*.bam" | sort > ${BAM_LIST}
    N_BAMS=$(wc -l < ${BAM_LIST})
    echo "Found ${N_BAMS} BAM files"
else
    N_BAMS=$(wc -l < ${BAM_LIST})
    echo "Using existing BAM list (${N_BAMS} files)"
fi

if [ ${N_BAMS} -eq 0 ]; then
    echo "ERROR: No BAM files found!"
    exit 1
fi

# ============================================
# OUTPUT FILES
# ============================================
OUT_VCF="${OUT_DIR}/ont_${CHR}.raw.vcf.gz"
OUT_INDEX="${OUT_VCF}.csi"

# ============================================
# RESUME LOGIC
# ============================================
if [ -f "${OUT_VCF}" ] && [ -f "${OUT_INDEX}" ]; then
    if bcftools view --no-version -h "${OUT_VCF}" > /dev/null 2>&1; then
        echo "SKIPPING: ${CHR} already complete and validated"
        exit 0
    else
        echo "WARNING: Existing VCF failed validation, rerunning..."
        rm -f "${OUT_VCF}" "${OUT_INDEX}"
    fi
fi

# ============================================
# VARIANT CALLING
# ============================================
echo "Starting variant calling for ${CHR}..."
echo "$(date): Running bcftools mpileup..."

bcftools mpileup \
    --threads ${SLURM_CPUS_PER_TASK} \
    -B \
    -I \
    --min-MQ 20 \
    --min-BQ 20 \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --fasta-ref ${REF} \
    --regions ${CHR} \
    --bam-list ${BAM_LIST} \
    --output-type u | \
bcftools call \
    --threads ${SLURM_CPUS_PER_TASK} \
    --multiallelic-caller \
    --variants-only \
    --format-fields GQ,GP \
    --output-type z \
    --output ${OUT_VCF}

echo "$(date): Indexing VCF..."
bcftools index --threads ${SLURM_CPUS_PER_TASK} ${OUT_VCF}

# ============================================
# VALIDATION & STATS
# ============================================
echo "$(date): Validating VCF..."
if bcftools view --no-version -h "${OUT_VCF}" > /dev/null 2>&1; then
    N_VARIANTS=$(bcftools view -H ${OUT_VCF} | wc -l)
    echo "SUCCESS: ${CHR} complete with ${N_VARIANTS} variants"
else
    echo "ERROR: VCF validation failed!"
    exit 1
fi

echo "=========================================="
echo "Finished: ${CHR}"
echo "Output: ${OUT_VCF}"
echo "=========================================="
