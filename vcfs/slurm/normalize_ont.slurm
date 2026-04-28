#!/usr/bin/env bash
#SBATCH --job-name=norm_ont
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cjperkins1@arizona.edu

# Load bcftools
module load bcftools

# Paths
REF="/xdisk/rpalaniv/cjperkins1/dp_alignment/reference/S_lycopersicum_chromosomes.4.00.fa"
ONT_DIR="/xdisk/yadegari/cjperkins1/dp_vcfs/ont_raw_vcfs_2026-04-20"
OUTDIR="/xdisk/yadegari/cjperkins1/dp_vcfs/ont_raw_vcfs_2026-04-20/normalized"

mkdir -p ${OUTDIR}

# Map array task ID to chromosome
case ${SLURM_ARRAY_TASK_ID} in
    1)  CHR="ch01" ;;
    2)  CHR="ch02" ;;
    3)  CHR="ch03" ;;
    4)  CHR="ch04" ;;
    5)  CHR="ch05" ;;
    6)  CHR="ch06" ;;
    7)  CHR="ch07" ;;
    8)  CHR="ch08" ;;
    9)  CHR="ch09" ;;
    10) CHR="ch10" ;;
    11) CHR="ch11" ;;
    12) CHR="ch12" ;;
esac

echo "========================================"
echo "Job: ${SLURM_JOB_NAME}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Chromosome: ${CHR}"
echo "Node: ${SLURMD_NODENAME}"
echo "Start time: $(date)"
echo "========================================"
echo ""

# Process this chromosome
echo "Normalizing ont_SL4.0${CHR}.raw.vcf.gz..."
bcftools norm -f ${REF} -m -both ${ONT_DIR}/ont_SL4.0${CHR}.raw.vcf.gz | \
bcftools view -v snps -Oz -o ${OUTDIR}/ont_SL4.0${CHR}.snps.vcf.gz

echo "Indexing..."
bcftools index --threads ${SLURM_CPUS_PER_TASK} ${OUTDIR}/ont_SL4.0${CHR}.snps.vcf.gz

# Count SNPs
COUNT=$(bcftools view -H ${OUTDIR}/ont_SL4.0${CHR}.snps.vcf.gz | wc -l)
echo ""
echo "✓ Complete: ${CHR} has ${COUNT} SNPs"
echo "End time: $(date)"
