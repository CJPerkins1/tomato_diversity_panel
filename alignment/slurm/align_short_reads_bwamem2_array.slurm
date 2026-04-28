#!/usr/bin/env bash
#SBATCH --job-name=tomato_bwamem2_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --array=1-152
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=cjperkins1@arizona.edu

REFERENCE="/xdisk/rpalaniv/cjperkins1/dp_alignment/reference/S_lycopersicum_chromosomes.4.00.fa"
READS_DIR="/xdisk/rpalaniv/cjperkins1/git/pollen_quantitative_genetics/nextflow/genomes/dp/fastq"
OUTPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/aligned_short_reads_bwamem2"

source /home/u28/cjperkins1/.bashrc
conda activate bwa-mem2
module load samtools

mkdir -p "${OUTPUT_DIR}"

# Build list of R1 fastq files
mapfile -t R1_FASTQS < <(ls ${READS_DIR}/*_1.fastq.gz ${READS_DIR}/*_1.fq.gz 2>/dev/null | sort)

fq1=${R1_FASTQS[$SLURM_ARRAY_TASK_ID-1]}

if [ -z "$fq1" ]; then
    echo "No FASTQ file for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Get R2
fq2="${fq1/_1./_2.}"
if [ ! -f "$fq2" ]; then
    echo "ERROR: No R2 pair found for $fq1"
    exit 1
fi

base=$(basename "$fq1" .fastq.gz)
base=${base%.fq.gz}
base=${base%_1}

output_bam="${OUTPUT_DIR}/${base}_bwamem2.bam"
output_bai="${output_bam}.bai"

echo "Sample:    $base"
echo "R1:        $fq1"
echo "R2:        $fq2"
echo "Output:    $output_bam"
echo "Start time: $(date)"

# Resume logic
if [ -f "${output_bam}" ] && [ -f "${output_bai}" ]; then
    echo "Alignment already complete for ${base}, skipping."
    exit 0
fi

if [ -f "${output_bam}" ] && [ ! -f "${output_bai}" ]; then
    echo "BAM exists but index missing, creating index."
    samtools index "${output_bam}"
    exit 0
fi

# Alignment
bwa-mem2 mem \
    -t ${SLURM_CPUS_PER_TASK} \
    -R "@RG\tID:${base}\tSM:${base}\tPL:ILLUMINA\tLB:${base}" \
    "${REFERENCE}" \
    "$fq1" "$fq2" | \
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "${output_bam}" -

samtools index "${output_bam}"

echo "Finished ${base}"
echo "End time: $(date)"
