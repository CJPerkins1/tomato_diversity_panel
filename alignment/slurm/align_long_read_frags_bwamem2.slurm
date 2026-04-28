#!/usr/bin/env bash
#SBATCH --job-name=tomato_bwamem2_150bp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --array=1-63
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=cjperkins1@arizona.edu

REFERENCE="/xdisk/rpalaniv/cjperkins1/dp_alignment/reference/S_lycopersicum_chromosomes.4.00.fa"

READS_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/filtered_long_reads_150bp"
OUTPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/aligned_long_reads_150bp_bwamem2"

source /home/u28/cjperkins1/.bashrc
conda activate bwa-mem2
module load samtools

mkdir -p "${OUTPUT_DIR}"

# Build list of FASTQ files
mapfile -t FASTQS < <(ls ${READS_DIR}/*.fastq ${READS_DIR}/*.fq ${READS_DIR}/*.fastq.gz ${READS_DIR}/*.fq.gz 2>/dev/null)

fq=${FASTQS[$SLURM_ARRAY_TASK_ID-1]}

if [ -z "$fq" ]; then
    echo "No FASTQ file for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

base=$(basename "$fq" .fastq.gz)
base=${base%.fq.gz}
base=${base%.fastq}
base=${base%.fq}

output_bam="${OUTPUT_DIR}/${base}_bwamem2.bam"
output_bai="${output_bam}.bai"

echo "Processing file: $fq"
echo "Output BAM: $output_bam"
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

# Check reference index exists
if [ ! -f "${REFERENCE}.bwt.2bit.64" ]; then
    echo "ERROR: BWA-MEM2 index not found for reference."
    echo "Run: bwa-mem2 index ${REFERENCE}"
    exit 1
fi

# Alignment
bwa-mem2 mem \
    -t ${SLURM_CPUS_PER_TASK} \
    -R "@RG\tID:${base}\tSM:${base}\tPL:ONT\tLB:${base}" \
    "${REFERENCE}" \
    "$fq" | \
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "${output_bam}" -

samtools index "${output_bam}"

echo "Finished ${base}"
echo "End time: $(date)"
