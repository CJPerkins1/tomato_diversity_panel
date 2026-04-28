#!/usr/bin/env bash
#SBATCH --job-name=fragment_long_reads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=48:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --array=1-63%5
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu

source /home/u28/cjperkins1/.bashrc
conda activate seqkit

# ============== PATHS ==============
READS_DIR="/xdisk/rpalaniv/cjperkins1/git/pollen_quantitative_genetics/nextflow/genomes/dp/fastq"
FRAG_OUTPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/fragmented_long_reads_150bp"
# ===================================

mkdir -p "${FRAG_OUTPUT_DIR}"

# Get array of FASTQ files (excluding paired-end _1/_2)
mapfile -t fq_files < <(
  for fq in "${READS_DIR}"/*.fastq "${READS_DIR}"/*.fq "${READS_DIR}"/*.fastq.gz "${READS_DIR}"/*.fq.gz; do
    [ -f "$fq" ] || continue
    base=$(basename "$fq" .fastq.gz)
    base=${base%.fq.gz}
    base=${base%.fastq}
    base=${base%.fq}
    [[ "$base" == *_1 ]] || [[ "$base" == *_2 ]] && continue
    echo "$fq"
  done | sort
)

# Get the file for this array task
fq="${fq_files[$((SLURM_ARRAY_TASK_ID - 1))]}"

if [ -z "$fq" ] || [ ! -f "$fq" ]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID}: No file to process"
    exit 0
fi

base=$(basename "$fq" .fastq.gz)
base=${base%.fq.gz}
base=${base%.fastq}
base=${base%.fq}

output="${FRAG_OUTPUT_DIR}/${base}_150bp.fastq.gz"

echo "============================================================"
echo "SLURM Array Task ${SLURM_ARRAY_TASK_ID}: Fragmenting $base"
echo "============================================================"
echo "Start time: $(date)"
echo "Input:  $fq"
echo "Output: $output"
echo ""

if [[ "$fq" == *.gz ]]; then
    gzip -dc "$fq" | seqkit sliding -W 150 -s 150 | gzip > "$output"
else
    cat "$fq" | seqkit sliding -W 150 -s 150 | gzip > "$output"
fi

if [ $? -eq 0 ]; then
    output_size=$(du -h "$output" | cut -f1)
    echo ""
    echo "✅ SUCCESS: Fragmentation completed for $base"
    echo "   Output size: $output_size"
    echo "End time: $(date)"
    exit 0
else
    echo ""
    echo "❌ FAILED: Fragmentation failed for $base"
    echo "End time: $(date)"
    exit 1
fi
