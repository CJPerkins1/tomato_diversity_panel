#!/usr/bin/env bash
#SBATCH --job-name=fastp_frag_qc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --array=0-62   # <-- UPDATE after counting files
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu

source /home/u28/cjperkins1/.bashrc

module load fastp

# Directories
INPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/fragmented_long_reads_150bp"
OUTPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/filtered_long_reads_150bp"

mkdir -p "${OUTPUT_DIR}"

# Get list of input files
FILES=(${INPUT_DIR}/*_150bp.fastq.gz)

# Select file for this array task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Safety check
if [[ -z "$FILE" ]]; then
  echo "No file found for task ID ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

# Extract sample name
BASENAME=$(basename "$FILE")
SAMPLE=${BASENAME%.fastq.gz}

# Output files
OUT_FASTQ="${OUTPUT_DIR}/${SAMPLE}_filtered.fastq.gz"
HTML="${OUTPUT_DIR}/${SAMPLE}.html"
JSON="${OUTPUT_DIR}/${SAMPLE}.json"

# Skip if already processed (resumable)
if [[ -f "$OUT_FASTQ" ]]; then
  echo "Skipping $SAMPLE (already processed)"
  exit 0
fi

echo "Running fastp on $SAMPLE"

fastp \
  -i "$FILE" \
  -o "$OUT_FASTQ" \
  --disable_adapter_trimming \
  --length_required 140 \
  -q 7 \
  -u 40 \
  -n 10 \
  -h "$HTML" \
  -j "$JSON" \
  -w ${SLURM_CPUS_PER_TASK}

echo "Finished $SAMPLE"
