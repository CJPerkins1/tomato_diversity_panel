#!/usr/bin/env bash
#SBATCH --job-name=vcf2mvf_long_reads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=01:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --array=1-63
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu

# Load environment
source /home/u28/cjperkins1/.bashrc
conda activate mvftools

# VCF to MVF Conversion Script (SLURM Job Array Version)
# Uses job arrays to process multiple samples in parallel

# Configuration
BASE_DIR="/xdisk/yadegari/cjperkins1/dp_vcfs/output_vcfs_long_reads"
MVFTOOLS_PATH="/home/u28/cjperkins1/git/mvftools/mvftools.py"

# MVFtools thresholds
LOW_QUAL=10
LOW_DEPTH=5
MASK_DEPTH=4
OUT_FLAVOR="dna"

# Get array of sample directories sorted
mapfile -t vcf_dirs < <(find "$BASE_DIR" -maxdepth 1 -type d -name "*minimap2" | sort)

# Get the sample directory for this array task
SAMPLE_DIR="${vcf_dirs[$((SLURM_ARRAY_TASK_ID - 1))]}"

# If no sample directory for this task ID, exit gracefully
if [ -z "$SAMPLE_DIR" ] || [ ! -d "$SAMPLE_DIR" ]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID}: No sample to process (beyond total count)"
    exit 0
fi

sample_name=$(basename "$SAMPLE_DIR" | sed 's/_minimap2//')
vcf_file="${SAMPLE_DIR}/merge_output.vcf.gz"
mvf_file="${SAMPLE_DIR}/merge_output.mvf"

echo "============================================================"
echo "SLURM Array Task ${SLURM_ARRAY_TASK_ID}: Processing $sample_name"
echo "============================================================"
echo "Start time: $(date)"
echo "Hostname: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo ""

# Check if VCF file exists
if [ ! -f "$vcf_file" ]; then
    echo "❌ FAILED: VCF file not found"
    echo "   Expected: $vcf_file"
    exit 1
fi

# Check if conversion is already complete
if [ -f "$mvf_file" ]; then
    echo "✅ SKIPPED: Conversion already complete for $sample_name"
    mvf_size=$(du -h "$mvf_file" | cut -f1)
    echo "   Output file exists: $mvf_size"
    echo "   File: $mvf_file"
    echo "End time: $(date)"
    exit 0
fi

echo "🔄 Starting conversion..."
echo "   Input:  $vcf_file"
echo "   Output: $mvf_file"
echo ""

# Run conversion with detailed output
if python3 "$MVFTOOLS_PATH" ConvertVCF2MVF \
    --vcf "$vcf_file" \
    --out "$mvf_file" \
    --out-flavor "$OUT_FLAVOR" \
    --sample-replace "SAMPLE:${sample_name}" \
    --low-qual "$LOW_QUAL" \
    --low-depth "$LOW_DEPTH" \
    --mask-depth "$MASK_DEPTH"; then
    
    echo ""
    echo "✅ SUCCESS: Conversion completed for $sample_name"
    
    # Verify output file
    if [ -f "$mvf_file" ]; then
        mvf_size=$(du -h "$mvf_file" | cut -f1)
        echo "   Output file size: $mvf_size"
        
        # Check file is not empty
        if [ -s "$mvf_file" ]; then
            echo "   File integrity: OK"
        else
            echo "   ⚠️  WARNING: Output file is empty!"
            rm -f "$mvf_file"
            echo "End time: $(date)"
            exit 1
        fi
    fi
    
    echo "End time: $(date)"
    exit 0
else
    echo ""
    echo "❌ FAILED: Conversion failed for $sample_name"
    echo "End time: $(date)"
    exit 1
fi
