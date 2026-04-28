#!/usr/bin/env bash
#SBATCH --job-name=popgen_filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cjperkins1@arizona.edu

#=========================================================================
# Population Genomics Filtering Pipeline
# Input: Concordant SNPs (already normalized, biallelic, SNPs only)
# Output: High-quality filtered VCF for PCA, phylogenetics, structure
#=========================================================================

set -eo pipefail  # Removed 'u' flag to handle unset variables gracefully

# Load required modules
module load bcftools

# Input/Output paths
WORKDIR="/xdisk/yadegari/cjperkins1/dp_vcfs/concordant_snps"
INPUT_VCF="${WORKDIR}/merged_concordant_snps.vcf.gz"
OUTDIR="${WORKDIR}/popgen_filtered"
PREFIX="concordant_popgen"

# Create output directory
mkdir -p ${OUTDIR}

# Filtering parameters
MIN_QUAL=100          # Minimum variant quality score
MAX_MISS=0.05         # Maximum 5% missing data per SNP
MIN_MAF=0.05          # Minimum minor allele frequency (5%)
MIN_DEPTH=10          # Minimum total depth
MAX_DEPTH=10000       # Maximum total depth (removes extreme outliers)

# Number of threads (use SLURM allocation)
THREADS=${SLURM_CPUS_PER_TASK:-8}

# Safe handling of SLURM variables
MEM_ALLOCATED=${SLURM_MEM_PER_NODE:-"64GB"}

echo "========================================================================="
echo "Population Genomics Filtering Pipeline"
echo "========================================================================="
echo ""
echo "Job ID: ${SLURM_JOB_ID:-'LOCAL'}"
echo "Running on: $(hostname)"
echo "Cores: ${THREADS}"
echo "Memory: ${MEM_ALLOCATED}"
echo ""
echo "Input VCF: ${INPUT_VCF}"
echo "Output directory: ${OUTDIR}"
echo ""
echo "Filter Parameters:"
echo "  - Min Quality (QUAL): ${MIN_QUAL}"
echo "  - Max Missing per SNP: ${MAX_MISS} ($(echo "${MAX_MISS} * 100" | bc)% missing allowed)"
echo "  - Min MAF: ${MIN_MAF}"
echo "  - Depth range: ${MIN_DEPTH} - ${MAX_DEPTH}x"
echo ""
echo "========================================================================="

# Step 1: Get initial stats
echo ""
echo "[Step 1/5] Collecting input statistics..."
echo "Start time: $(date)"

TOTAL_SNPS=$(bcftools view -H ${INPUT_VCF} | wc -l)
TOTAL_SAMPLES=$(bcftools query -l ${INPUT_VCF} | wc -l)

echo "  ✓ Input SNPs: $(printf "%'d" ${TOTAL_SNPS})"
echo "  ✓ Input Samples: ${TOTAL_SAMPLES}"

# Step 2: Apply quality and depth filters
echo ""
echo "[Step 2/5] Applying quality and depth filters..."
echo "Start time: $(date)"

bcftools view \
  --threads ${THREADS} \
  -i "QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DEPTH} && INFO/DP<=${MAX_DEPTH}" \
  -Oz -o ${OUTDIR}/${PREFIX}.qual_depth.vcf.gz \
  ${INPUT_VCF}

bcftools index --threads ${THREADS} ${OUTDIR}/${PREFIX}.qual_depth.vcf.gz

AFTER_QUAL=$(bcftools view -H ${OUTDIR}/${PREFIX}.qual_depth.vcf.gz | wc -l)
echo "  ✓ SNPs after quality/depth filter: $(printf "%'d" ${AFTER_QUAL})"
echo "  ✓ Removed: $(printf "%'d" $((TOTAL_SNPS - AFTER_QUAL))) SNPs"
echo "End time: $(date)"

# Step 3: Filter by missing data
echo ""
echo "[Step 3/5] Filtering by missing data (max ${MAX_MISS})..."
echo "Start time: $(date)"

bcftools view \
  --threads ${THREADS} \
  -i "F_MISSING<=${MAX_MISS}" \
  -Oz -o ${OUTDIR}/${PREFIX}.qual_depth_miss.vcf.gz \
  ${OUTDIR}/${PREFIX}.qual_depth.vcf.gz

bcftools index --threads ${THREADS} ${OUTDIR}/${PREFIX}.qual_depth_miss.vcf.gz

AFTER_MISS=$(bcftools view -H ${OUTDIR}/${PREFIX}.qual_depth_miss.vcf.gz | wc -l)
echo "  ✓ SNPs after missing data filter: $(printf "%'d" ${AFTER_MISS})"
echo "  ✓ Removed: $(printf "%'d" $((AFTER_QUAL - AFTER_MISS))) SNPs"
echo "End time: $(date)"

# Step 4: Filter by MAF
echo ""
echo "[Step 4/5] Filtering by MAF (min ${MIN_MAF})..."
echo "Start time: $(date)"

bcftools view \
  --threads ${THREADS} \
  -q ${MIN_MAF}:minor \
  -Oz -o ${OUTDIR}/${PREFIX}.filtered.vcf.gz \
  ${OUTDIR}/${PREFIX}.qual_depth_miss.vcf.gz

bcftools index --threads ${THREADS} ${OUTDIR}/${PREFIX}.filtered.vcf.gz

AFTER_MAF=$(bcftools view -H ${OUTDIR}/${PREFIX}.filtered.vcf.gz | wc -l)
echo "  ✓ SNPs after MAF filter: $(printf "%'d" ${AFTER_MAF})"
echo "  ✓ Removed: $(printf "%'d" $((AFTER_MISS - AFTER_MAF))) SNPs"
echo "End time: $(date)"

# Step 5: Generate summary statistics
echo ""
echo "[Step 5/5] Generating summary statistics..."
echo "Start time: $(date)"

# Per-chromosome counts
echo "  Computing per-chromosome statistics..."
bcftools query -f '%CHROM\n' ${OUTDIR}/${PREFIX}.filtered.vcf.gz | \
  sort | uniq -c | \
  awk '{printf "    %s: %'"'"'d SNPs\n", $2, $1}' > ${OUTDIR}/${PREFIX}.chr_counts.txt

# Calculate actual statistics
echo "  Computing quality and depth statistics..."
bcftools query -f '%QUAL\t%INFO/DP\n' ${OUTDIR}/${PREFIX}.filtered.vcf.gz | \
  awk '{
    qual[NR]=$1; depth[NR]=$2;
    qual_sum+=$1; depth_sum+=$2;
  } 
  END {
    # Sort for median
    asort(qual, qual_sorted);
    asort(depth, depth_sorted);
    n=NR;
    mid=int(n/2);
    
    print "Quality Statistics:" > "'${OUTDIR}'/'${PREFIX}'.stats.txt";
    print "  Mean: " qual_sum/n >> "'${OUTDIR}'/'${PREFIX}'.stats.txt";
    print "  Median: " (n%2 ? qual_sorted[mid+1] : (qual_sorted[mid]+qual_sorted[mid+1])/2) >> "'${OUTDIR}'/'${PREFIX}'.stats.txt";
    print "" >> "'${OUTDIR}'/'${PREFIX}'.stats.txt";
    print "Depth Statistics:" >> "'${OUTDIR}'/'${PREFIX}'.stats.txt";
    print "  Mean: " depth_sum/n >> "'${OUTDIR}'/'${PREFIX}'.stats.txt";
    print "  Median: " (n%2 ? depth_sorted[mid+1] : (depth_sorted[mid]+depth_sorted[mid+1])/2) >> "'${OUTDIR}'/'${PREFIX}'.stats.txt";
  }'

echo "End time: $(date)"

# Summary report
cat > ${OUTDIR}/${PREFIX}.filtering_report.txt << REPORT
========================================================================
POPULATION GENOMICS FILTERING REPORT
========================================================================

Job ID: ${SLURM_JOB_ID:-'LOCAL'}
Date: $(date)
Input: ${INPUT_VCF}
Output: ${OUTDIR}/${PREFIX}.filtered.vcf.gz

------------------------------------------------------------------------
RESOURCES USED
------------------------------------------------------------------------
Cores: ${THREADS}
Memory allocated: ${MEM_ALLOCATED}

------------------------------------------------------------------------
FILTERING PARAMETERS
------------------------------------------------------------------------
Min Quality (QUAL):        ${MIN_QUAL}
Max Missing per SNP:       ${MAX_MISS} ($(echo "${MAX_MISS} * 100" | bc)% missing allowed)
Min MAF:                   ${MIN_MAF}
Depth range:               ${MIN_DEPTH} - ${MAX_DEPTH}x

------------------------------------------------------------------------
FILTERING RESULTS
------------------------------------------------------------------------
Input SNPs:                $(printf "%'12d" ${TOTAL_SNPS})
Samples:                   $(printf "%'12d" ${TOTAL_SAMPLES})

After quality/depth:       $(printf "%'12d" ${AFTER_QUAL}) (removed: $(printf "%'d" $((TOTAL_SNPS - AFTER_QUAL))))
After missing data:        $(printf "%'12d" ${AFTER_MISS}) (removed: $(printf "%'d" $((AFTER_QUAL - AFTER_MISS))))
After MAF filter:          $(printf "%'12d" ${AFTER_MAF}) (removed: $(printf "%'d" $((AFTER_MISS - AFTER_MAF))))

Final SNPs:                $(printf "%'12d" ${AFTER_MAF})
Retention rate:            $(echo "scale=2; ${AFTER_MAF} * 100 / ${TOTAL_SNPS}" | bc)%

------------------------------------------------------------------------
SNPs PER CHROMOSOME
------------------------------------------------------------------------
$(cat ${OUTDIR}/${PREFIX}.chr_counts.txt)

------------------------------------------------------------------------
VARIANT STATISTICS (filtered dataset)
------------------------------------------------------------------------
$(cat ${OUTDIR}/${PREFIX}.stats.txt)

========================================================================
REPORT

echo ""
echo "========================================================================="
echo "FILTERING COMPLETE!"
echo "========================================================================="
echo ""
echo "Summary:"
echo "  Input SNPs:     $(printf "%'12d" ${TOTAL_SNPS})"
echo "  Final SNPs:     $(printf "%'12d" ${AFTER_MAF})"
echo "  Retention:      $(echo "scale=2; ${AFTER_MAF} * 100 / ${TOTAL_SNPS}" | bc)%"
echo ""
echo "Output files:"
echo "  Main VCF:       ${OUTDIR}/${PREFIX}.filtered.vcf.gz"
echo "  Index:          ${OUTDIR}/${PREFIX}.filtered.vcf.gz.csi"
echo "  Report:         ${OUTDIR}/${PREFIX}.filtering_report.txt"
echo ""
echo "Intermediate files (can be removed if space needed):"
echo "  - ${PREFIX}.qual_depth.vcf.gz"
echo "  - ${PREFIX}.qual_depth_miss.vcf.gz"
echo ""
echo "========================================================================="
echo ""
echo "Next steps:"
echo "  1. Review: cat ${OUTDIR}/${PREFIX}.filtering_report.txt"
echo "  2. LD pruning (optional): for PCA/ADMIXTURE"
echo "  3. Convert to PLINK format for PCA"
echo ""
echo "Total runtime: $SECONDS seconds"
echo "End time: $(date)"
echo ""

