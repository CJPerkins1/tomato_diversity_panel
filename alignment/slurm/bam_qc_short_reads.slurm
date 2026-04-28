#!/usr/bin/env bash
#SBATCH --job-name=tomato_bam_qc_illumina
#SBATCH --array=1-152
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=cjperkins1@arizona.edu

set -euo pipefail

module load samtools/1.19

# ----------------------------
# INPUT / OUTPUT
# ----------------------------
INPUT_DIR="/xdisk/yadegari/cjperkins1/dp_alignment/aligned_short_reads_bwamem2"
PARENT_DIR=$(dirname "$INPUT_DIR")
OUTDIR="${PARENT_DIR}/Bam_QC/short_reads"
mkdir -p "${OUTDIR}"

DATA_TYPE="Illumina_shortreads"
QC_PIPELINE="samtools_bam_qc_v3"
REF_NAME="S_lycopersicum_4.00"

# Shared output file — header written once via flock
COMBINED="${OUTDIR}/Illumina_all_samples_bam_qc.tsv"
HEADER="sample\tdata_type\tqc_pipeline\treference\ttotal_reads\tmapped_reads\tmapped_pct\tproper_pairs\tproper_pairs_pct\tmean_mapq\tmedian_mapq\tmq0\tmq0_pct\tmq1_19\tmq1_19_pct\tmq20_29\tmq20_29_pct\tmq30_59\tmq30_59_pct\tmq60\tmq60_pct\tcoverage_mean\tcoverage_median\tpct_zero_coverage"

flock "${COMBINED}.lock" bash -c '
    if [ ! -f "$1" ]; then
        echo -e "$2" > "$1"
    fi
' _ "${COMBINED}" "${HEADER}"

# ----------------------------
# GET BAM FILE
# ----------------------------
BAMS=(${INPUT_DIR}/*.bam)
BAM=${BAMS[$SLURM_ARRAY_TASK_ID-1]}
SAMPLE=$(basename "${BAM}" .bam)

echo "[$(date)] Starting QC for ${SAMPLE}"

# ----------------------------
# 1. FLAGSTAT
# ----------------------------
FLAGSTAT=$(samtools flagstat "${BAM}")

TOTAL_READS=$(echo "$FLAGSTAT"   | awk '/in total/             {print $1}')
PROPER_PAIRS=$(echo "$FLAGSTAT"  | awk '/properly paired/      {print $1}')
MAPPED_READS=$(echo "$FLAGSTAT"  | awk '/\+ [0-9]+ mapped/     {print $1; exit}')
MAPPED_PCT=$(echo "$FLAGSTAT"    | awk '/\+ [0-9]+ mapped/     {gsub(/[()%]/,"",$5); print $5; exit}')

if [[ -z "${MAPPED_READS}" || "${MAPPED_READS}" -eq 0 ]]; then
    echo "[$(date)] WARNING: MAPPED_READS is 0 or empty for ${SAMPLE} — check BAM integrity"
    MAPPED_READS=0
    MAPPED_PCT="0"
    PROPER_PAIRS_PCT="0"
else
    PROPER_PAIRS_PCT=$(awk "BEGIN {printf \"%.4f\", ${PROPER_PAIRS} / ${TOTAL_READS} * 100}")
fi

echo "[$(date)] Flagstat done — ${MAPPED_READS} mapped reads (${MAPPED_PCT}%), ${PROPER_PAIRS} properly paired"

# ----------------------------
# 2. MAPQ
# ----------------------------
MAPQ_STATS=$(samtools view -F 4 "${BAM}" | awk '
{
    mq = $5
    freq[mq]++
    sum += mq
    total++

    if      (mq == 0)  mq0++
    else if (mq < 20)  mq1_19++
    else if (mq < 30)  mq20_29++
    else if (mq < 60)  mq30_59++
    else if (mq == 60) mq60++
}
END {
    if (total == 0) {
        print "0\t0\t0\t0\t0\t0\t0"
        exit
    }

    mean = sum / total

    target = (total + 1) / 2
    cumsum = 0
    median = 0
    for (mq = 0; mq <= 60; mq++) {
        if (mq in freq) {
            cumsum += freq[mq]
            if (cumsum >= target && median == 0)
                median = mq
        }
    }

    printf "%.4f\t%d\t%d\t%d\t%d\t%d\t%d\n",
        mean, median,
        (mq0     ? mq0     : 0),
        (mq1_19  ? mq1_19  : 0),
        (mq20_29 ? mq20_29 : 0),
        (mq30_59 ? mq30_59 : 0),
        (mq60    ? mq60    : 0)
}')

MEAN_MAPQ=$(echo   "$MAPQ_STATS" | cut -f1)
MEDIAN_MAPQ=$(echo "$MAPQ_STATS" | cut -f2)
MQ0=$(echo         "$MAPQ_STATS" | cut -f3)
MQ1_19=$(echo      "$MAPQ_STATS" | cut -f4)
MQ20_29=$(echo     "$MAPQ_STATS" | cut -f5)
MQ30_59=$(echo     "$MAPQ_STATS" | cut -f6)
MQ60=$(echo        "$MAPQ_STATS" | cut -f7)

if [[ "${MAPPED_READS}" -eq 0 ]]; then
    MQ0_PCT="0"; MQ1_19_PCT="0"; MQ20_29_PCT="0"; MQ30_59_PCT="0"; MQ60_PCT="0"
else
    MQ0_PCT=$(awk     "BEGIN {printf \"%.4f\", ${MQ0}     / ${MAPPED_READS} * 100}")
    MQ1_19_PCT=$(awk  "BEGIN {printf \"%.4f\", ${MQ1_19}  / ${MAPPED_READS} * 100}")
    MQ20_29_PCT=$(awk "BEGIN {printf \"%.4f\", ${MQ20_29} / ${MAPPED_READS} * 100}")
    MQ30_59_PCT=$(awk "BEGIN {printf \"%.4f\", ${MQ30_59} / ${MAPPED_READS} * 100}")
    MQ60_PCT=$(awk    "BEGIN {printf \"%.4f\", ${MQ60}    / ${MAPPED_READS} * 100}")
fi

echo "[$(date)] MAPQ stats done — mean ${MEAN_MAPQ}, median ${MEDIAN_MAPQ}"

# ----------------------------
# 3. COVERAGE
# ----------------------------
COVERAGE_STATS=$(samtools depth -a "${BAM}" | awk '
{
    depth = $3
    sum  += depth
    count++
    freq[depth]++
    if (depth == 0) zero++
}
END {
    if (count == 0) {
        print "0\t0\t0"
        exit
    }

    mean = sum / count

    target  = (count + 1) / 2
    cumsum  = 0
    median  = 0
    for (d = 0; d <= 100000; d++) {
        if (d in freq) {
            cumsum += freq[d]
            if (cumsum >= target && median == 0)
                median = d
        }
    }

    pct_zero = zero / count * 100
    printf "%.4f\t%d\t%.4f\n", mean, median, pct_zero
}')

COVERAGE_MEAN=$(echo   "$COVERAGE_STATS" | cut -f1)
COVERAGE_MEDIAN=$(echo "$COVERAGE_STATS" | cut -f2)
PCT_ZERO_COV=$(echo    "$COVERAGE_STATS" | cut -f3)

echo "[$(date)] Coverage stats done — mean ${COVERAGE_MEAN}x, ${PCT_ZERO_COV}% zero coverage"

# ----------------------------
# 4. WRITE OUTPUT ROW
# ----------------------------
ROW="${SAMPLE}\t${DATA_TYPE}\t${QC_PIPELINE}\t${REF_NAME}\t${TOTAL_READS}\t${MAPPED_READS}\t${MAPPED_PCT}\t${PROPER_PAIRS}\t${PROPER_PAIRS_PCT}\t${MEAN_MAPQ}\t${MEDIAN_MAPQ}\t${MQ0}\t${MQ0_PCT}\t${MQ1_19}\t${MQ1_19_PCT}\t${MQ20_29}\t${MQ20_29_PCT}\t${MQ30_59}\t${MQ30_59_PCT}\t${MQ60}\t${MQ60_PCT}\t${COVERAGE_MEAN}\t${COVERAGE_MEDIAN}\t${PCT_ZERO_COV}"

flock "${COMBINED}.lock" echo -e "${ROW}" >> "${COMBINED}"

echo "[$(date)] QC complete for ${SAMPLE} — results written to ${COMBINED}"
