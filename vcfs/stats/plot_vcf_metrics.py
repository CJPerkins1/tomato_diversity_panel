import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='Plot VCF quality metrics from bcftools query output')
parser.add_argument('-i', '--input', required=True, help='Input file (output from bcftools query)')
parser.add_argument('-o', '--output', default=None, help='Output prefix for plot (default: based on input filename)')
args = parser.parse_args()

# Set output prefix
if args.output:
    output_prefix = args.output
else:
    output_prefix = os.path.splitext(os.path.basename(args.input))[0]
    if output_prefix.endswith('.txt'):
        output_prefix = os.path.splitext(output_prefix)[0]

# Read the data with proper handling of missing values
print(f"Loading variant metrics from: {args.input}")
df = pd.read_csv(args.input, 
                 sep='\t', 
                 names=['CHROM', 'POS', 'QUAL', 'DP', 'MQ', 'MQ0F', 'VDB', 'RPBZ', 'MQBZ', 'BQBZ'],
                 na_values='.')  # Treat '.' as missing values

print(f"Total variants: {len(df):,}")
print("\nBasic statistics:")
print(df.describe())

# Set up the plotting style
sns.set_style("whitegrid")
fig, axes = plt.subplots(3, 3, figsize=(15, 12))
chrom_name = df['CHROM'].iloc[0] if len(df) > 0 else 'Unknown'
fig.suptitle(f'Variant Quality Metrics - {chrom_name}', fontsize=16, y=1.00)

# 1. QUAL distribution
ax = axes[0, 0]
df['QUAL'].hist(bins=100, ax=ax, edgecolor='black', alpha=0.7)
ax.axvline(df['QUAL'].quantile(0.05), color='red', linestyle='--', label='5th percentile')
ax.axvline(df['QUAL'].quantile(0.10), color='orange', linestyle='--', label='10th percentile')
ax.set_xlabel('QUAL Score')
ax.set_ylabel('Count')
ax.set_title('Variant Quality (QUAL)')
ax.legend()
ax.set_yscale('log')

# 2. DP distribution
ax = axes[0, 1]
df['DP'].hist(bins=100, ax=ax, edgecolor='black', alpha=0.7)
median_dp = df['DP'].median()
ax.axvline(median_dp * 0.5, color='red', linestyle='--', label='0.5x median')
ax.axvline(median_dp * 2, color='red', linestyle='--', label='2x median')
ax.set_xlabel('Total Depth (DP)')
ax.set_ylabel('Count')
ax.set_title(f'Total Depth (median: {median_dp:.0f})')
ax.legend()
ax.set_yscale('log')

# 3. MQ distribution
ax = axes[0, 2]
df['MQ'].hist(bins=60, ax=ax, edgecolor='black', alpha=0.7)
ax.axvline(40, color='red', linestyle='--', label='MQ=40 (common cutoff)')
ax.set_xlabel('Mapping Quality (MQ)')
ax.set_ylabel('Count')
ax.set_title('Mapping Quality')
ax.legend()
ax.set_yscale('log')

# 4. MQ0F distribution
ax = axes[1, 0]
df['MQ0F'].hist(bins=100, ax=ax, edgecolor='black', alpha=0.7)
ax.axvline(0.1, color='red', linestyle='--', label='MQ0F=0.1')
ax.set_xlabel('MQ0 Fraction')
ax.set_ylabel('Count')
ax.set_title('Fraction of MQ0 Reads (lower is better)')
ax.legend()
ax.set_yscale('log')

# 5. VDB distribution
ax = axes[1, 1]
vdb_clean = df['VDB'].dropna()
if len(vdb_clean) > 0:
    vdb_clean.hist(bins=100, ax=ax, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Variant Distance Bias (VDB)')
    ax.set_ylabel('Count')
    ax.set_title(f'VDB (higher is better, n={len(vdb_clean):,})')
    ax.set_yscale('log')
else:
    ax.text(0.5, 0.5, 'No VDB data', ha='center', va='center')

# 6. RPBZ distribution
ax = axes[1, 2]
rpbz_clean = df['RPBZ'].dropna()
if len(rpbz_clean) > 0:
    rpbz_clean.hist(bins=100, ax=ax, edgecolor='black', alpha=0.7)
    ax.axvline(-2, color='red', linestyle='--', alpha=0.5)
    ax.axvline(2, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Read Position Bias Z-score')
    ax.set_ylabel('Count')
    ax.set_title(f'RPBZ (closer to 0 is better, n={len(rpbz_clean):,})')
    ax.set_yscale('log')
else:
    ax.text(0.5, 0.5, 'No RPBZ data', ha='center', va='center')

# 7. MQBZ distribution
ax = axes[2, 0]
mqbz_clean = df['MQBZ'].dropna()
if len(mqbz_clean) > 0:
    mqbz_clean.hist(bins=100, ax=ax, edgecolor='black', alpha=0.7)
    ax.axvline(-2, color='red', linestyle='--', alpha=0.5)
    ax.axvline(2, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Mapping Quality Bias Z-score')
    ax.set_ylabel('Count')
    ax.set_title(f'MQBZ (closer to 0 is better, n={len(mqbz_clean):,})')
    ax.set_yscale('log')
else:
    ax.text(0.5, 0.5, 'No MQBZ data', ha='center', va='center')

# 8. BQBZ distribution
ax = axes[2, 1]
bqbz_clean = df['BQBZ'].dropna()
if len(bqbz_clean) > 0:
    bqbz_clean.hist(bins=100, ax=ax, edgecolor='black', alpha=0.7)
    ax.axvline(-2, color='red', linestyle='--', alpha=0.5)
    ax.axvline(2, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Base Quality Bias Z-score')
    ax.set_ylabel('Count')
    ax.set_title(f'BQBZ (closer to 0 is better, n={len(bqbz_clean):,})')
    ax.set_yscale('log')
else:
    ax.text(0.5, 0.5, 'No BQBZ data', ha='center', va='center')

# 9. QUAL vs DP scatter (sampled)
ax = axes[2, 2]
sample_size = min(10000, len(df))
sample = df.sample(n=sample_size)
ax.scatter(sample['DP'], sample['QUAL'], alpha=0.1, s=1)
ax.set_xlabel('Total Depth (DP)')
ax.set_ylabel('Count')
ax.set_title(f'QUAL vs DP (n={sample_size:,})')
ax.set_xscale('log')
ax.set_yscale('log')

plt.tight_layout()
output_file = f'{output_prefix}_quality_distributions.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPlot saved as: {output_file}")

# Print suggested filtering thresholds
print("\n" + "="*60)
print("SUGGESTED FILTERING THRESHOLDS (based on data distribution)")
print("="*60)
print(f"\nQUAL > {df['QUAL'].quantile(0.10):.1f}  (removes bottom 10%)")
print(f"DP between {median_dp * 0.5:.0f} - {median_dp * 2:.0f}  (0.5x to 2x median depth)")
print(f"MQ > 40  (standard cutoff for good mapping)")
print(f"MQ0F < 0.1  (< 10% reads with MQ=0)")

if len(rpbz_clean) > 0:
    print(f"abs(RPBZ) < 3  (remove extreme read position bias)")
if len(mqbz_clean) > 0:
    print(f"abs(MQBZ) < 3  (remove extreme mapping quality bias)")
if len(bqbz_clean) > 0:
    print(f"abs(BQBZ) < 3  (remove extreme base quality bias)")
    
print("\nNote: These are starting points - adjust based on your plots!")
