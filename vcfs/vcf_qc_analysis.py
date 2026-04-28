import argparse
import subprocess
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from collections import defaultdict
import re
import os
import sys

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Quality Control Analysis for VCF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic usage with input VCF
  python vcf_qc_analysis.py -i merged_snps.vcf.gz
  
  # Specify custom output directory
  python vcf_qc_analysis.py -i merged_snps.vcf.gz -o /path/to/output/
  
  # Specify custom output prefix
  python vcf_qc_analysis.py -i merged_snps.vcf.gz -o qc_results
  
  # Adjust sampling rate (default is 100, meaning every 100th SNP)
  python vcf_qc_analysis.py -i merged_snps.vcf.gz --sample-rate 50
        '''
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input VCF file (can be .vcf or .vcf.gz)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default=None,
        help='Output directory or prefix for results (default: same directory as input with input basename)'
    )
    
    parser.add_argument(
        '--sample-rate',
        type=int,
        default=100,
        help='Sample every Nth SNP for plotting (default: 100). Lower values = more detail but slower.'
    )
    
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='DPI for output plot (default: 300)'
    )
    
    return parser.parse_args()

def determine_output_paths(input_vcf, output_arg):
    """Determine output file paths based on input and output arguments"""
    
    if output_arg is None:
        # Default: use input directory and basename
        input_dir = os.path.dirname(os.path.abspath(input_vcf))
        input_base = os.path.basename(input_vcf)
        # Remove .vcf.gz or .vcf extension
        prefix = input_base.replace('.vcf.gz', '').replace('.vcf', '')
        output_prefix = os.path.join(input_dir, f"{prefix}_qc")
    
    elif os.path.isdir(output_arg):
        # Output is a directory
        input_base = os.path.basename(input_vcf)
        prefix = input_base.replace('.vcf.gz', '').replace('.vcf', '')
        output_prefix = os.path.join(output_arg, f"{prefix}_qc")
    
    else:
        # Output is a prefix (could include directory)
        output_dir = os.path.dirname(output_arg)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        output_prefix = output_arg
    
    plot_file = f"{output_prefix}_report.png"
    stats_file = f"{output_prefix}_stats.txt"
    
    return plot_file, stats_file

def main():
    args = parse_arguments()
    
    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Determine output paths
    output_plot, output_stats = determine_output_paths(args.input, args.output)
    
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (14, 10)
    
    print("=" * 70)
    print("Quality Control Analysis: VCF SNPs")
    print("=" * 70)
    print(f"\nInput VCF:      {args.input}")
    print(f"Output plot:    {output_plot}")
    print(f"Output stats:   {output_stats}")
    print(f"Sampling rate:  Every {args.sample_rate} SNPs\n")
    
    # Initialize data structures
    chromosomes = []
    positions = []
    qualities = []
    filters = []
    depths = []
    allele_freqs = []
    missing_rates = []
    chr_counts = defaultdict(int)
    
    print("[1/3] Parsing VCF and collecting statistics...")
    if args.sample_rate > 1:
        print(f"      (sampling every {args.sample_rate}th SNP for memory efficiency)")
    
    line_count = 0
    snp_count = 0
    sample_count = 0
    
    # Determine if file is gzipped
    if args.input.endswith('.gz'):
        open_func = lambda f: gzip.open(f, 'rt')
    else:
        open_func = lambda f: open(f, 'r')
    
    with open_func(args.input) as f:
        for line in f:
            if line.startswith('##'):
                continue
            
            if line.startswith('#CHROM'):
                # Get sample names
                samples = line.strip().split('\t')[9:]
                sample_count = len(samples)
                print(f"      ✓ Found {sample_count} samples")
                continue
            
            line_count += 1
            
            # Sample every Nth SNP to keep memory manageable
            if line_count % args.sample_rate != 0:
                snp_count += 1
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            qual = fields[5]
            filt = fields[6]
            info = fields[7]
            genotypes = fields[9:]
            
            # Store data
            chromosomes.append(chrom)
            positions.append(pos)
            chr_counts[chrom] += 1
            filters.append(filt)
            
            # Parse quality
            if qual != '.':
                qualities.append(float(qual))
            
            # Parse depth from INFO
            dp_match = re.search(r'DP=(\d+)', info)
            if dp_match:
                depths.append(int(dp_match.group(1)))
            
            # Calculate allele frequency and missing rate
            gt_calls = [g.split(':')[0] for g in genotypes]
            total_alleles = 0
            alt_alleles = 0
            missing = 0
            
            for gt in gt_calls:
                if gt in ['./.', '.|.', '.']:
                    missing += 1
                else:
                    alleles = gt.replace('|', '/').split('/')
                    for a in alleles:
                        if a != '.':
                            total_alleles += 1
                            if a != '0':
                                alt_alleles += 1
            
            if total_alleles > 0:
                allele_freqs.append(alt_alleles / total_alleles)
            missing_rates.append(missing / len(genotypes))
            
            snp_count += 1
            
            if snp_count % 100000 == 0:
                print(f"      Processed {snp_count:,} SNPs...")
    
    print(f"\n      ✓ Total SNPs processed: {snp_count:,}")
    print(f"      ✓ SNPs sampled for plots: {len(chromosomes):,}")
    
    print("\n[2/3] Computing summary statistics...")
    
    # Summary statistics
    stats = {
        'Total SNPs': snp_count,
        'Samples': sample_count,
        'Mean Quality': np.mean(qualities) if qualities else 'N/A',
        'Median Quality': np.median(qualities) if qualities else 'N/A',
        'Mean Depth': np.mean(depths) if depths else 'N/A',
        'Median Depth': np.median(depths) if depths else 'N/A',
        'Mean Missing Rate': np.mean(missing_rates) if missing_rates else 'N/A',
        'PASS filter rate': (filters.count('PASS') / len(filters) * 100) if filters else 'N/A'
    }
    
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"{key:.<40} {value:.2f}")
        else:
            print(f"{key:.<40} {value}")
    
    print("\n" + "=" * 70)
    print("SNPs PER CHROMOSOME")
    print("=" * 70)
    chr_order = sorted(chr_counts.keys(), key=lambda x: (len(x), x))
    for chrom in chr_order:
        count = chr_counts[chrom] * args.sample_rate
        pct = (count / snp_count) * 100
        print(f"{chrom:.<20} {count:>12,} ({pct:>5.2f}%)")
    
    print("\n[3/3] Generating QC plots...")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    
    # 1. SNP distribution across chromosomes
    ax1 = plt.subplot(3, 3, 1)
    chr_counts_plot = {k: v * args.sample_rate for k, v in chr_counts.items()}
    chr_sorted = sorted(chr_counts_plot.items(), key=lambda x: (len(x[0]), x[0]))
    chrs, counts = zip(*chr_sorted)
    bars = ax1.bar(range(len(chrs)), counts, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.set_xticks(range(len(chrs)))
    ax1.set_xticklabels(chrs, rotation=45, ha='right')
    ax1.set_ylabel('SNP Count')
    ax1.set_title('SNPs per Chromosome', fontweight='bold', fontsize=12)
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x/1000)}K' if x >= 1000 else f'{int(x)}'))
    ax1.grid(axis='y', alpha=0.3)
    
    # 2. Quality score distribution
    ax2 = plt.subplot(3, 3, 2)
    if qualities:
        ax2.hist(qualities, bins=50, color='coral', edgecolor='black', alpha=0.7)
        ax2.axvline(np.median(qualities), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(qualities):.1f}')
        ax2.set_xlabel('Quality Score')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Quality Score Distribution', fontweight='bold', fontsize=12)
        ax2.legend()
        ax2.grid(axis='y', alpha=0.3)
    
    # 3. Depth distribution
    ax3 = plt.subplot(3, 3, 3)
    if depths:
        depth_clip = [min(d, 1000) for d in depths]  # Clip at 1000 for visualization
        ax3.hist(depth_clip, bins=50, color='lightgreen', edgecolor='black', alpha=0.7)
        ax3.axvline(np.median(depths), color='darkgreen', linestyle='--', linewidth=2, label=f'Median: {np.median(depths):.1f}')
        ax3.set_xlabel('Read Depth (DP)')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Sequencing Depth Distribution', fontweight='bold', fontsize=12)
        ax3.legend()
        ax3.grid(axis='y', alpha=0.3)
    
    # 4. Allele frequency spectrum
    ax4 = plt.subplot(3, 3, 4)
    if allele_freqs:
        ax4.hist(allele_freqs, bins=50, color='purple', edgecolor='black', alpha=0.7)
        ax4.set_xlabel('Alternate Allele Frequency')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Allele Frequency Spectrum', fontweight='bold', fontsize=12)
        ax4.grid(axis='y', alpha=0.3)
    
    # 5. Missing data rate
    ax5 = plt.subplot(3, 3, 5)
    if missing_rates:
        ax5.hist(missing_rates, bins=50, color='orange', edgecolor='black', alpha=0.7)
        ax5.axvline(np.median(missing_rates), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(missing_rates):.3f}')
        ax5.set_xlabel('Missing Genotype Rate')
        ax5.set_ylabel('Frequency')
        ax5.set_title('Missing Data Distribution', fontweight='bold', fontsize=12)
        ax5.legend()
        ax5.grid(axis='y', alpha=0.3)
    
    # 6. Filter status
    ax6 = plt.subplot(3, 3, 6)
    if filters:
        filter_counts = pd.Series(filters).value_counts()
        colors_filter = ['green' if x == 'PASS' else 'red' for x in filter_counts.index]
        ax6.bar(range(len(filter_counts)), filter_counts.values, color=colors_filter, edgecolor='black', alpha=0.7)
        ax6.set_xticks(range(len(filter_counts)))
        ax6.set_xticklabels(filter_counts.index, rotation=45, ha='right')
        ax6.set_ylabel('Count')
        ax6.set_title('Filter Status', fontweight='bold', fontsize=12)
        ax6.grid(axis='y', alpha=0.3)
    
    # 7. Quality vs Depth scatter
    ax7 = plt.subplot(3, 3, 7)
    if qualities and depths:
        sample_indices = np.random.choice(len(qualities), min(10000, len(qualities)), replace=False)
        qual_sample = [qualities[i] for i in sample_indices]
        depth_sample = [min(depths[i], 1000) for i in sample_indices]
        ax7.scatter(depth_sample, qual_sample, alpha=0.3, s=1, color='navy')
        ax7.set_xlabel('Read Depth')
        ax7.set_ylabel('Quality Score')
        ax7.set_title('Quality vs Depth', fontweight='bold', fontsize=12)
        ax7.grid(alpha=0.3)
    
    # 8. Cumulative allele frequency
    ax8 = plt.subplot(3, 3, 8)
    if allele_freqs:
        sorted_af = sorted(allele_freqs)
        cumulative = np.arange(1, len(sorted_af) + 1) / len(sorted_af)
        ax8.plot(sorted_af, cumulative, color='teal', linewidth=2)
        ax8.set_xlabel('Alternate Allele Frequency')
        ax8.set_ylabel('Cumulative Proportion')
        ax8.set_title('Cumulative Allele Frequency', fontweight='bold', fontsize=12)
        ax8.grid(alpha=0.3)
    
    # 9. SNP density along chromosomes (for first few chromosomes)
    ax9 = plt.subplot(3, 3, 9)
    chr_data = defaultdict(list)
    for c, p in zip(chromosomes, positions):
        chr_data[c].append(p)
    
    # Plot density for first 5 chromosomes
    plot_chrs = sorted(list(chr_data.keys()), key=lambda x: (len(x), x))[:5]
    for i, chrom in enumerate(plot_chrs):
        pos_chr = sorted(chr_data[chrom])
        if len(pos_chr) > 1:
            bins = np.linspace(min(pos_chr), max(pos_chr), 100)
            counts, _ = np.histogram(pos_chr, bins=bins)
            bin_centers = (bins[:-1] + bins[1:]) / 2
            ax9.plot(bin_centers / 1e6, counts, label=chrom, alpha=0.7)
    
    ax9.set_xlabel('Position (Mb)')
    ax9.set_ylabel('SNP Count per Bin')
    ax9.set_title('SNP Density Along Chromosomes', fontweight='bold', fontsize=12)
    ax9.legend(fontsize=8)
    ax9.grid(alpha=0.3)
    
    plt.suptitle(f'VCF Quality Control Report\n{snp_count:,} SNPs | {sample_count} Samples', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_plot, dpi=args.dpi, bbox_inches='tight')
    print(f"\n✓ QC plot saved: {output_plot}")
    
    # Save summary statistics
    with open(output_stats, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("VCF QUALITY CONTROL REPORT\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"VCF File: {args.input}\n")
        f.write(f"Sampling Rate: Every {args.sample_rate} SNPs\n\n")
        
        f.write("=" * 70 + "\n")
        f.write("SUMMARY STATISTICS\n")
        f.write("=" * 70 + "\n")
        for key, value in stats.items():
            if isinstance(value, float):
                f.write(f"{key:.<45} {value:.2f}\n")
            else:
                f.write(f"{key:.<45} {value}\n")
        
        f.write("\n" + "=" * 70 + "\n")
        f.write("SNPs PER CHROMOSOME\n")
        f.write("=" * 70 + "\n")
        for chrom in chr_order:
            count = chr_counts[chrom] * args.sample_rate
            pct = (count / snp_count) * 100
            f.write(f"{chrom:.<25} {count:>12,} ({pct:>5.2f}%)\n")
        
        if qualities:
            f.write("\n" + "=" * 70 + "\n")
            f.write("QUALITY SCORE PERCENTILES\n")
            f.write("=" * 70 + "\n")
            percentiles = [10, 25, 50, 75, 90, 95, 99]
            for p in percentiles:
                val = np.percentile(qualities, p)
                f.write(f"{p}th percentile{'':.>30} {val:.2f}\n")
        
        if depths:
            f.write("\n" + "=" * 70 + "\n")
            f.write("DEPTH PERCENTILES\n")
            f.write("=" * 70 + "\n")
            for p in percentiles:
                val = np.percentile(depths, p)
                f.write(f"{p}th percentile{'':.>30} {val:.2f}\n")
    
    print(f"✓ Statistics saved: {output_stats}")
    
    print("\n" + "=" * 70)
    print("QC ANALYSIS COMPLETE!")
    print("=" * 70)
    print("\nKey Findings:")
    print(f"  • Total SNPs analyzed: {snp_count:,}")
    print(f"  • Total samples: {sample_count}")
    if qualities:
        print(f"  • Median quality score: {np.median(qualities):.1f}")
    if depths:
        print(f"  • Median depth: {np.median(depths):.1f}x")
    if missing_rates:
        print(f"  • Median missing rate: {np.median(missing_rates):.3f}")
    if filters:
        pass_rate = (filters.count('PASS') / len(filters) * 100)
        print(f"  • PASS filter rate: {pass_rate:.1f}%")
    
    print("\n" + "=" * 70)

if __name__ == "__main__":
    main()
