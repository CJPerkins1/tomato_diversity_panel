import sys
from Bio import SeqIO
import numpy as np

records = list(SeqIO.parse(sys.argv[1], "fasta"))

print(f"Total samples: {len(records)}")
print(f"Alignment length: {len(records[0].seq)} sites\n")

missing_chars = set('NnXx-?.')

results = []
for r in records:
    seq = str(r.seq).upper()
    n_missing = sum(1 for b in seq if b in missing_chars)
    pct_missing = n_missing / len(seq) * 100
    results.append((r.id, n_missing, pct_missing))

results.sort(key=lambda x: x[2], reverse=True)

print("=== Samples with >50% missing data ===")
bad = [r for r in results if r[2] > 50]
if bad:
    for sample_id, n_miss, pct in bad:
        print(f"  {sample_id}: {pct:.1f}% missing ({n_miss} sites)")
else:
    print("  None found.")

print("\n=== Completely empty samples (100% missing) ===")
empty = [r for r in results if r[2] == 100.0]
if empty:
    for sample_id, _, _ in empty:
        print(f"  {sample_id}")
else:
    print("  None found.")

pcts = [r[2] for r in results]
print(f"\n=== Missing data summary across all samples ===")
print(f"  Mean:   {np.mean(pcts):.1f}%")
print(f"  Median: {np.median(pcts):.1f}%")
print(f"  Max:    {np.max(pcts):.1f}%  ({results[0][0]})")
print(f"  Min:    {np.min(pcts):.1f}%")

print(f"\n=== Top 10 samples by missing data ===")
for sample_id, n_miss, pct in results[:10]:
    bar = '#' * int(pct / 2)
    print(f"  {sample_id:<30} {pct:5.1f}%  {bar}")
