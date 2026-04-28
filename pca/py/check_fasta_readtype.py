import sys
from Bio import SeqIO
import numpy as np

records = list(SeqIO.parse(sys.argv[1], "fasta"))

missing_chars = set('NnXx-?.')

results = []
for r in records:
    if r.id == "REF":
        continue
    seq = str(r.seq).upper()
    n_missing = sum(1 for b in seq if b in missing_chars)
    pct_missing = n_missing / len(seq) * 100
    read_type = "long_fragmented" if "150bp" in r.id else "short"
    results.append((r.id, pct_missing, read_type))

short = [r for r in results if r[2] == "short"]
long_frag = [r for r in results if r[2] == "long_fragmented"]

bins = [(0,20), (20,50), (50,80), (80,99), (99,100), (100,101)]
bin_labels = ["0-20%", "20-50%", "50-80%", "80-99%", "99-100%", "100% (empty)"]

def print_distribution(samples, label):
    pcts = [s[1] for s in samples]
    print(f"\n=== {label} (n={len(samples)}) ===")
    print(f"  Mean:   {np.mean(pcts):.1f}%")
    print(f"  Median: {np.median(pcts):.1f}%")
    print(f"  Max:    {np.max(pcts):.1f}%")
    print(f"  Min:    {np.min(pcts):.1f}%")
    print(f"\n  Missing data distribution:")
    for (lo, hi), lbl in zip(bins, bin_labels):
        count = sum(1 for p in pcts if lo <= p < hi)
        bar = '#' * count
        print(f"    {lbl:<12} {count:>3} samples  {bar}")
    print(f"\n  Worst 5 samples:")
    for sid, pct, _ in sorted(samples, key=lambda x: x[1], reverse=True)[:5]:
        print(f"    {sid:<55} {pct:.1f}%")
    print(f"\n  Best 5 samples:")
    for sid, pct, _ in sorted(samples, key=lambda x: x[1])[:5]:
        print(f"    {sid:<55} {pct:.1f}%")

print_distribution(short, "Short reads")
print_distribution(long_frag, "Long reads (fragmented to 150bp)")

print(f"\n=== Overall (excluding REF) ===")
all_pcts = [r[1] for r in results]
for (lo, hi), lbl in zip(bins, bin_labels):
    count = sum(1 for p in all_pcts if lo <= p < hi)
    short_count = sum(1 for r in results if lo <= r[1] < hi and r[2] == "short")
    long_count  = sum(1 for r in results if lo <= r[1] < hi and r[2] == "long_fragmented")
    bar = '#' * count
    print(f"  {lbl:<12} {count:>3} samples  (short={short_count}, long={long_count})  {bar}")
