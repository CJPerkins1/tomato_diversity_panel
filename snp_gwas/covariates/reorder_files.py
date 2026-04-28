#!/usr/bin/env python3

# Read FAM file to get correct order
fam_order = []
with open('/xdisk/yadegari/cjperkins1/snp_gwas/data/concordant_gwas.fam', 'r') as f:
    for line in f:
        fam_order.append(line.split()[1])  # Second column

# Read phenotype file into dictionary
pheno_dict = {}
with open('gemma_phenotypes_unordered.txt', 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        sample_id = parts[1]  # Second column is IID
        pheno_dict[sample_id] = line.strip()

# Read covariate file into dictionary
covar_dict = {}
with open('gemma_covariates_unordered.txt', 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        sample_id = parts[1]
        covar_dict[sample_id] = line.strip()

# Write ordered files
with open('gemma_phenotypes.txt', 'w') as f:
    for sample_id in fam_order:
        if sample_id in pheno_dict:
            f.write(pheno_dict[sample_id] + '\n')

with open('gemma_covariates.txt', 'w') as f:
    for sample_id in fam_order:
        if sample_id in covar_dict:
            f.write(covar_dict[sample_id] + '\n')

print("Done! Files reordered.")
