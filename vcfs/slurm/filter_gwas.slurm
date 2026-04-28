#!/usr/bin/env bash
#SBATCH --job-name=prep_gwas_snps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=4:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load plink/1.9

VCF_FILE="/xdisk/yadegari/cjperkins1/dp_vcfs/concordant_snps/merged_concordant_snps.vcf.gz"
OUTPUT_DIR="/xdisk/yadegari/cjperkins1/snp_gwas/data"

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

echo "Creating GWAS SNP set with MAF ≥ 0.01..."

plink --vcf ${VCF_FILE} \
  --double-id \
  --allow-extra-chr \
  --geno 0.2 \
  --maf 0.01 \
  --make-bed \
  --out concordant_gwas

echo "Done! GWAS SNP set ready: ${OUTPUT_DIR}/concordant_gwas.{bed,bim,fam}"
