#!/bin/bash
#SBATCH --job-name=align_genomes_mafft_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u28/cjperkins1/.bashrc

#load mafft module
module load mafft

#get to working directory
cd /xdisk/rpalaniv/cjperkins1/mafft/test

mkdir -p results

#define the input TSV file containing genome SRR IDs
input_tsv="/xdisk/rpalaniv/cjperkins1/mafft/test/input/genomes.tsv"

# Read the TSV file and process each paired-end read
while IFS=$'\t' read -r srx_id sr1 sr2; do
    echo "Aligning $sr1 and $sr2..."

    # Perform alignment using MAFFT with --merge
    mafft --merge --thread 4 "${sr1}" "${sr2}" > "results/${srx_id}_aligned.fasta"

done < "$input_tsv"

#mafft --auto genome1.fasta genome2.fasta > alignment_output.fasta

