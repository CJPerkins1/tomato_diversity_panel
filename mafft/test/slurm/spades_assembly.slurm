#!/bin/bash
#SBATCH --job-name=assemble_genome_spades
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u28/cjperkins1/.bashrc

#load module
module load spades

cd /xdisk/rpalaniv/cjperkins1/mafft/test/genomes

spades.py --pe1-1 SRX4183308_SRR7279483_1.fastq.gz --pe1-2 SRX4183308_SRR7279483_2.fastq.gz \
          --pe2-1 SRX4183310_SRR7279481_1.fastq.gz --pe2-2 SRX4183310_SRR7279481_2.fastq.gz \
          -o spades_output
