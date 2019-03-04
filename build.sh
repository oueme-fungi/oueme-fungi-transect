#!/bin/bash -l

# Runs bioinformatics pipeline to demultiplex, quality filter, and denoise
# sequences from IonTorrent and PacBio runs in Oueme Fungi transect project

#SBATCH -A snic2018-8-131
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J oueme-fungi-transect
#SBATCH -C usage_mail
#SBATCH --mail-type=ALL
#SBATCH --output="logs/snakemake-%j.log"
#SBATCH --error="logs/snakemake-%j.log"

## texlive from conda isn't functional
module load texlive &&
conda activate oueme1 &&

snakemake -pr --jobs 999\
  --rerun-incomplete\
  --profile config/UPPMAX

#exit $out
