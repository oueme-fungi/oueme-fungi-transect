#!/bin/bash -l

# Runs bioinformatics pipeline to demultiplex, quality filter, and denoise
# sequences from IonTorrent and PacBio runs in Oueme Fungi transect project

#SBATCH -A snic2018-8-131
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH -J oueme-fungi-transect
#SBATCH -C usage_mail
#SBATCH -M snowy
#SBATCH --mail-type=ALL
#SBATCH --output="logs/snakemake-%j.log"
#SBATCH --error="logs/snakemake-%j.log"

## texlive from conda isn't functional
module load texlive &&

snakemake -pr --jobs $SLURM_JOB_CPUS_PER_NODE\
  --keep-going\
  --use-conda\
  --shadow-prefix /scratch\
  .preDADA

#exit $out
