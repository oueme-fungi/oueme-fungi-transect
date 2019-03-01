#!/bin/bash -l

# Runs bioinformatics pipeline to demultiplex, quality filter, and denoise
# sequences from IonTorrent and PacBio runs in Oueme Fungi transect project

#SBATCH -A snic2018-8-131
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -J oueme-fungi-transect
#SBATCH -C usage_mail
#SBATCH --mail-type=ALL

## texlive from conda isn't functional
module load texlive &&
conda activate oueme1 &&

snakemake --immediate-submit\
  --cluster "sbatch --account {account} --time {time} --cpu {threads} --dependency {dependencies}"\
  --cluster-config config/UPPMAX.yaml\
  --jobscript UPPMAX-jobscript.sh

exit $out
