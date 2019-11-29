#!/bin/bash -l

# Runs bioinformatics pipeline to demultiplex, quality filter, and denoise
# sequences from IonTorrent and PacBio runs in Oueme Fungi transect project

#SBATCH -A snic2018-8-131
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -J oueme-fungi-transect-repair
#SBATCH -C usage_mail
#SBATCH -M snowy
#SBATCH --mail-type=ALL
#SBATCH --output="logs/repair-%j.log"
#SBATCH --error="logs/repair-%j.log"

snakemake --use-conda repair
