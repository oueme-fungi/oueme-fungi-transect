#!/bin/bash -l

# Runs bioinformatics pipeline to make region-specific ribsomal RNA ITS and LSU databases

#SBATCH -A snic2018-8-131
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -J oueme-fungi-transect
#SBATCH -C usage_mail
#SBATCH --mail-type=ALL
#SBATCH --output="logs/snakemake-%j.log"
#SBATCH --error="logs/snakemake-%j.log"

snakemake -pr --jobs $SLURM_JOB_CPUS_PER_NODE\
  --keep-going\
  --use-conda\
  --shadow-prefix /scratch\
  reference/unite.ITS.fasta.gz\
  reference/unite.ITS1.fasta.gz\
  reference/unite.ITS2.fasta.gz\
  reference/rdp.LSU.fasta.gz\
  reference/silva.LSU.fasta.gz
