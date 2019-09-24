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

for x in $(grep -l '"incomplete": true' .snakemake/metadata/*)
do 
 echo "removing file $x"
 rm $x
done

snakemake -pr --jobs 999\
  --keep-going\
  --profile config/UPPMAX\
  --use-conda\
  --immediate-submit\
  --notemp\
  all_taxonomy

#exit $out
