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

# Disabled -- not sure how to make it work when calling new SLURM jobs from the main one.

# parallel rsync from https://wiki.ncsa.illinois.edu/display/~wglick/2013/11/01/Parallel+Rsync
# seems to speed it up by almost a factor of 10 by using 20 processes.
# more important at the start (so we can find out quickly if we're going to fail)
# than at the end (if the build was successful, another 10 minutes won't hurt anything)
#startdir=$(pwd) &&
#echo "Copying project directory to node scratch drive..." &&
#echo "Copying directory structure..." &&
#time rsync -r -f"+ */" -f"- *" ./ $SNIC_TMP/ &&
#echo "Copying files..." &&
#time find . ! -type d -print0 |
#  xargs -0 -n1 -P$SLURM_JOB_CPUS_PER_NODE -I% \
#  rsync -a % $SNIC_TMP/% &&
#cd $SNIC_TMP &&

# always redo this makefile, in case something (number of cores!) changes.
rm -f demux.make &&

## texlive from conda isn't functional
module load texlive &&
conda activate oueme1 &&

# convert PacBio RSII files to .bam format
# these are single-threaded options
make -j$SLURM_JOB_CPUS_PER_NODE convert-pacbio ccs ion-raw pb-demux ion-trim &&

# finish the analysis in drake.
make drake
out=$?

# copying back should be much faster, so no need to be parallel
# The parallel version also causes loop directory links in directories which are links
# (e.g., raw_data on UPPMAX)
#echo "Copying changes back to project directory..." &&
#rsync -aiu $SNIC_TMP/ $startdir/ &&
#cd $startdir

exit $out
