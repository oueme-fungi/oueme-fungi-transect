#!/bin/bash -l

# Runs bioinformatics pipeline to demultiplex, quality filter, and denoise
# sequences from IonTorrent and PacBio runs in Oueme Fungi transect project

#SBATCH -A snic2018-8-131
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH -J oueme-fungi-transect

# always redo this makefile, in case something (number of cores!) changes.
rm demux.make

module load R/3.5.0 gcc bioinfo-tools blast samtools ITSx gnuparallel

# demultiplex and quality filter
# For these targets, operations can be done in parallel on many files,
# so run parallel make.
make -j$(nproc) trim 

# denoise to find amplicon sequence variants.
# For this target, many files are processed together, but the dada2 library
# which is used is already multithreaded, so run a serial make.
make data/demux.counts dada data/fastq.counts
