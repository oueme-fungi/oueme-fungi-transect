#!/bin/bash -l

# Runs bioinformatics pipeline to demultiplex, quality filter, and denoise
# sequences from IonTorrent and PacBio runs in Oueme Fungi transect project

#SBATCH -A snic2018-8-131
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH -J oueme-fungi-transect
#SBATCH -C usage_mail
#SBATCH --mail-type=ALL

# always redo this makefile, in case something (number of cores!) changes.
rm demux.make

module load R/3.5.0 gcc bioinfo-tools blast/2.7.1+ samtools/1.9 ITSx/1.1-beta \
            gnuparallel/20180822 SMRT/5.0.1 cutadapt/1.16 Fastx/0.0.14

# Make sure we are using the versions of packages that we are supposed to use
# make packrat &&

# convert PacBio RSII files to .bam format
# these are single-threaded options
make -j$SLURM_JOB_CPUS_PER_NODE convert-pacbio ccs ion-raw pb-demux # &&

# demultiplex Pacbio BAM files
# this is multithreaded per-file, so run serial make (with parallel tasks)
#make demux-pacbio &&

# demultiplex and quality filter
# For these targets, operations can be done in parallel on many files,
# so run parallel make.
#make -j$SLURM_JOB_CPUS_PER_NODE trim &&

# denoise to find amplicon sequence variants.
# the dada2 library is already multithreaded, so run a serial make.
#make data/demux.counts data/fastq.counts dada taxonomy
