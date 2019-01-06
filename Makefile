# Here is my makefile!

SHELL=/bin/bash

# Set up directory and file names
BASEDIR := $(shell pwd)

DATADIR := ${BASEDIR}/data
SEQDIR := ${BASEDIR}/raw_data
LABDIR := $(DATADIR)/lab_setup
PACBIODATA := ${DATADIR}/PacBio
REF_ROOT := ${DATADIR}/reference/reference
TAG_ROOT := ${LABDIR}/tags
TAGS := gits7 gits7_ion its1 lr5 its4
TAG_FILES := $(addsuffix .fasta,$(addprefix $(TAG_ROOT),$(TAGS)))
ALL_PREFIX := ${DATADIR}
DATASET := $(LABDIR)/datasets.csv
PREREQLIST = $^

GITS7_TAGFILE := $(LABDIR)/Hectors_tag_primer_plates.xlsx
LR5_TAGFILE := $(LABDIR)/Brendan_soil2.xlsx

# Options for R
ROPT := --no-save --no-restore -q
# default variables for R
RVARS = options(repos = c(CRAN = 'https://ftp.acc.umu.se/mirror/CRAN/')); \
        base.dir <- '${BASEDIR}'; \
        data.dir <- '$(DATADIR)'; \
        seq.dir <- '$(SEQDIR)'; \
        lab.dir <- '$(LABDIR)'; \
        gits7.file <- '$(GITS7_TAGFILE)'; \
        its1.lr5.file <- '$(LR5_TAGFILE)'; \
        tags.dir <- '$(TAG_ROOT)'; \
        dataset.file <- '$(DATASET)'; \
        stem <- '$*'; \
        target <- '$@'; \
        prereqs <- unlist(strsplit('$(PREREQLIST)', ' ')); \
        splits <- unlist(strsplit('$(SPLITS)', ' '));
# Command to knit an Rmarkdown file
RMD = cd $(<D) && R $(ROPT) -e "require(rmarkdown); $(RVARS) render('$(<F)', output_format = 'pdf_document', output_dir = '../output')"
# command to run an R script
R = cd $(<D) && R $(ROPT) -e "$(RVARS) source('$(<F)', echo = TRUE)" &>"$(patsubst %.R,%.Rout,$@.$(<F))"
# shell commands to convert fastq.gz to fasta
FASTAQ_A = gzip -dc $< | paste - - - - | cut -f 1,2 | tr "\t" "\n" | sed "s/^@/>/"
# columns to use in blasting
TAGBLAST_COLS := qseqid sseqid length qcovs nident pident bitscore evalue sstart send qstart qend
# shell commands to blast a fastq.gz against a tag database
TAGBLAST = $(FASTAQ_A) | blastn \
-task blastn-short \
-db $(TAG_FILE) \
-outfmt '6 $(TAGBLAST_COLS)'\
-culling_limit 10 \
-num_threads 3 \
>$@


# define paths to look for fast[aq] files
vpath %.fasta ${REF_ROOT}:${TAGROOT}
vpath %.R $(BASEDIR)/R

.PHONY: all demultiplex dada trim

all: trim

# Create a makefile to demultiplex the files.
# This will read the plate definitions
# and look in the data directory to see what raw files are present
# and how they need to be demultiplexed, so that we don't need to specify
# every file by hand here.
demux.make: make_demux.R $(DATASET)
	$(R)

include demux.make

# Rule to make a .fasta.gz from a .bam
%.fastq.gz: %.bam
	samtools fastq $< -0 $@

# How many cores do we have?
# Note: this is used to decide how many pieces to split files into for
# demultiplexing. It does not in itself determine how many cores will be used
# in the pipeline.  That is determined by the -j flag to make! However, in
# general, this will allow all cores to be usefully utilized.
NCORES := $(shell nproc)
SPLITS := $(shell echo x{a..z}{a..z} | cut -d' ' -f 1-$(NCORES))

# Rules to split a fasta.gz into one file per processor
# .split is an empty flag file. Its recipe generates all the other files.
%.split : %.fastq.gz
	touch $@.prestamp
	zcat $< | \
	paste - - - - | \
	split -n r/4 --filter 'sed "s/\t/\n/g" | gzip -c >$*-$$FILE.fastq.gz'
	mv -f $@.prestamp $@
# each split file depends on the index
define SPLITRECIPE =
%-$(1).fastq.gz : %.split ;

endef

$(foreach INDEX,$(SPLITS),$(eval $(call SPLITRECIPE,$(INDEX))))

# rule to put the demultiplexed files back together
UNSPLIT = zcat $^ | gzip -c >$@ && rm -f $^

# Rule to make a blast database from a fasta file
%.nhr %.nin %.nsq: %.fasta
	makeblastdb -in $< -out $* -dbtype nucl

# Rule to blast a fastq.gz file against the reference (taxonomy) database
%.blast: %.fastq.gz $(REF_ROOT).nhr \
          $(REF_ROOT).nin $(REF_ROOT).nsq
	$(TAGBLAST)

# Rule to generate fasta files of the (tagged) primers from the various files
# they are stored in.
$(TAG_ROOT)/gits7.fasta: $(GITS7_TAGFILE)
$(TAG_ROOT)/its4.fasta: $(GITS7_TAGFILE)
$(TAG_ROOT)/gits7_ion.fasta: $(GITS7_TAGFILE)
$(TAG_ROOT)/lr5.fasta: $(LR5_TAGFILE)
$(TAG_ROOT)/its1.fasta: $(LR5_TAGFILE)
$(TAG_ROOT)/%.fasta: %.extract.R
	$(R)

# Rule to filter and trim a (demultilpexed) fastq file
# The name of the source file is generated by make.demux.R and found in
# demux.make (in order to place them in seperate directories).

%.trim.fastq.gz : quality_trim.R
	mkdir -p $(@D)
	$(R)

# Rule to find amplicon sequence variants (ASV).
# source files needed to be added as additional prerequisites
# (done in demux.make)
# The prereq list gets too long to pass to R for this command, so we override
# it with just the name of the directory.
%.asv.RDS : PREREQLIST=$(sort $(dir $^))
%.asv.RDS : dada.R
	$(R)

.INTERMEDIATE: $(TAG_ROOT)/gits7 \
               $(TAG_ROOT)/its4 \
               $(TAG_ROOT)/gits7_ion \
               $(TAG_ROOT)/lr5 \
               $(TAG_ROOT)/its1

$(TAG_ROOT)/%: $(TAG_ROOT)/%.nsq \
               $(TAG_ROOT)/%.nhr \
               $(TAG_ROOT)/%.nin ;

# the script for creating fasta files of the primers is the same for every set,
# it just looks at its own name to know what to do.
%.extract.R: tags.extract.R
	cd $(BASEDIR)/R && ln -srf $(<F) $(@F)

.PRECIOUS: %.fasta %.nin %.nhr %.nsq