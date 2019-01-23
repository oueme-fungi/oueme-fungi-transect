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
ROPT := --vanilla --args --bootstrap-packrat
# default variables for R
RVARS = base.dir <- '${BASEDIR}'; \
        data.dir <- '$(DATADIR)'; \
        seq.dir <- '$(SEQDIR)'; \
        lab.dir <- '$(LABDIR)'; \
        gits7.file <- '$(GITS7_TAGFILE)'; \
        its1.lr5.file <- '$(LR5_TAGFILE)'; \
        tags.dir <- '$(TAG_ROOT)'; \
        dataset.file <- '$(DATASET)'; \
        stem <- '$*'; \
        target <- '$@'; \
        splits <- unlist(strsplit('$(SPLITS)', ' '));
# Command to knit an Rmarkdown file
RMD = cd $(<D) && R $(ROPT) -e "require(rmarkdown); $(RVARS) render('$(<F)', output_format = 'pdf_document', output_dir = '../output')"
# command to run an R script
R = cd $(<D) && echo $^ | R -e "source('.Rprofile', echo = TRUE); $(RVARS) source('$(<F)', echo = TRUE)" $(ROPT) &>"$(patsubst %.R,%.Rout,$@.$(<F))"
# shell commands to convert fastq.gz to fasta
FASTQ_A = zcat $< | paste - - - - | cut -f 1,2 | tr "\t" "\n" | sed "s/^@/>/" >$@
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
vpath %.fasta ${TAGROOT}
vpath %_ref.fastq.gz $(REF_ROOT)
vpath %.R $(BASEDIR)/R

# Find all raw PacBio RSII files
RAWPACBIOFILES=$(shell find $(SEQDIR) -name *.bas.h5)
# What directories are they in?
RAWPACBIOPATH=$(dir $(RAWPACBIOFILES))
# names of .bam files that will be generated
RAWPACBIOBAM=$(RAWPACBIOFILES:.bas.h5=.subreads.bam)
vpath %.subreads.bam $(subst  ,:,$(RAWPACBIOPATH))
# Make a pacbio .bam from older .h5 files
%.subreads.bam %.scraps.bam : %.bas.h5 %.1.bax.h5 %.2.bax.h5 %.3.bax.h5
	bax2bam $(filter %.bax.h5,$^)
convert-pacbio : $(RAWPACBIOBAM)
# Demultiplex bam files (CCS)
%.demux.subreads.bam %.demux.scraps.fastq.gz : %.subreads.bam %.scraps.bam
	bam2bam $(filter %.subreads.bam,$^) \
	        $(filter %.scraps.bam,$^) \
	        -j $(NCORES) \
	        -o $*.demux \
	        --barcodes=$(filter %.fasta,$^) \
	        --scoremode=asymmetric

.PHONY: all demultiplex dada trim convert-pacbio taxonomy

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

# make a .fasta from a .fastq.gz
%.fasta : %.fastq.gz
	$(FASTQ_A)

# use ITSx to find ITS and LSU sequences
%.positions.txt : %.fasta
	ITSx -t f -i $< --summary F --fasta F --graphical F \
	--not-found F --complement F --save_regions "" -o $*

# split a fastq.gz file into ITS1, ITS2, and LSU
%.ITS1.fastq.gz %.ITS2.fastq.gz %.LSU.fastq.gz : extract_regions.R %.fastq.gz %.positions.txt
	$(R)

# How many cores do we have?
# Note: this is used to decide how many pieces to split files into before
# demultiplexing. It does not in itself determine how many cores will be used
# in the pipeline.  That is determined by the -j flag to make! However, in
# general, this will allow all cores to be usefully utilized.
NCORES := $(shell nproc)
SPLITS := $(shell echo x{a..z}{a..z} | cut -d' ' -f 1-$(NCORES))
$(info "ncores: $(NCORES)")
$(info "splits: $(SPLITS)")

# Rules to split a fasta.gz into one file per processor
# .split is an empty flag file. Its recipe generates all the other files.
%.split : %.fastq.gz
	touch $@.prestamp
	zcat $< | \
	paste - - - - | \
	split -n r/$(NCORES) --filter 'sed "s/\t/\n/g" | gzip -c >$*-$$FILE.fastq.gz'
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

# Commands to demultiplex a fastq.gz file.  Depends on other information
# defined in demux.make
define DEMUX=
	[ -e $(@D) ] || mkdir -p $(@D)
	rm -f $@.prestamp
	touch $@.prestamp
	$(R)
	mv -f $@.prestamp $@"

endef

# Rule to filter and trim a (demultiplexed) fastq file
# The name of the source file is generated by make.demux.R and found in
# demux.make (in order to place them in seperate directories).

%.trim.fastq.gz : quality_trim.R
	mkdir -p $(@D)
	$(R)

# Rule to find amplicon sequence variants (ASV).
# source files needed to be added as additional prerequisites
# (done in demux.make)
%.dada.Rdata : dada.R
	$(R)

# Rule to assign taxonomy to ASVs
# The correct reference file needs to be added as a separate rule.
%.dada.taxonomy.rds : assign_taxonomy.R %.dada.seqtable.rds
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

# look through the log files for demultiplexing to find out how many sequences
# were present at each stage
data/demux.counts : demultiplex
	for file in $$(find raw_data -name *demultiplex_all.Rout) ; \
	  do cat $$file | grep sequences | sed 's@^@'"$$file"': @'; done >$@

# count sequences in fastq.gz files generated at different stages.
# files to count are added as prerequisites in demux.make
define filecho =
@echo $(1) >>$@.temp

endef

data/fastq.counts :
	rm -f $@.temp
	touch $@.temp
	$(foreach f,$^,$(call filecho,$f))
	parallel 'echo {}, $$(zcat {} | grep -c "^@")' <$@.temp >$@
	rm $@.temp

# dada.R produced a map from sequences to ASVs as well as an ASV matrix
# for the community.  Because it's easiest if the make target has only one output,
# it puts them into a single Rdata file.  This splits them into two .rds files.
%.dada.dadamap.rds %.dada.seqtable.rds : split_rdata.R %.dada.Rdata 
	$(R)

.PHONY: clean
clean :
	for f in $$(find $(SEQDIR) -name demultiplex -type d); do rm -r $$f; done
	for f in $$(find $(SEQDIR) -name trim -type d); do rm -r $$f; done
	for f in $$(find $(SEQDIR) -name *-x??.fastq.gz); do rm $$f; done
	for f in $$(find . -name *.Rout); do rm $$f; done
	for f in $$(find $(DATADIR) -name *.dada.Rdata); do rm $$f; done
	for f in $$(find $(DATADIR) -name *.dada.asv.rds); do rm $$f; done
	for f in $$(find $(DATADIR) -name *.dada.derep.rds); do rm $$f; done
	rm -f $(DATADIR)/demux.counts
	rm -f $(DATADIR)/fastq.counts
