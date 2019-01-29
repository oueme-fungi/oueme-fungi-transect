###############################################################################
# Makefile to perform analysis for Oueme soil fungi transect project
# 
# Including base calling on raw RSII output files, demultiplexing,
# quality filtering, denoising with dada2, and assigning sequences to
# taxonomy and guilds.
# 
# Early processing steps rely heavily on command-line tools, which are called
# directly from this makefile, with some help from demux.make, which includes
# dependencies which are created based on the contents of datasets.csv.
# 
# Later processing steps are performed in R, and job scheduling is performed
# by the drake package.
# 
###############################################################################

####################### parallel processing setup #############################
# Check if we are running on a SLURM cluster, and set up the number of cores
# to use accordingly.
ifdef SLURM_JOB_ID
$(info Using SLURM.)
SHELL=srun
SRUNFLAGS=-n 1
.SHELLFLAGS= -c $(CORES_PER_TASK) $(SRUNFLAGS) --exclusive /bin/bash -c
NCORES := $(SLURM_JOB_CPUS_PER_NODE)# assume we have the whole node.
else
$(info Not using SLURM.)
SHELL=/bin/bash -c
# If we're not doing this on the cluster (i.e., testing),
# be nice and leave one extra CPU
NCORES := $(shell expr $$(nproc) - 1)
endif
# This is a default value that we will override when we do a build step that is
# inherently multithreaded.
CORES_PER_TASK := 1

# For "embarrasingly" parallel computations, we need to split the input files
# into one piece per core.  These will be labeled by the three-character strings
# in SPLITS (which mimic the default labelling scheme of the "split"" utility)
export SPLITS := $(shell echo x{a..z}{a..z} | cut -d' ' -f 1-$(NCORES))
$(info ncores: $(NCORES))
$(info splits: $(SPLITS))

###################### Set up directory and file names ########################
export BASEDIR := $(shell pwd)

RDIR := ${BASEDIR}/R# R scripts
vpath %.R $(BASEDIR)/R

export RAWDIR := ${BASEDIR}/raw_data# Data received from the sequencing center

SEQDIR := ${BASEDIR}/sequences# sequence data produced by this pipeline
MOVIEDIR := ${SEQDIR}/rawmovie# PacBio movies in *.bam format
DEMUXMOVIEDIR := ${SEQDIR}/demuxmovie#Demultiplexed PacBio movies
CCSDIR := ${SEQDIR}/ccs#Circular consensus BAM files
DEMUXDIR := ${SEQDIR}/demux#Demultiplexed .fastq.gz files
TRIMDIR := ${SEQDIR}/trim
vpath %.trim.fastq.gz $(TRIMDIR)

DATADIR := ${BASEDIR}/data# non-sequence data
export LABDIR := $(DATADIR)/lab_setup# data defining experimental setup
export TAG_ROOT := ${LABDIR}/tags
export DATASET := $(LABDIR)/datasets.csv
export GITS7_TAGFILE := $(LABDIR)/Hectors_tag_primer_plates.xlsx
export LR5_TAGFILE := $(LABDIR)/Brendan_soil2.xlsx

REF_ROOT := ${BASEDIR}/reference# reference databases
vpath %_ref.fasta.gz $(REF_ROOT)

########################## high-level targets ##########################
.PHONY: all demultiplex dada trim convert-pacbio taxonomy clean
all: taxonomy

############################## Running R scripts #############################
# Options for R
ROPT := --no-save --no-environ --no-restore
RARGS := --bootstrap-packrat
# R scripts will use this to know what file(s) they are supposed to be writing.
export TARGETLIST = $@

# Command to knit an Rmarkdown file
RMD = cd $(<D) &&\
      R $(ROPT) -e "require(rmarkdown);\
      $(RVARS)\
      render('$(<F)', output_format = 'pdf_document', output_dir = '../output')"
# command to run an R script
# The list of prerequisites is passed on stdin.  Not all scripts use this.
R = cd $(<D) &&\
    echo $^ |\
    Rscript $(ROPT) $(<F) $(RARGS) &>"$(patsubst %.R,%.Rout,$@.$(<F))"

# R package management is done by packrat.  Before we start any R scripts, we
# should make sure the packrat library is current.
.packrat: packrat/packrat.lock
	Rscript $(ROPT) -e 'packrat::restore()' --args $(RARGS) &> packrat.Rout
	touch $@

#######################  basecalling PacBio RSII files ########################

# Find all raw PacBio RSII files
PB.h5=$(shell find $(RAWDIR) -name *.bas.h5)

# What directories are they in?
PB.PATH=$(sort $(dir $(PB.h5)))
vpath %.h5 $(subst  ,:,$(PB.PATH))

# Get the names of the movie files
PB.movies=$(sort $(notdir $(PB.h5)))

# And the names of the plates (these will be in the paths)
PB.plates=$(sort $(shell echo $(PB.h5) |\
  grep -o -E 'pb_[0-9]{3}_[0-9]{3}'))
$(info PacBio plate names:)
$(info $(PB.plates))

# Make a pacbio .bam from older .h5 files
$(MOVIEDIR)/%.subreads.bam $(MOVIEDIR)/%.scraps.bam : %.bas.h5 %.1.bax.h5 %.2.bax.h5 %.3.bax.h5
	mkdir -p $(@D)
	bax2bam $(filter %.bax.h5,$^) -o $(MOVIEDIR)/$*

# Demultiplex bam files
$(DEMUXMOVIEDIR)/%.demux.subreads.bam\
$(DEMUXMOVIEDIR)/%.demux.scraps.bam: CORES_PER_TASK := $(NCORES)
$(DEMUXMOVIEDIR)/%.demux.subreads.bam\
$(DEMUXMOVIEDIR)/%.demux.scraps.bam: $(MOVIEDIR)/%.subreads.bam\
                                     $(MOVIEDIR)/%.scraps.bam
	mkdir -p $(@D)
	bam2bam $(filter %.subreads.bam,$^) \
	        $(filter %.scraps.bam,$^) \
	        -j $(NCORES) \
	        -o $*.demux \
	        --barcodes=$(filter %.fasta,$^) \
	        --scoremode=asymmetric

# The circular consensus sequence files should have the whole plate in
# one file, named after the plate.

# Function to find the movie names that match a certain plate
matchplates=$(shell echo $(PB.h5) |\
  tr " " "\n" |\
  sed -n -r '/$(1)/ { s@.*/([^/]+)\.bas\.h5@$(DEMUXMOVIEDIR)/\1.demux.subreads.bam@ p}' |\
  tr "\n" " ")
# Recipe to make a CCS for all the reads from one plate,
# using the demultiplexed BAM files for that plate
vpath %.ccs.bam $(CCSDIR)
define CCSRULE=
$(CCSDIR)/$1.ccs.bam: CORES_PER_TASK := $(NCORES)
$(CCSDIR)/$1.ccs.bam: $(call matchplates,$(1))
	mkdir -p $(@D)
	ccs $$@ $$+

endef
# Make all the plates
$(foreach plate,$(PB.plates),\
  $(eval $(call CCSRULE,$(plate))))

convert-pacbio : $(addsuffix .ccs.bam,$(PB.plates))

# Create a makefile to demultiplex the files.
# This will read the plate definitions
# and look in the data directory to see what raw files are present
# and how they need to be demultiplexed, so that we don't need to specify
# every file by hand here.
demux.make: make_demux.R $(DATASET) .packrat
	$(R)

include demux.make

########################## Bioinformatics recipes ##############################

# make a .fastq.gz from a .bam
%.fastq.gz: %.bam
	samtools fastq $< -0 $@

# make a .fasta from a .fastq.gz
%.fasta: %.fastq.gz
	zcat $< | paste - - - - | cut -f 1,2 | tr "\t" "\n" | sed "s/^@/>/" >$@

# use ITSx to find ITS and LSU sequences
%.positions.txt: %.fasta
	ITSx -t f -i $< --summary F --fasta F --graphical F \
	--not-found F --complement F --save_regions "" -o $*

# split a fastq.gz file into ITS1, ITS2, and LSU
%.ITS1.fastq.gz %.ITS2.fastq.gz %.LSU.fastq.gz: extract_regions.R\
                                                %.fastq.gz\
                                                %.positions.txt\
                                                .packrat
	$(R)

# Split a fastq.gz into one file per processor.
# .split is an empty flag file. Its recipe generates all the other files.
%.split: %.fastq.gz
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

# Put the demultiplexed files back together
# This recipe line is used in demux.make.
UNSPLIT = zcat $^ | gzip -c >$@ && rm -f $^

# Rule to make a blast database from a fasta file
vpath %.fasta ${TAGROOT}
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
$(TAG_ROOT)/%.fasta: tags.extract.R .packrat
	$(R)

# Recipe to demultiplex a fastq.gz file.
# Called in rules in demux.make
define DEMUX=
	[ -e $(@D) ] || mkdir -p $(@D)
	rm -f $@.prestamp
	touch $@.prestamp
	$(R)
	mv -f $@.prestamp $@

endef

# Rule to filter and trim a (demultiplexed) fastq file
# The name of the source file is generated by make.demux.R and found in
# demux.make.

$(TRIMDIR)/%.trim.fastq.gz: quality_trim.R $(DEMUXDIR)/%.fastq.gz .packrat
	mkdir -p $(@D)
	$(R)

# Rule to find amplicon sequence variants (ASV).
# source files needed to be added as additional prerequisites
# (done in demux.make)
.INTERMEDIATE: %.dada.Rdata
%.dada.Rdata: CORES_PER_TASK := NCORES
%.dada.Rdata: dada.R .packrat
	$(R)

# Rule to assign taxonomy to ASVs
# The correct reference file needs to be added as a separate rule.
%.dada.taxonomy.rds %.dada.taxonomy.csv: CORES_PER_TASK := NCORES
%.dada.taxonomy.rds %.dada.taxonomy.csv: assign_taxonomy.R %.dada.nochim.rds .packrat
	$(R)

.INTERMEDIATE: $(TAG_ROOT)/gits7 \
               $(TAG_ROOT)/its4 \
               $(TAG_ROOT)/gits7_ion \
               $(TAG_ROOT)/lr5 \
               $(TAG_ROOT)/its1

$(TAG_ROOT)/%: $(TAG_ROOT)/%.nsq \
               $(TAG_ROOT)/%.nhr \
               $(TAG_ROOT)/%.nin ;

.PRECIOUS: %.fasta %.nin %.nhr %.nsq

# look through the log files for demultiplexing to find out how many sequences
# were present at each stage
data/demux.counts: demultiplex
	find $(DEMUXDIR) -name *demultiplex_all.Rout | \
	xargs awk '/sequences/ {print FILENAME " : " $$0 }' >$@


# count sequences in fastq.gz files generated at different stages.
# files to count are added as prerequisites in demux.make
define filecho =
@echo $(1) >>$@.temp

endef

# count the number of sequences in different fastq files
data/fastq.counts: CORES_PER_TASK := NCORES
data/fastq.counts:
	rm -f $@.temp
	touch $@.temp
	$(foreach f,$^,$(call filecho,$f))
	parallel 'echo {}, $$(zcat {} | grep -c "^@")' <$@.temp >$@
	rm $@.temp

# dada.R produced a map from sequences to ASVs as well as an ASV matrix
# for the community.  Because it's easiest if the make target has only one output,
# it puts them into a single Rdata file.  This splits them into two .rds files.
%.dada.dadamap.rds %.dada.seqtable.rds %.dada.nochim.rds: split_rdata.R %.dada.Rdata .packrat
	$(R)

.PHONY: clean
clean:
	for f in $$(find $(SEQDIR) -name demultiplex -type d); do rm -r $$f; done
	for f in $$(find $(SEQDIR) -name trim -type d); do rm -r $$f; done
	for f in $$(find $(SEQDIR) -name *-x??.fastq.gz); do rm $$f; done
	for f in $$(find . -name *.Rout); do rm $$f; done
	for f in $$(find $(DATADIR) -name *.dada.Rdata); do rm $$f; done
	for f in $$(find $(DATADIR) -name *.dada.asv.rds); do rm $$f; done
	for f in $$(find $(DATADIR) -name *.dada.derep.rds); do rm $$f; done
	rm -f $(DATADIR)/demux.counts
	rm -f $(DATADIR)/fastq.counts
