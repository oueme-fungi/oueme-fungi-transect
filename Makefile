# Here is my makefile!

SHELL=/bin/bash

# Set up directory and file names
BASEDIR != pwd

DATADIR=${BASEDIR}/data
DEMUXDIR=${DATADIR}/demultiplex
LABDIR=$(DATADIR)/lab_setup
PACBIODATA=${DATADIR}/PacBio
REF_ROOT=${DATADIR}/reference/reference
TAG_ROOT=${DATADIR}/lab_setup/tags
TAGS=gits7 its1 lr5 its4
TAG_FILES=$(addsuffix .fasta,$(addprefix $(TAG_ROOT),$(TAGS)))
ALL_PREFIX=${DATADIR}
DATASET=$(LABDIR)/datasets.csv

GITS7_TAGFILE=$(LABDIR)/Hectors_tag_primer_plates.xlsx
LR5_TAGFILE=$(LABDIR)/Brendan_soil2.xlsx

# Options for R
OPT = --vanilla -q
# default variables for R
RVARS = base.dir <- '${BASEDIR}'; \
        data.dir <- '$(DATADIR)'; \
        lab.dir <- '$(LABDIR)'; \
        gits7.file <- '$(GITS7_TAGFILE)'; \
        its1.lr5.file <- '$(LR5_TAGFILE)'; \
        tags.dir <- '$(TAG_ROOT)'; \
        dataset.file <- '$(DATASET)';
# Command to knit an Rmarkdown file
RMD = cd $(<D) && R $(ROPT) -e "require(rmarkdown); $(RVARS) render('$(<F)', output_format = 'pdf_document', output_dir = '../output')"
# command to run an R script
R = cd $(<D) && R $(ROPT) -e "$(RVARS) target <- '$@'; prereqs <- strsplit('$^', ' '); source('$(<F)', echo = TRUE)" &>"$(notdir $(patsubst %.R,%.Rout,$*.$(<F)))"
# shell commands to convert fastq.gz to fasta
FASTAQ_A = gzip -dc $< | paste - - - - | cut -f 1,2 | tr "\t" "\n" | sed "s/^@/>/"
# columns to use in blasting
TAGBLAST_COLS = qseqid sseqid length qcovs nident pident bitscore evalue sstart send qstart qend
# shell commands to blast a fastq.gz against a tag database
TAGBLAST = $(FASTAQ_A) | blastn \
-task blastn-short \
-db $(TAG_FILE) \
-outfmt '6 $(TAGBLAST_COLS)'\
-culling_limit 10 \
-num_threads 3 \
>$@

# find raw PacBio files, generate target names for processed versions
SHORTPACBIO = $(shell find ${PACBIODATA} -iregex ".*pb_483_[0-9]+\.reads_of_insert\.fastq\.gz")
LONGPACBIO = $(shell find ${PACBIODATA} -iregex ".*pb_500_[0-9]+\.reads_of_insert\.fastq\.gz")
PACBIOFASTA = $(patsubst %.fastq.gz,%.fasta,$(SHORTPACBIO) $(LONGPACBIO))
PACBIOTAGBLAST_GITS7 = $(patsubst %.fastq.gz,%.gits7.blast,$(SHORTPACBIO))
PACBIOTAGBLAST_ITS4 = $(patsubst %.fastq.gz,%.its4.blast,$(SHORTPACBIO))
PACBIOTAGBLAST_LR5 = $(patsubst %.fastq.gz,%.lr5.blast,$(LONGPACBIO))
PACBIOTAGBLAST_ITS1 = $(patsubst %.fastq.gz,%.its1.blast,$(LONGPACBIO))
PACBIOTAGBLAST = $(PACBIOTAGBLAST_GITS7) $(PACBIOTAGBLAST_ITS4) $(PACBIOTAGBLAST_LR5) $(PACBIOTAGBLAST_ITS1)
PACBIOLONGDEMUX  = $(patsubst %.fastq.gz,%.groups,$(LONGPACBIO))
PACBIOSHORTDEMUX = $(patsubst %.fastq.gz,%.groups,$(SHORTPACBIO))
PACBIOTRIMDEMUX = $(patsubst %.fastq.gz,%.demux.trim.fastq.gz,$(LONGPACBIO) $(SHORTPACBIO))


# define paths to look for fast[aq] files
vpath %.fasta ${REF_ROOT}:${PACBIODATA}
vpath %.fastq.gz ${PACBIODATA}
vpath %.R $(BASEDIR)/R

.PHONY: all

all : demultiplex

# create a makefile to demultiplex the files
demux.make : make_demux.R $(DATASET)
	$(R)

include demux.make

# make a blast database from a fasta file
%.nhr %.nin %.nsq : %.fasta
	makeblastdb -in $< -out $* -dbtype nucl

# blast fasta files against the tags databases
%.gits7.blast : TAG_FILE := $(TAG_ROOT)/gits7
%.gits7.blast : %.fastq.gz $(TAG_ROOT)/gits7.nhr $(TAG_ROOT)/gits7.nin $(TAG_ROOT)/gits7.nsq
	$(TAGBLAST)
%.its1.blast : TAG_FILE := $(TAG_ROOT)/its1
%.its1.blast : %.fastq.gz $(TAG_ROOT)/its1.nhr $(TAG_ROOT)/its1.nin $(TAG_ROOT)/its1.nsq
	$(TAGBLAST)
%.lr5.blast : TAG_FILE := $(TAG_ROOT)/lr5
%.lr5.blast : %.fastq.gz $(TAG_ROOT)/lr5.nhr $(TAG_ROOT)/lr5.nin $(TAG_ROOT)/lr5.nsq
	$(TAGBLAST)
%.its4.blast : TAG_FILE := $(TAG_ROOT)/its4
%.its4.blast : %.fastq.gz $(TAG_ROOT)/its4.nhr $(TAG_ROOT)/its4.nin $(TAG_ROOT)/its4.nsq
	$(TAGBLAST)

# demultiplex 

%.groups %.demux.fastq.gz : demultiplex_all.R %.fastq.gz
	$(R)

$(DEMUXDIR)/%.fasta.gz : demultiplex_one.R $(DATASET)
	mkdir -p $(dir $@)
	$(R)

# quality trim

%.trim.fastq.gz %.scrap.fastq.gz : %.fastq.gz
	vsearch --fastq_filter $< --fastq_qmax 50 --fastq_maxee_rate 0.001 --fastq_minlen 150 --fastq_maxns 0 --fastqout_discarded $*.scrap.fastq.gz --fastqout $*.trim.fastq.gz

# blast a fastq.gz file against the reference (taxonomy) database
%.blast : %.fastq.gz $(REF_ROOT).nhr \
          $(REF_ROOT).nin $(REF_ROOT).nsq
	gzip -dc $< | \
	paste - - - - | \
	cut -f 1,2 | \
	tr "\t" "\n" | \
	sed "s/^@/>/" | \
	$(TAGBLAST)

# convert a fastq.gz file to fasta; (don't) rename the sequences to be the file name plus an index
#%.fasta : %.fastq.gz
#	gzip -dc $< | \
#	paste - - - - | \
#	cut -f 1,2 | \
#	tr "\t" "\n" | \
#	sed "s/^@/>/" >$@
#awk 'BEGIN {x=1} /^@/ {printf(">$(basename $(basename $(@F)))-%08d\n", x++)} !/^@/ {print}'

# generate fasta files of the (tagged) primers
$(TAG_ROOT)/gits7.fasta : $(GITS7_TAGFILE)
$(TAG_ROOT)/its4.fasta : $(GITS7_TAGFILE)
$(TAG_ROOT)/lr5.fasta : $(LR5_TAGFILE)
$(TAG_ROOT)/its1.fasta : $(LR5_TAGFILE)
$(TAG_ROOT)/%.fasta : %.extract.R
	$(R)

# the script for creating fasta files of the primers is the same for every set, it just looks at its own name to know what to do.
$(BASEDIR)/R/%.extract.R : tags.extract.R
	cd $(BASEDIR)/R && ln -s $< $@

.PRECIOUS : %.fasta %.nin %.nhr %.nsq