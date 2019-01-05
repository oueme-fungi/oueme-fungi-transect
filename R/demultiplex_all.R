source("install_packages.R")

library(magrittr)
library(tidyverse)
# library(multidplyr)
library(readxl)
library(seqinr)
library(glue)
library(Biostrings)
library(ShortRead)
library(rBLAST)

seq_count <- function(data, description) {
  if ("sample" %in% names(data)){
    data %<>% group_by(sample)
    format_string <- "{n} sequences from {sample} {description}."
  } else {
    data %<>% ungroup
    format_string <- "{n} sequences {description}."
  }
  data %>% summarize(n = n_distinct(qseqid)) %>%
    glue_data(format_string) %>%
    paste(collapse = "\n") %>%
    cat("\n")
  return(data)
}

# for scripted use, these are specified in the makefile
if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- file.path(data.dir, "lab_setup")
  seq.dir <- file.path(base.dir, "raw_data")
  dataset <- "short-ion"
  seq.run <- "is_057"
  fullstem <- file.path(seq.dir, dataset, seq.run)
  in.fastq <- file.path(fullstem, "rawdata", "no_bc-subset", "rawlib.basecaller.fastq.gz")
  demux.dir <- file.path(fullstem, "demultiplex")
  
  gits7.file <- file.path(lab.dir, "Hectors tag primer plates.xlsx")
  its1.lr5.file <- file.path(lab.dir, "Brendan soil2.xlsx")
  tag.files <- list.files(file.path(lab.dir, "tags"), ".+\\.fasta", full.names = TRUE)
  
  blastdb.fwd <- file.path(lab.dir, "tags", "gits7_ion")
  blastdb.rev <- file.path(lab.dir, "tags", "its4")
  
  # in.fwd <- str_replace(in.fastq, fixed(".fastq.gz"), ".gits7_ion.blast")
  # in.rev <- str_replace(in.fastq, fixed(".fastq.gz"), ".its4.blast")
  # in.mid <- str_replace(in.fastq, fixed(".fastq.gz"), ".its4.blast")
  #out.fastq <- str_replace(in.fastq, fixed(".fastq.gz"), ".demux.fastq.gz")
  #out.groups <- str_replace(in.fastq, fixed(".fastq.gz"), ".groups")
  platekey.file <- file.path(lab.dir, "gITS7_platekey.csv")
  plate <- "short_ion_001"
  shard <- 'xaa'
} else {
  tag.files <- str_subset(prereqs, "tags.+\\.fasta")
  glue("tag file: {tag.files}") %>% glue_collapse(sep = "\n") %>% cat("\n", ., "\n")
  in.fastq <- str_subset(prereqs, "\\.fastq\\.gz")
  glue("in fastq: {in.fastq}") %>% glue_collapse(sep = "\n") %>% cat("\n", ., "\n")
  platekey.file <- str_subset(prereqs, "_platekey\\.csv")
  glue("platekey: {platekey.file}") %>% glue_collapse(sep = "\n") %>% cat("\n", ., "\n")
}

stopifnot(file.exists(tag.files),
          file.exists(in.fastq),
          file.exists(platekey.file),
          file.exists(paste0(blastdb.fwd, c(".nin", ".nsq", ".nhr"))),
          file.exists(paste0(blastdb.rev, c(".nin", ".nsq", ".nhr"))),
          !exists("in.mid") || file.exists(in.mid))

if (!dir.exists(demux.dir)) {
  dir.create(demux.dir)
}

blast_opts = c("-task blastn-short",
               "-culling_limit 10",
               "-num_threads 3")
blast_cols = "qseqid sseqid length qcovs nident pident bitscore evalue sstart send qstart qend"

tags <- map(tag.files, read.fasta) %>%
  flatten %>%
  {tibble(seq = getSequence(., as.string = TRUE) %>% unlist,
          tag = getName(.),
          length = str_length(seq))} %>%
  unique

platekey <- read_csv(platekey.file) %>%
  mutate(plate = plate) %>%
  unite("group", plate, well, sep = "_") %>%
  mutate(out.file = glue("{file.path(demux.dir, group)}-{shard}.fastq.gz"))

blankread <- ShortReadQ()

for (f in platekey$out.file) {
  if (file.exists(f)) file.remove(f)
  writeFastq(blankread, file = f)
}

fastq <- FastqStreamer(in.fastq)

blastlist <- list()
blastlist$F <- blast(db = blastdb.fwd, type = "blastn")
blastlist$R <- blast(db = blastdb.rev, type = "blastn")

# clust <- multidplyr::get_default_cluster()

# cluster_library(clust, c("stringr", "dplyr", "rBLAST"))

while (length(fq <- yield(fastq))) {
  names(fq@sread) <- fq@id
  cat("Blasting...\n")
  tagblast <- lapply(X = blastlist,
                     FUN = predict,
                     newdata = fq@sread,
                     BLAST_args = blast_opts,
                     custom_format = blast_cols) %>%
    tibble(direction = names(.), data = .) %>%
    unnest(data) %>%
    left_join(tibble(qseqid = as.character(fq@id),
                     qlength = width(fq@sread),
                     idx = seq_along(qseqid)))%>%
    left_join(tags %>%
                select(sseqid = tag, slength = length))
  
cat("Demultiplexing...\n")

tagblast2 <- tagblast  %>%
  # gITS7 and ITS1 are forward primers; LR5 is reverse
  mutate(primer = str_extract(sseqid, "(gITS7|ITS1|LR5|ITS4)")) %>%
  tidyr::separate(sseqid, c("tag", "v"), sep = "_v") %>%
  # partition(qseqid, cluster = clust) %>%
  seq_count("with at least one tag match")%>%
  # remove matches with less than 90% coverage of the tag
  filter(length / slength >= 0.9) %>%
  seq_count("match at least 90% of tag length") %>%
  # filter(length - nident <= 3) %>%
  seq_count("with at most 3 tag mismatches") %>%
  # take away matches to a degenerate variant with lower score
  group_by(qseqid, tag) %>%
  summarize_all(dplyr::first) %>%
  # look at results for forward and reverse seperately
  group_by(qseqid, direction) %>%
  # if there is a perfect match, then it must be right
  filter(pident == 100 | max(pident) < 100) %>%
  # remove matches that are at least 10x worse
  filter(evalue < 10*min(evalue)) %>%
  # remove sequences that still have more than one hit in each direction
  filter(n() == 1) %>%
  seq_count("with no more than one distinct tag in each direction") %>% 
  group_by() %>%
  select(idx, qseqid, tag, direction, sstart:qend, length, qlength, slength) %>%
  # find sequences which are reversed
  mutate(rev.comp = if_else(direction == "F",
                            sstart > send & qend > qstart,
                            sstart < send & qstart < qend)
  ) %>%
  select(-sstart, -send, -slength, -length) %>%
  # collect() %>%
#  {midtags <<- filter(., direction == "M")} %>%
  # put the forward and reverse tag on the same line
  {inner_join(filter(., direction == "F"),
             filter(., direction == "R"),
             by = c("idx", "qseqid", "qlength", "rev.comp"),
             suffix = c(".fwd", ".rev"))} %>%
  # take only sequences with both tags
  # filter(complete.cases(.)) %>%
  seq_count("with tags in both directions") %>%
  # {if (exists("in.mid")) semi_join(., midtags, by = "idx") %>%
  #     seq_count("with a single middle site") else
  #       .} %>%
  # calculate where to trim to remove the tags
  mutate(trimstart = if_else(rev.comp,
                             qlength - qstart.fwd + 2L,
                             qend.fwd + 1L),
         trimend = if_else(rev.comp,
                           qlength - qend.rev,
                           qstart.rev - 1L)) %>%
  # the tags must be in the correct order
  filter(trimend > trimstart) %>%
  seq_count("with tags appearing in correct order") %>%
  select(idx, qseqid, tag.fwd, tag.rev, rev.comp, trimstart, trimend) %>%
  arrange(idx) %>%
  left_join(platekey) %>%
  group_by(group)

revcomps <- filter(tagblast2, rev.comp)$idx
fq[revcomps] <- reverseComplement(fq[revcomps])

do(tagblast2,
   out = {
     fqsub = fq[.$idx]
     fqsub <- narrow(fqsub, start = .$trimstart, end = .$trimend)
     writeFastq(object = fqsub, file = .$out.file[1],
                mode = "a")
   })
}