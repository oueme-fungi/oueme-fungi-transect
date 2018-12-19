source("install_packages.R")

library(magrittr)
library(tidyverse)
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
  in.fastq <- file.path(fullstem, "rawdata", "bc-subset", "IonXpress_001_rawlib.basecaller.fastq.gz")
  
  gits7.file <- file.path(lab.dir, "Hectors tag primer plates.xlsx")
  its1.lr5.file <- file.path(lab.dir, "Brendan soil2.xlsx")
  tag.files <- list.files(file.path(lab.dir, "tags"), ".+\\.fasta", full.names = TRUE)
  
  in.fwd <- str_replace(in.fastq, fixed(".fastq.gz"), ".gits7_ion.blast")
  in.rev <- str_replace(in.fastq, fixed(".fastq.gz"), ".its4.blast")
  # in.mid <- str_replace(in.fastq, fixed(".fastq.gz"), ".its4.blast")
  out.fastq <- str_replace(in.fastq, fixed(".fastq.gz"), ".demux.fastq.gz")
  out.groups <- str_replace(in.fastq, fixed(".fastq.gz"), ".groups")
  platekey <- file.path(lab.dir, "gITS7_platekey.csv")
  plate <- "short_001"
} else {
  tag.files <- str_subset(prereqs, "tags.+\\.fasta")
  in.fastq <- str_subset(prereqs, "\\.fastq\\.gz")
  platekey <- str_subset(prereqs, "_platekey\\.csv")
  out.fastq <- paste0(stem, ".demux.fastq.gz")
  out.groups <- paste0(stem, ".groups")
}

stopifnot(file.exists(tag.files),
          file.exists(in.fastq),
          file.exists(in.fwd),
          file.exists(in.rev),
          file.exists(platekey),
          !exists("in.mid") || file.exists(in.mid))

tags <- map(tag.files, read.fasta) %>%
  flatten %>%
  {tibble(seq = getSequence(., as.string = TRUE) %>% unlist,
          tag = getName(.),
          length = str_length(seq))} %>%
  unique

fastq <- readFastq(in.fastq)

tagblast <-
  # Read the files
  c(F = in.fwd, R = in.rev, M = if (exists("in.mid")) in.mid else NULL) %>%
  map(read_delim,
          delim = "\t",
          col_names = c("qseqid", "sseqid", "length", "qcovs",
                        "nident", "pident", "bitscore", "evalue",
                        "sstart", "send", "qstart", "qend")) %>%
  tibble(direction = names(.), data = .) %>%
  unnest(data) %>%
  left_join(tibble(qseqid = as.character(fastq@id),
                   qlength = width(fastq@sread),
                   idx = seq_along(qseqid)))%>%
  left_join(tags %>%
              select(sseqid = tag, slength = length))

midtags <- tibble()

tagblast2 <- tagblast  %>%
  # gITS7 and ITS1 are forward primers; LR5 is reverse
  mutate(primer = str_extract(sseqid, "(gITS7|ITS1|LR5|ITS4)")) %>%
  tidyr::separate(sseqid, c("tag", "v"), sep = "_v") %>%
  seq_count("with at least one tag match")%>%
  # remove matches with less than 90% coverage of the tag
  filter(length / slength >= 0.9) %>%
  seq_count("match at least 90% of tag length") %>%
  # filter(length - nident <= 3) %>%
  # seq_count("with at most 3 tag mismatches") %>%
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
  ungroup %>%
  select(idx, qseqid, tag, direction, sstart:qend, length, qlength, slength) %>%
  # find sequences which are reversed
  mutate(rev.comp = if_else(direction == "F",
                            sstart > send & qend > qstart,
                            sstart < send & qstart < qend)
  ) %>%
  select(-sstart, -send, -slength, -length) %T>%
  {midtags <<- filter(., direction == "M")} %>%
  # put the forward and reverse tag on the same line
  {full_join(filter(., direction == "F"),
             filter(., direction == "R"),
             by = c("idx", "qseqid", "qlength", "rev.comp"),
             suffix = c(".fwd", ".rev"))} %>%
  # take only sequences with both tags
  filter(complete.cases(.)) %>%
  seq_count("with tags in both directions") %>%
  {if (exists("in.mid")) semi_join(., midtags, by = "idx") %>%
      seq_count("with a single middle site") else
        .} %>%
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
  arrange(idx)

revcomps <- filter(tagblast2, rev.comp)$idx

fastq[revcomps] <- reverseComplement(fastq[revcomps])
fastq <- fastq[tagblast2$idx]
tagblast2 %<>% mutate(newseqid = glue("{in.name}_{idx}"))
fastq@id <- BStringSet(tagblast2$newseqid)


fastq <- narrow(fastq, start = tagblast2$trimstart, end = tagblast2$trimend)

if (file.exists(out.fastq)) file.remove(out.fastq)
writeFastq(fastq, file = out.fastq)

tagblast2 %>%
  select(seqid = newseqid, tag.fwd, tag.rev) %>%
  left_join(read_csv(platekey)) %>%
  mutate(group = paste(plate, well, sep = "_")) %>%
  select(seqid, group) %>%
  write_csv(out.groups)