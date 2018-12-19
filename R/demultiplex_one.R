source("install_packages.R")

library(magrittr)
library(tidyverse)
library(ShortRead)

# for scripted use, these are specified in the makefile
if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- file.path(data.dir, "lab_setup")
  in.fastq <- file.path(pacbio.dir, "pb_483_001.reads_of_insert.demux.fastq.gz")
  out.fastq <- str_replace(in.fastq, fixed(".fastq.gz"), ".demux.fastq.gz")
  in.groups <- str_replace(in.fastq, fixed(".demux.fastq.gz"), ".groups")
  platekey <- file.path(lab.dir, "gITS7_platekey.csv")
  dset <- "short"
  plate <- paste0(dset, "_001")
  demux.dir <- file.path(data.dir, "demultiplex")
  plate.dir <- file.path(demux.dir, dset)
  target <- file.path(plate.dir, paste0(plate, "_A1.fasta.gz"))
} else {
  in.groups <- str_subset(prereqs, "\\.groups")
  in.fastq <- str_subset(prereqs, "\\.fastq\\.gz")
}
stopifnot(file.exists(in.groups),
          file.exists(in.fastq))

the_group <- str_replace(target, fixed(".fasta.gz"), "") %>%
  basename
groups <- read_csv(in.groups) %>%
  filter(group == the_group)
if (!dir.exists(dirname(target))) dir.create(dirname(target), recursive = TRUE)
if (file.exists(target)) {
  file.remove(target)
}
file.create(target)
fastq <- FastqStreamer(in.fastq)
while (length(fq <- yield(fastq))) {
  fq %>% magrittr::extract(as.character(.@id) %in% groups$qseqid) %>%
    writeFastq(file = target, mode = "a")
}