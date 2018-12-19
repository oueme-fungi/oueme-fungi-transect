source("install_packages.R")

library(magrittr)
library(tidyverse)
library(glue)

if (interactive()) {
  base.dir <- str_extract(getwd(), ".+oueme-fungi-transect")
  seq.dir <- file.path(base.dir, "raw_data")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- file.path(data.dir, "lab_setup")
  dataset.file <- file.path(lab.dir, "datasets.csv")
  target <- file.path(base.dir, "demux.make")
}

datasets <- read_csv(dataset.file) %>%
  mutate(rootdir = file.path(seq.dir, Dataset, Seq.Run), 
         file = map2(rootdir,
                     glue("({Seq.Run}_\\d+\\.reads_of_insert\\.fastq\\.gz|IonXpress_\\d+_rawlib.basecaller.bam)"),
                     list.files,
                     recursive = TRUE)) %>%
  unnest(file) %>%
  mutate(Stem = str_replace(file, fixed(".reads_of_insert.fastq.gz"), ""),
         Stem = str_replace(Stem, "_?rawlib\\.basecaller\\.bam", ""),
         FullStem = file.path(rootdir, Stem) %>%
           str_replace(fixed(seq.dir), "$(SEQDIR)"),
         InFile = str_replace(file, fixed(".bam"), ".fastq.gz"),
         Plate = str_match(Stem, glue("{Seq.Run}_(\\d+)$"))[,2],
         Plate = replace_na(Plate, "001"),
         Plate = glue("{Dataset}-{Plate}"))

#Instructions to make the demultiplex key.
#This needs to be done once per pacbio sequence library.
pb.datasets %>%
  glue_data("
    {FullStem}groups : RVARS += in.fwd <- '{FullStem}{Forward}.blast'); |
                           in.rev <- '{FullStem}{Reverse}.blast'); |
                           plate <- '{Plate}';
    $(SEQDIR)/{Stem}groups : {FullStem}{Forward}.blast |
                           {FullStem}{Reverse}.blast |
                           $(TAG_ROOT)/{Forward}.fasta |
                           $(TAG_ROOT)/{Reverse}.fasta |
                           $(LABDIR)/{PlateKey}") %>%
  paste0(collapse = "\n") %>%
  str_replace_all(fixed("|"), "\\") %>%
  cat(file = target)

# List of groups which need to be extracted.
pb.datasets %<>%
  mutate(PlateKey = file.path(lab.dir, PlateKey),
         PlateKey = map(PlateKey, read_csv)) %>%
  unnest(PlateKey) %>%
  mutate(demux.dir = file.path(rootdir, "demultiplex") %>%
           str_replace(fixed(seq.dir), "$(SEQDIR)"),
         OutFile = glue("{file.path(demux.dir, Plate)}_{well}.group.fasta.gz"))

pb.datasets %$%
  paste0(demux.dir, " : ", OutFile, collapse = "\n") %>%
  cat("\n", ., file = target, append = TRUE, sep = "")

pb.datasets %>%
  mutate(Stem = file.path(rootdir, Stem) %>%
           str_replace(fixed(seq.dir), "$(SEQDIR)"),
         FASTQ = glue("{Stem}demux.fasta.gz"),
         group.file = glue("{Stem}groups")) %>%
  glue_data("{OutFile} : {FASTQ} {group.file}") %>%
  paste0(collapse = "\n") %>%
  cat("\n", ., file = target, append = TRUE, sep = "")

ion.datasets <- datasets %>%
  mutate(File = map2(rootdir,
                     "IonXpress_\\d+_rawlib.basecaller.bam",
                     list.files,
                     recursive = TRUE)) %>%
  unnest(File) %>%
  mutate(Stem = str_replace(File, fixed("_rawlib.basecaller.bam"), ""),
         demux.dir = file.path(rootdir, "demultiplex") %>%
           str_replace(fixed(seq.dir), "$(SEQDIR)"),
         tag.fwd = str_sub(Stem, -3, -1) %>%
           as.integer %>%
           paste0("gITS7mod-", .)) %>%
  left_join(read_csv(file.path(lab.dir, unique(.$PlateKey)))) %>%
  mutate(out.file = glue("{Dataset}-001_{well}.fasta.gz") %>%
           file.path(demux.dir, .)) %>%
  glue_data("{OutFile} : ")
  
