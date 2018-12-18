source("install_packages.R")

library(magrittr)
library(tidyverse)
library(glue)

if (interactive()) {
  base.dir <- str_extract(getwd(), ".+oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- file.path(data.dir, "lab_setup")
  dataset.file <- file.path(lab.dir, "datasets.csv")
}

datasets <- read_csv(dataset.file) %>%
  mutate(File = map(glue("{Seq.Run}_\\d+\\.reads_of_insert\\.fastq\\.gz"),
                    list.files, path = data.dir, recursive = TRUE)) %>%
  unnest(File) %>%
  mutate(Stem = str_replace(File, fixed("fastq.gz"), ""),
         Plate = str_match(Stem, glue("{Seq.Run}_(\\d+)\\."))[,2],
         Plate = glue("{Dataset}_{Plate}"))


#Instructions to make the demultiplex key.
#This needs to be done once per sequence library.
datasets %>%
  glue_data("
    $(DATADIR)/{Stem}groups : RVARS += in.fwd <- file.path(data.dir, '{Stem}{Forward}.blast'); |
                           in.rev <- file.path(data.dir, '{Stem}{Reverse}.blast'); |
                           in.fastq <- file.path(data.dir, '{Stem}fastq.gz'); |
                           out.fastq <- file.path(data.dir, '{Stem}demux.fastq.gz'); |
                           out.groups <- file.path(data.dir, '{Stem}groups'); |
                           platekey <- file.path(lab.dir, '{PlateKey}'); |
                           tag.files <- file.path(tags.dir, c('{Forward}.fasta',|
                                          '{Reverse}.fasta')); |
                           plate <- '{Plate}';
    $(DATADIR)/{Stem}groups : $(DATADIR)/{Stem}{Forward}.blast |
                           $(DATADIR)/{Stem}{Reverse}.blast |
                           $(TAG_ROOT)/{Forward}.fasta |
                           $(TAG_ROOT)/{Reverse}.fasta |
                           $(LABDIR)/{PlateKey}") %>%
  paste0(collapse = "\n") %>%
  str_replace_all(fixed("|"), "\\") %>%
  cat(file = file.path(base.dir, "demux.make"))

# List of groups which need to be extracted.
# Each of these will be created by a general rule.
datasets %<>%
  mutate(PlateKey = file.path(lab.dir, PlateKey),
         PlateKey = map(PlateKey, read_csv)) %>%
  unnest(PlateKey) %>%
  mutate(OutFile = glue("{file.path(data.dir, 'demultiplex', Dataset, Plate)}_{well}.fasta.gz"))

datasets %$%
  paste0(OutFile, collapse = " \\\n              ") %>%
  paste0("\ndemultiplex : ", .) %>%
  cat(file = file.path(base.dir, "demux.make"), append = TRUE)

datasets %>%
  glue_data("{OutFile} : {file.path(data.dir, Stem, 'demux.fasta.gz')}") %>%
  paste0(collapse = "\n")
