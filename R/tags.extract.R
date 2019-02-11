library(magrittr)
#library(plyr)
library(tidyverse)
library(readxl)
library(seqinr)
library(here)
library(assertthat)

# for scripted use, these are specified in the makefile
if (interactive()) {
  data.dir <- here("data")
  lab.dir <- here("config")
  gits7.file <- file.path(lab.dir, "Hectors_tag_primer_plates.xlsx")
  its1.lr5.file <- file.path(lab.dir, "Brendan_soil2.xlsx")
  tags.dir <- file.path(lab.dir, "tags")
  dataset.file <- file.path(lab.dir, "datasets.csv")
} else {
  lab.dir <- Sys.getenv("LABDIR")
  tags.dir <- Sys.getenv("TAG_ROOT")
  gits7.file <- Sys.getenv('GITS7_TAGFILE')
  its1.lr5.file <- Sys.getenv('LR5_TAGFILE')
  dataset.file <- Sys.getenv("DATASET")
}

assert_that(file.exists(gits7.file),
            file.exists(its1.lr5.file),
            file.exists(dataset.file))

tags = list()

# read the gITS7 tags
tags$gits7_tag <- read_xlsx(gits7.file, skip = 1) %>%
  select(name = oligoname, object = sequence) %>%
  filter(str_detect(name, "gITS7mod"))

# read the gITS7 tags, and remove sequencing adapter
ionA <- read_xlsx(gits7.file, skip = 1) %>% filter(oligoname == "Ion-A") %$% sequence
tags$gits7_iontag <- tags$gits7_tag %>%
  mutate_at("object", str_replace, ionA, "")

# read the (non-tagged) gITS7 primer
tags$gits7 <- read_xlsx(gits7.file, skip = 1) %>%
  select(name = oligoname, object = sequence) %>%
  filter(!str_detect(name, "gITS7mod"),
         str_detect(name, "gITS7"))

# read the (non-tagged) ITS4 primers
tags$its4 <- read_xlsx(gits7.file, skip = 1) %>%
  select(name = oligoname, object = sequence) %>%
  filter(str_detect(name, "ITS4"), !str_detect(name, "-IT"))

# read the ITS1 tags
tags$its1_tag <- read_xlsx(its1.lr5.file, sheet = "Taggar ITS1 and LR5", range = "B2:E14") %>%
  dplyr::rename(primer = `Forward primer`) %>%
  mutate(object = paste0(pad, barcode, primer)) %>%
  select(object, name = oligoname)

# read the LR5 tags
tags$lr5_tag <- read_xlsx(its1.lr5.file, sheet = "Taggar ITS1 and LR5", range = "J2:M10") %>%
  dplyr::rename(primer = `Reverse primer`) %>%
  mutate(object = paste0(pad, barcode, primer)) %>%
  select(object, name = oligoname)

# function to take the reverse complement of DNA sequence(s), represented as a string.
revcomp <- function(s)
  map_chr(s, ~ c2s(rev(comp(s2c(.), ambiguous = TRUE, forceToLower = FALSE))))

# read the dataset definitions
dataset <- read_csv(dataset.file) %>%
  mutate(
    tagfile.name = file.path(tags.dir, paste0(Seq.Run, ".fasta")),
    PlateKey = map(file.path(lab.dir, PlateKey), read_csv),
    Forward = tags[Forward],
    Reverse = tags[Reverse],
    out.fasta = pmap(
      list(PlateKey, Forward, Reverse),
      function(PlateKey, Forward, Reverse)
        crossing(select(Forward, tag.fwd = name, seq.fwd = object),
                 select(Reverse, tag.rev = name, seq.rev = object)) %>%
        mutate_at("seq.rev", revcomp) %>%
        left_join(select(PlateKey, starts_with("tag"), well)) %>%
        transmute(name = replace_na(well, "unnamed"),
                  object = paste0(seq.fwd, "...", seq.rev))
    ))
# Write the outputs
dataset %>%
  pwalk(function(tagfile.name, out.fasta, ...) {
    write.fasta(sequences = as.list(out.fasta$object),
                names = out.fasta$name,
                file.out = tagfile.name,
                as.string = TRUE)
  })
