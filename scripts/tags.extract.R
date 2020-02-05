# for scripted use, these are specified in the makefile
if (interactive()) {
  library(here)
  data.dir <- here("data")
  lab.dir <- here("config")
  gits7.file <- file.path(lab.dir, "Hectors_tag_primer_plates.xlsx")
  its1.lr5.file <- file.path(lab.dir, "Brendan_soil2.xlsx")
  tags.dir <- file.path(lab.dir, "tags")
  dataset.file <- file.path(lab.dir, "datasets.csv")
} else if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  lab.dir <- snakemake@config$labdir
  tags.dir <- snakemake@config$tagdir
  gits7.file <- snakemake@input$gits7_tags
  its1.lr5.file <- snakemake@input$lr5_tags
  dataset.file <- snakemake@input$dataset
} else {
  lab.dir <- Sys.getenv("LABDIR")
  tags.dir <- Sys.getenv("TAG_ROOT")
  gits7.file <- Sys.getenv('GITS7_TAGFILE')
  its1.lr5.file <- Sys.getenv('LR5_TAGFILE')
  dataset.file <- Sys.getenv("DATASET")
}

library(magrittr)
library(tidyverse)
library(readxl)
library(seqinr)
library(assertthat)

assert_that(file.exists(gits7.file),
            file.exists(its1.lr5.file),
            file.exists(dataset.file))

tags = list()

# read the gITS7 tags
tags$gits7_tag <- read_xlsx(gits7.file, skip = 1) %>%
  select(name = oligoname, object = sequence) %>%
  filter(str_detect(name, "gITS7mod"))

# read the gITS7 tags, and remove sequencing adapter
ionA <- read_xlsx(gits7.file, skip = 1) %>%
  filter(oligoname == "Ion-A") %$%
  sequence
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
  select(name = oligoname, object)

# read the LR5 tags
tags$lr5_tag <- read_xlsx(its1.lr5.file, sheet = "Taggar ITS1 and LR5", range = "J2:M10") %>%
  dplyr::rename(primer = `Reverse primer`) %>%
  mutate(object = paste0(pad, barcode, primer)) %>%
  select(name = oligoname, object)

# function to take the reverse complement of DNA sequence(s), represented as a string.
revcomp <- function(s)
  map_chr(s, ~ c2s(rev(comp(s2c(.), ambiguous = TRUE, forceToLower = FALSE))))

# read the dataset definitions
dataset <- read_csv(dataset.file, col_types = "ccccicccccicc") %>%
  mutate(
    tagfile_name = file.path(tags.dir, paste0(seq_run, ".fasta")),
    forward_name = file.path(tags.dir, paste0(seq_run, "_R1.fasta")),
    reverse_name = file.path(tags.dir, paste0(seq_run, "_R2.fasta")),
    plate_key = map(
      file.path(lab.dir, plate_key),
      read_csv,
      col_types = cols(
        row = col_integer(),
        col = col_integer(),
        .default = col_character()
      )
    ),
    forward = tags[forward],
    reverse = tags[reverse],
    out.fasta = pmap(
      list(plate_key, forward, reverse),
      function(plate_key, forward, reverse)
        crossing(select(forward, tag_fwd = name, seq_fwd = object),
                 select(reverse, tag_rev = name, seq_rev = object)) %>%
        mutate_at("seq_rev", revcomp) %>%
        left_join(select(plate_key, starts_with("tag"), well)) %>%
        transmute(name = replace_na(well, "unnamed"),
                  object = paste0("^", seq_fwd, "...", seq_rev))
    ),
    reverse.fasta = pmap(
      list(plate_key, forward, reverse),
      function(plate_key, forward, reverse)
        crossing(select(forward, tag_fwd = name, seq_fwd = object),
                 select(reverse, tag_rev = name, seq_rev = object)) %>%
        mutate_at("seq_fwd", revcomp) %>%
        left_join(select(plate_key, starts_with("tag"), well)) %>%
        transmute(name = replace_na(well, "unnamed"),
                  object = paste0("^", seq_rev, ";optional...", seq_fwd))
    )
  )
# Write the outputs
pwalk(
  dataset,
  function(tagfile_name, forward_name, reverse_name,
           tech, out.fasta, reverse.fasta, ...) {
    if (tech == "Illumina") {
      write.fasta(
        sequences = as.list(out.fasta$object),
        names = out.fasta$name,
        file.out = forward_name,
        as.string = TRUE
      )
      write.fasta(
        sequences = as.list(reverse.fasta$object),
        names = reverse.fasta$name,
        file.out = reverse_name,
        as.string = TRUE
      )
    } else {
      write.fasta(
        sequences = as.list(out.fasta$object),
        names = out.fasta$name,
        file.out = tagfile_name,
        as.string = TRUE
      )
    }
  }
)
