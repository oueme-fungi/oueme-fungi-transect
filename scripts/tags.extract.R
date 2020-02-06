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
tags$ionA <- tibble(name = "Ion-A", object = ionA)
tags$gits7_iontag <- tags$gits7_tag %>%
  mutate_at("object", str_replace, ionA, "")

# read the (non-tagged) gITS7 primer
tags$gits7 <- read_xlsx(gits7.file, skip = 1) %>%
  select(name = oligoname, object = sequence) %>%
  filter(!str_detect(name, "gITS7mod"),
         str_detect(name, "gITS7"))

# drop the gits7 from the ion barcodes
tags$iontag <- mutate_at(tags$gits7_iontag, "object", str_replace, tags$gits7$object, "")

# read the (non-tagged) ITS4 primers
tags$its4 <- read_xlsx(gits7.file, skip = 1) %>%
  select(name = oligoname, object = sequence) %>%
  filter(str_detect(name, "ITS4"), !str_detect(name, "-IT"))

# read the ITS1 tags
its1tag_raw <- read_xlsx(its1.lr5.file, sheet = "Taggar ITS1 and LR5", range = "B2:E14")

# whole tag, including pad, barcode, and forward primer
tags$its1_tag <- its1tag_raw %>%
  mutate(object = paste0(pad, barcode, `Forward primer`)) %>%
  select(name = oligoname, object)

# just the pad
tags$its1_pad <- its1tag_raw %>%
  select(object = pad) %>%
  unique() %>%
  mutate(name = "pad")

#just the barcode
tags$its1_barcode <- its1tag_raw %>%
  select(name = oligoname, object = barcode)

# just the primer
tags$its1 <- its1tag_raw %>%
  select(object = `Forward primer`) %>%
  unique() %>%
  mutate(name = "its1")

# read the LR5 tags
lr5tag_raw <- read_xlsx(its1.lr5.file, sheet = "Taggar ITS1 and LR5", range = "J2:M10")

# whole tag, including pad, barcode, and reverse primer
tags$lr5_tag <- lr5tag_raw %>%
  mutate(object = paste0(pad, barcode, `Reverse primer`)) %>%
  select(name = oligoname, object)

# just the pad
tags$lr5_pad <- lr5tag_raw %>%
  select(object = pad) %>%
  unique() %>%
  mutate(name = "pad")

#just the barcode
tags$lr5_barcode <- lr5tag_raw %>%
  select(name = oligoname, object = barcode)

# just the primer
tags$lr5 <- lr5tag_raw %>%
  select(object = `Reverse primer`) %>%
  unique() %>%
  mutate(name = "lr5")

# function to take the reverse complement of DNA sequence(s), represented as a string.
revcomp <- function(s)
  map_chr(s, ~ c2s(rev(comp(s2c(.), ambiguous = TRUE, forceToLower = FALSE))))

# read the dataset definitions
dataset <- read_csv(dataset.file, col_types = "cccccccicccccccccccicc") %>%
  pivot_longer(
    forward:primer_reverse, 
    names_to = c("type", ".value"),
    names_pattern = "(.+_)?(forward|reverse)"
  ) %>%
  mutate_at("type", str_replace, "^(.+)_", "_\\1") %>%
  mutate_at("type", replace_na, "") %>%
  mutate(
    file_name = file.path(tags.dir, paste0(seq_run, type, ".fasta")),
    forward_name = file.path(tags.dir, paste0(seq_run, "_R1", type, ".fasta")),
    reverse_name = file.path(tags.dir, paste0(seq_run, "_R2", type, ".fasta")),
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
      function(plate_key, forward, reverse) {
        if (is.null(forward)) forward <- tibble(name = "unnamed", object = "")
        if (is.null(reverse)) {
          reverse <- tibble(name = "unnamed", object = "")
        } else {
          reverse <- mutate_at(reverse, "object", revcomp)
        }
        crossing(select(forward, tag_fwd = name, seq_fwd = object),
                 select(reverse, tag_rev = name, seq_rev = object)) %>%
        left_join(select(plate_key, starts_with("tag"), well)) %>%
        transmute(name = replace_na(well, "unnamed"),
                  object = paste0("^", seq_fwd, "...", seq_rev) %>%
                    str_replace("\\.\\.\\.$", ""))
      }
    ),
    reverse.fasta = pmap(
      list(plate_key, forward, reverse),
      function(plate_key, forward, reverse) {
        
        if (is.null(reverse)) reverse <- tibble(name = "unnamed", object = "")
        if (is.null(forward)) {
          forward <- tibble(name = "unnamed", object = "")
        } else {
          forward <- mutate_at(forward, "object", revcomp)
        }
        crossing(select(forward, tag_fwd = name, seq_fwd = object),
                 select(reverse, tag_rev = name, seq_rev = object)) %>%
        left_join(select(plate_key, starts_with("tag"), well)) %>%
        transmute(name = replace_na(well, "unnamed"),
                  object = paste0("^", seq_rev, "...", seq_fwd) %>%
                    str_replace("^\\^\\.\\.\\.", ""))
      }
    )
  )
# Write the outputs
pwalk(
  dataset,
  function(file_name, forward_name, reverse_name,
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
        file.out = file_name,
        as.string = TRUE
      )
    }
  }
)
