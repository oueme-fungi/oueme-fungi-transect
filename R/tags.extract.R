library(magrittr)
library(plyr)
library(tidyverse)
library(readxl)
library(seqinr)

# for scripted use, these are specified in the makefile
if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- file.path(data.dir, "lab_setup")
  gits7.file <- file.path(lab.dir, "Hectors_tag_primer_plates.xlsx")
  its1.lr5.file <- file.path(lab.dir, "Brendan_soil2.xlsx")
  tags.dir <- file.path(lab.dir, "tags")
  which.tags = c("gits7", "its1", "lr5", "its5")
} else {
  which.tags = basename(sys.frame(1)$ofile) %>%
    str_extract("(gits7|its1|lr5|its4)")
}

# given a DNA sequence including IUPAC ambiguous bases, give all possible realizations
all_ambiguous_bases <- function(s) {
  if (str_length(s) == 0) return(s)
  base <-
    revalue(substring(s, 1, 1),
            c(R = "GA",
              Y = "TC",
              S = "GC",
              W = "AT",
              K = "GT",
              M = "AC",
              D = "GTA",
              H = "TAC",
              B = "GTC",
              V = "GAC",
              N = "GATC"),
            FALSE)
  base <- strsplit(base, "") %>% unlist
  if (str_length(s) == 1) return(base)
  s <- substring(s, 2)
  s <- expand.grid(base, all_ambiguous_bases(s), stringsAsFactors = FALSE)
  do.call(paste0, s)
}


tags = list()

# read the gITS7 tags
if ("gits7" %in% which.tags) {
  tags$gits7 <- read_xlsx(gits7.file, skip = 1) %>%
    select(name = oligoname, object = sequence) %>%
    filter(str_detect(name, "gITS7mod"))
}
# read the (non-tagged) primers
if ("its4" %in% which.tags) {
  tags$its4 <- read_xlsx(gits7.file, skip = 1) %>%
    select(name = oligoname, object = sequence) %>%
    filter(str_detect(name, "ITS4"), !str_detect(name, "-IT"))
}
# read the ITS1 tags
if ("its1" %in% which.tags) {
  tags$its1 <- read_xlsx(its1.lr5.file, sheet = "Taggar ITS1 and LR5", range = "B2:E14") %>%
    dplyr::rename(primer = `Forward primer`) %>%
    mutate(object = paste0(pad, barcode, primer)) %>%
    select(object, name = oligoname)
}
# read the LR5 tags
if ("lr5" %in% which.tags) {
  tags$lr5 <- read_xlsx(its1.lr5.file, sheet = "Taggar ITS1 and LR5", range = "J2:M10") %>%
    dplyr::rename(primer = `Reverse primer`) %>%
    mutate(object = paste0(pad, barcode, primer)) %>%
    select(object, name = oligoname)
}

# make the ambiguous bases explicit
tags <- map(tags,
            . %>% mutate_at("object", map, all_ambiguous_bases) %>%
              unnest(object) %>%
              group_by(name) %>%
              mutate(v = seq_along(object)) %>%
              ungroup %>%
              tidyr::unite( "name", name, v, sep = "_v"))

tags$primer

# parse them and write output file(s)
if (!dir.exists(tags.dir)) dir.create(tags.dir)
tags %>%
  map(pmap, as.SeqFastadna) %>%
  imap(~ write.fasta(.x, file.out = file.path(tags.dir, paste0(.y, ".fasta")), names = getName(.x)))
