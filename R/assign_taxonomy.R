library(magrittr)
library(tidyverse)
library(glue)
library(assertthat)
library(dada2)

if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  ref.dir <- file.path(base.dir, "reference")
  
  # choose a dataset and run for testing.
  dataset <- "long-pacbio"
  seq.run <- "pb_500"
  region = ".LSU"
  in.file <- glue("{file.path(data.dir, dataset)}_{seq.run}{region}.dada.seqtable.rds")
  target.rds <- str_replace(in.file, fixed(".seqtable.rds"), ".taxonomy.rds")
  target.csv <- str_replace(in.file, fixed(".seqtable.rds"), ".taxonomy.csv")
  reference <- file.path(ref.dir, "lsu_ref.fasta.gz")
} else {
  # Get the necessary information from environmental variables (set by make)
  # or from stdin.
  targets <- Sys.getenv("TARGETLIST") %>% str_split(" ") %>% unlist
  assert_that(length(targets) > 0)
  target.rds <- str_subset(targets, ".*\\.rds")
  target.csv <- str_subset(targets, ".*\\.csv")
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  prereqs <- str_split(prereqs, " ") %>% unlist
  in.file <- str_subset(prereqs, fixed(".nochim.rds"))
  reference <- str_subset(prereqs, fixed(".fasta.gz"))
}

assert_that(file.exists(in.file),
            is.readable(in.file),
            file.exists(reference),
            is.readable(reference),
            is.dir(dirname(target.rds)),
            is.dir(dirname(target.csv)),
            is.writeable(dirname(target.rds)),
            is.writeable(dirname(target.csv)))
seq.table <- readRDS(in.file)

tax <- seq.table %>%
  colnames %>%
  assignTaxonomy(reference, multithread = TRUE) %>%
  as_tibble(rownames = "seq") %>%
  # remove taxon rank prefixed from Unite reference
  mutate_at(vars(-seq), str_replace, "^[kpcofgs]__", "") %>%
  mutate(
    # add Species to RDP reference
    Species = if ("Species" %in% names(.)) Species else NA_character_,
    Species = ifelse(is.na(Genus) | is.na(Species),
                     NA_character_,
                     paste(Genus, Species)))

tax <- seq.table %>%
  t %>%
  as_tibble(rownames = "seq") %>%
  left_join(tax, by = "seq") %>%
  mutate(Taxonomy = paste(Kingdom, Phylum, Class, Order,
                          Family, Genus, Species,
                          sep = ";") %>%
           str_replace_all(fixed(";NA"), ""))

saveRDS(tax, file = target)
write_csv(tax, file = str_replace(target, "\\.rds", ".csv"))
