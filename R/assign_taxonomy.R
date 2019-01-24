library(magrittr)
library(tidyverse)
library(glue)
library(dada2)

if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  ref.dir <- file.path(base.dir, "reference")
  lab.dir <- file.path(data.dir, "lab_setup")
  seq.dir <- file.path(base.dir, "raw_data")
  
  # choose a dataset and run for testing.
  dataset <- "long-pacbio"
  seq.run <- "pb_500"
  region = ".LSU"
  in.file <- glue("{file.path(data.dir, dataset)}_{seq.run}{region}.dada.seqtable.rds")
  target <- str_replace(in.file, fixed(".seqtable.rds"), ".taxonomy.rds")
  reference <- file.path(ref.dir, "lsu_ref.fasta.gz")
} else {
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  prereqs <- str_split(prereqs, " ") %>% unlist
  in.file <- str_subset(prereqs, fixed(".nochim.rds"))
  reference <- str_subset(prereqs, fixed(".fasta.gz"))
}


cat("prereqs:", prereqs,
  "\nin.file:", in.file,
  "\nreference:", reference,
  "\n")

stopifnot(file.exists(in.file),
          file.exists(reference))
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
