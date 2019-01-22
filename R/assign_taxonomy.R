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
  dataset <- "short-ion"
  seq.run <- "is_057"
  region = ""
  in.file <- glue("{file.path(data.dir, dataset)}_{seq.run}{region}.dada.seqtable.rds")
  target <- str_replace(in.file, fixed(".seqtable.rds"), ".taxonomy.rds")
  reference <- file.path(ref.dir, "sh_general_release_dynamic_s_01.12.2017.fasta.gz")
} else {
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  prereqs <- str_split(prereqs, " ") %>% unlist
  in.file <- str_extract(prereqs, fixed(".seqtable.rds"))
  reference <- str_extract(prereqs, fixed(".fasta.gz"))
}

stopifnot(file.exists(in.file),
          file.exists(reference))
seq.table <- readRDS(in.file)
tax <- seq.table %>%
  colnames %>%
  assignTaxonomy(reference, multithread = TRUE) %>%
  as_tibble(rownames = "seq") %>%
  mutate_at(vars(Kingdom:Species), str_replace, "^[kpcofgs]__", "") %>%
  mutate(Species = ifelse(is.na(Genus) | is.na(Species),
                          NA_character_,
                          paste(Genus, Species)))

seqtab <- makeSequenceTable(asv) %>%
  t %>%
  as_tibble(rownames = "seq") %>%
  left_join(tax, by = "seq") %>%
  mutate(Taxonomy = paste(Kingdom, Phylum, Class, Order,
                          Family, Genus, Species,
                          sep = ";") %>%
           str_replace_all(fixed(";NA"), ""))

saveRDS(seqtab, file = target)
