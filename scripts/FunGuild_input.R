library(magrittr)
library(tidyverse)
library(readxl)


base.dir <- str_extract(getwd(), ".+oueme-fungi-transect")
data.dir <- file.path(base.dir, "data")
pacbio.dir <- file.path(data.dir, "PacBio")
pacbio.otu.file <- file.path(pacbio.dir, "PACBIO_OTU_table.xlsx")
pacbio.funguild.input.file <- file.path(pacbio.dir, "PacBio_FUNGuild_in.txt")

iontorrent.dir <- file.path(data.dir, "IonTorrent")
iontorrent.otu.file <- file.path(iontorrent.dir, "higherLvl", "OTU.txt")
iontorrent.classification.file <- file.path(iontorrent.dir, "hiera_BLAST.txt")
iontorrent.funguild.input.file <- file.path(iontorrent.dir, "IonTorrent_FUNGuild.txt")

read_xlsx(pacbio.otu.file) %>%
  dplyr::bind_cols(
    stringr::str_split(.$X__2, stringr::fixed("|"), n = 5, simplify = TRUE) %>%
      tibble::as_data_frame() %>%
      set_names(c("Subject.Name", "Accession", "Species.Hypothesis",
                  "Ref.Type", "Classification"))) %>%
  tidyr::extract("Classification",
                 c("Kingdom", "Phylum", "Class", "Order",
                   "Family", "Genus", "Species"),
                 "k__(.+)_p__(.+)_c__(.+)_o__(.+)_f__(.+)_g__(.+)_s__(.+)") %>%
  dplyr::mutate(taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";")) %>%
  select(Group, matches("[A-H]\\d+_\\d+_\\d"), taxonomy) %>%
  write_tsv(path = pacbio.funguild.input.file)

read_tsv(iontorrent.otu.file) %>%
  left_join(read_tsv(iontorrent.classification.file)) %>%
  dplyr::mutate(taxonomy = paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep = ";")) %>%
  select(Group = OTU, matches("[A-H]\\d+"), taxonomy) %>%
  write_tsv(iontorrent.funguild.input.file)
