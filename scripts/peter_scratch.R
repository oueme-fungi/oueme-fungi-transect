library(drake)
library(magrittr)
library(tidyverse)

loadd(conseq)

rs2_LSU <- unique(conseq$LSU)

peter_LSU <- Biostrings::readDNAStringSet("../oueme-fungi-soils/data/peter/LSU.fasta")

new_LSU <- Biostrings::readDNAStringSet("../oueme-fungi-soils/data/peter/newLSU.fasta")

rs2_LSU <- gsub("U", "T", rs2_LSU)
rs2_LSU <- setdiff(rs2_LSU, as.character(peter_LSU))
rs2_LSU <- setdiff(rs2_LSU, as.character(new_LSU))
rs2_LSU <- rs2_LSU[!is.na(rs2_LSU)]

dplyr::last(names(new_LSU))

names(rs2_LSU) <- paste0("ASV_", seq_along(rs2_LSU) + 1678)
rs2_LSU <- Biostrings::DNAStringSet(rs2_LSU)
Biostrings::writeXStringSet(rs2_LSU, "../oueme-fungi-soils/data/peter/rs2LSU.fasta")
