library(magrittr)
library(tidyverse)
library(glue)
library(ShortRead)

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
  plate <- "001"
  well <- "F2"
  position.file <- glue("{file.path(seq.dir, dataset, seq.run, 'demultiplex', dataset)}-{plate}_{well}.positions.txt")
  seq.file <-  glue("{file.path(seq.dir, dataset, seq.run, 'demultiplex', dataset)}-{plate}_{well}.fastq.gz")
  stem <- str_replace(position.file, fixed(".positions.txt"), "")
} else {
  position.file <- str_extract(prereqs, fixed(".positions.txt"))
  seq.file <- str_extract(prereqs, fixed(".fasta.gz"))
}


stopifnot(file.exists(position.file),
          file.exists(seq.file))

pos <- read_tsv(position.file, col_names = c("seq", "length", "SSU", "ITS1",
                                             "5_8S", "ITS2", "LSU", "comment")) %>%
  gather(key = "region", value = "pos", SSU:LSU) %>%
  mutate_at("pos", str_extract, pattern = "\\d+-\\d+") %>%
  tidyr::extract(pos, into = c("start", "end"), regex = "(\\d+)-(\\d+)") %>%
  mutate_at(vars(start:end), as.integer)

for (r in c("ITS1", "ITS2", "LSU")) {
  out.file <- glue("{stem}.{r}.fastq.gz")
  if (file.exists(out.file)) file.remove(out.file)
  fastq <- ShortReadQ()
  writeFastq(fastq, out.file)
  p <- filter(pos, region == r)
  fastq <- readFastq(seq.file)
  fastq <- fastq[match(p$seq, fastq@id)] %>%
    narrow(start = p$start, end = p$end)
  writeFastq(fastq, out.file, mode = "a")
}
