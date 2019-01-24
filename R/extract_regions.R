library(magrittr)
library(tidyverse)
library(glue)
library(ShortRead)

if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  seq.dir <- file.path(base.dir, "raw_data")
  
  # choose a dataset and run for testing.
  dataset <- "long-pacbio"
  seq.run <- "pb_500"
  plate <- "001"
  well <- "F2"
  position.file <- glue("{file.path(seq.dir, dataset, seq.run, 'demultiplex', dataset)}-{plate}_{well}.positions.txt")
  seq.file <-  glue("{file.path(seq.dir, dataset, seq.run, 'demultiplex', dataset)}-{plate}_{well}.fastq.gz")
} else {
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  prereqs <- str_split(prereqs, " ") %>% unlist
  
  position.file <- str_subset(prereqs, fixed(".positions.txt"))
  seq.file <- str_subset(prereqs, fixed(".fastq.gz"))
}


stopifnot(file.exists(position.file),
          file.exists(seq.file))


stem <- str_replace(position.file, fixed(".positions.txt"), "")

pos <- tibble(seq = character(0),
              length = character(0),
              region = character(0),
              start = integer(0),
              end = integer(0))
if (file.size(position.file) > 0) {
  pos <- read_tsv(position.file, col_names = c("seq", "length", "SSU", "ITS1",
                                               "5_8S", "ITS2", "LSU", "comment")) %>%
    gather(key = "region", value = "pos", SSU:LSU) %>%
    mutate_at("pos", str_extract, pattern = "\\d+-\\d+") %>%
    tidyr::extract(pos, into = c("start", "end"), regex = "(\\d+)-(\\d+)") %>%
    mutate_at(vars(start:end), as.integer)
}

for (r in c("ITS1", "ITS2", "LSU")) {
  out.file <- glue("{stem}.{r}.fastq.gz")
  if (file.exists(out.file)) file.remove(out.file)
  writeFastq(ShortReadQ(), out.file)
  p <- filter(pos, region == r,
              !is.na(start),
              !is.na(end))
  fastq <- readFastq(seq.file)
  fastq <- fastq[match(p$seq, fastq@id)] %>%
    narrow(start = p$start, end = p$end)
  writeFastq(fastq, out.file, mode = "a")
}
