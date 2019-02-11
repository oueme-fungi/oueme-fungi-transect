library(magrittr)
library(tidyverse)
library(Biostrings)
library(here)

if (interactive()) {
  data.dir <- here("data")
  lab.dir <- here("config")
  seq.dir <- here("sequences")
  dataset.file <- file.path(lab.dir, "datasets.csv")
} else {
  dataset.file <- Sys.getenv("DATASET")
  target <- Sys.getenv("TARGETLIST")
  
  # get the prereq files on stdin
  # there are too many to pass as an environmental variable!
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  prereqs <- str_split(prereqs, " ") %>% unlist
  cat("prereqs: ", prereqs, "\n")
  
  in.files <- str_subset(prereqs,
                         pattern = "\\.trim\\.fastq\\.gz$")
  cat("in.files: ", in.files, "\n")
}