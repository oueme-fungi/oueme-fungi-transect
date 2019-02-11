library(magrittr)
library(tidyverse)
library(assertthat)
library(ShortRead)

# Works basically like a stripped-down version of dada2::filterAndTrim,
# but also has the option to filter on error rate.
quality_filter <- function(in.file, out.file,
                           MaxERate = Inf,
                           MaxEE = Inf,
                           MinLength = 0,
                           MaxLength = Inf) {
  
  # initialize the output file with an empty fastq.
  blankread <- ShortReadQ()
  if (!dir.exists(dirname(out.file))) {
    dir.create(dirname(out.file), recursive = TRUE)}
  if (file.exists(out.file)) file.remove(out.file)
  writeFastq(blankread, out.file)
  
  # stream the input file to keep memory usage manageable
  fqs <- FastqStreamer(source)
  while (length(fq <- yield(fqs)) > 0) {
    #calculate number of expected errors for each read
    eexp <- 10^(-1 * as(fq@quality, "matrix")/10) %>%
      rowSums(na.rm = TRUE)
    # get length of each read
    l <- width(fq)
    # calculate error rate
    erate <- eexp / l
    # filter out files which have too many errors or are too short or too long,
    # according to the dataset definition file
    fq <- fq[erate <= MaxERate &
               eexp <= MaxEE &
               l >= MinLength &
               l <= MaxLength]
    # write to the new file
    writeFastq(object = fq, 
               file = out.file,
               mode = "a")
  }
}

# for scripted use, these are specified in the makefile
if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- here("config")
  seq.dir <- file.path(base.dir, "raw_data")
  
  # choose a dataset and run for testing.
  dataset <- "short-ion"
  seq.run <- "is_057"
  fullstem <- file.path(seq.dir, dataset, seq.run)
  demux.dir <- file.path(fullstem, "demultiplex")
  
  # for interactive testing; choose the file in the dataset with the most reads.
  in.file <- list.files(path = demux.dir, pattern = "*.fastq.gz", full.names = TRUE)
  in.file <- in.file[!grepl(pattern = ".trim.", in.file, fixed = TRUE)]
  in.file <- in.file[1]
  
  out.file <- sub(pattern = ".fastq.gz", replacement = ".trim.fastq.gz",
                  x = in.file, fixed = TRUE)
  out.file <- sub(pattern = "demultiplex", replacement = "trim", x = out.file,
                  fixed = TRUE)
  in.file <- "/home/brendan/Documents/Uppsala/Projects/oueme-fungi-transect/raw_data/short-pacbio/pb_483/demultiplex/short-pacbio-002_H1.fastq.gz"
  out.file <- "/home/brendan/Documents/Uppsala/Projects/oueme-fungi-transect/raw_data/short-pacbio/pb_483/trim/short-pacbio-002_H1.trim.fastq.gz"
  dataset.file <- file.path(lab.dir, "datasets.csv")
} else {
  out.file <- Sys.getenv("TARGETLIST")
  dataset.file <- Sys.getenv("DATASET")
  assert_that(is.string(out.file))
  
  # read the prereqs from stdin
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  
  prereqs <- str_split(prereqs, " ") %>% unlist
  
  cat("prereqs: ", prereqs, "\n")
  in.file <- prereqs[grep(pattern = ".fastq.gz", 
                  prereqs)]
  cat("in.file: ", in.file, "\n")
}

# read the dataset definitions
datasets <- read_csv(dataset.file)
# figure out which dataset this file belongs to
dataset <- str_extract(in.file, datasets$Dataset) %>% na.omit()
cat("dataset: ", dataset, "\n")
# choose only the relevant row.
datasets %<>% filter(Dataset == dataset)

stopifnot(nrow(datasets) == 1,
          file.exists(in.file))

quality_filter(in.file = in.file, out.file = out.file,
               MaxERate = dataset$MaxERate,
               MaxEE = dataset$MaxEE,
               MaxLength = dataset$MaxLength,
               MinLength = dataset$MinLength)