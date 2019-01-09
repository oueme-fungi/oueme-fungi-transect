
# source("install_packages.R", echo = TRUE)
library(magrittr)
library(tidyverse)
library(dada2)

if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- file.path(data.dir, "lab_setup")
  seq.dir <- file.path(base.dir, "raw_data")
  
  # choose a dataset and run for testing.
  dataset <- "short-pacbio"
  seq.run <- "pb_483"
  fullstem <- file.path(seq.dir, dataset, seq.run)
  trim.dir <- file.path(fullstem, "trim")
  in.files <- list.files(path = trim.dir,
                        pattern = "\\.trim\\.fastq\\.gz$",
                        full.names = TRUE)
  
  target <- paste0(file.path(data.dir, dataset), "_", seq.run, ".RDS")
} else {
  cat("prereqs: ", prereqs, "\n")
  in.files <- list.files(path = prereqs,
                         pattern = "\\.trim\\.fastq\\.gz$",
                         full.names = TRUE)
  cat("in.files: ", in.files, "\n")
}

in.files %<>% keep(~ file.size(.) > 40)

# read the dataset definitions
datasets <- read_csv(dataset.file)
# figure out which dataset this file belongs to
dataset <- str_extract(in.files[1], datasets$Dataset) %>% na.omit()
cat("dataset: ", dataset, "\n")
seq.run <- str_extract(in.files[1], datasets$Seq.Run) %>% na.omit()
cat("seq.run: ", seq.run, "\n")
# choose only the relevant row.
datasets %<>% filter(Dataset == dataset,
                     Seq.Run == seq.run)
err.fun <- ifelse(datasets$Tech == "PacBio",
                  PacBioErrfun,
                  loessErrfun)

# from here, based on sample script at http://benjjneb.github.io/dada2/bigdata.html
# File parsing
sample.names <- sapply(strsplit(basename(in.files), ".", fixed = TRUE), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
names(in.files) <- sample.names
# Learn error rates
set.seed(100)
err <- learnErrors(in.files, nbases = 1e8, multithread=TRUE, randomize=TRUE,
                   errorEstimationFunction = err.fun,
                   verbose = TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
#for(sam in sample.names) {
#  cat("Processing:", sam, "\n")
  derep <- derepFastq(in.files)#[[sam]])
  dds <- dada(derep, err=err, multithread=TRUE, pool = TRUE)
#}
# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)
saveRDS(dds, target) # CHANGE ME to where you want sequence table saved