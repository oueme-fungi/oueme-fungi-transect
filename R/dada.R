
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
  dataset <- "short-ion"
  seq.run <- "is_057"
  fullstem <- file.path(seq.dir, dataset, seq.run)
  trim.dir <- file.path(fullstem, "trim")
  in.files <- list.files(path = trim.dir,
                        pattern = "\\.trim\\.fastq\\.gz$",
                        full.names = TRUE) 
  in.files <- in.files[1]
  
  target <- paste0(file.path(data.dir, dataset), "_", seq.run, ".dada.Rdata")
  dataset.file <- file.path(lab.dir, "datasets.csv")
} else {
  dataset.file <- Sys.getenv("DATASET")
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
derep <- derepFastq(in.files)#[[sam]])

dada.opt <- list(derep = derep, err = err, multithread = TRUE, pool = TRUE)
if (!is.na(datasets$DadaOpt)) {
  datasets$DadaOpt %>%
    str_split(",") %>%
    unlist() %>%
    trimws() %>%
    str_split("=") %>%
    walk(~ (dada.opt <<- replace(dada.opt, unlist(.)[1], eval(parse(text = unlist(.)[2])))))
}

asv <- do.call(dada, dada.opt)

object.size(asv)

map_dbl(asv, object.size)

asv.big <- map_dbl(asv, object.size) %>% which.max

asv.big <- asv[[asv.big]]

str(asv.big)

map_dbl(asv.big@.Data, object.size)

object.size(derep)

map_dbl(derep, object.size)

derep.big <- map_dbl(derep, object.size) %>% which.max
derep.big <- derep[[derep.big]]

str(derep.big)

map_dbl(derep.big, object.size)

dadamap <-
  map2(derep, asv,
       function(derep, asv) {
         m <- tibble(seq.idx = seq_along(derep$map),
                     derep.idx = derep$map,
                     derep.seq = names(derep$uniques)[derep.idx])
         m %<>%
           left_join(
             tibble(
               asv.idx = asv$map,
               derep.idx = seq_along(asv.idx),
               asv.seq = asv$sequence[asv.idx]))
         m
       })

seqtable <- makeSequenceTable(asv)

nochim <- removeBimeraDenovo(seqtable, method = "consensus",
                                      multithread = TRUE)

save(dadamap, seqtable, nochim, file = target)