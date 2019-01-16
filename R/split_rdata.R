library(magrittr)
library(tidyverse)

if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  
  # choose a dataset and run for testing.
  dataset <- "short-pacbio"
  seq.run <- "pb_483"
  in.files <- list.files(path = data.dir,
                        pattern = "\\.Rdata$",
                        full.names = TRUE) 
} else {
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  prereqs <- str_split(prereqs, " ") %>% unlist
  
  cat("prereqs: ", prereqs, "\n")
  in.files <- str_subset(prereqs,
                         pattern = "\\.Rdata$")
  cat("in.files: ", in.files, "\n")
}

for (f in in.files) {
  e <- new.env()
  load(f, envir = e)
  for (o in names(e)) {
    saveRDS(get(o, e), file = str_replace(f, "\\.Rdata$", paste0(".", o, ".rds")))
  }
  remove(e)
}
