if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  longASV_file <- snakemake@output$longasv
} else {
  load("drake.Rdata")
  longASV_file <- file.path(pasta.dir, "long_ASVs.fasta")
}

library(magrittr)
library(backports)
setup_log("consensus")

#### Taxonomy targets from DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)
targets <- c("lsualn", "cons_tax")

dada_cpus <- local_cpus()

if (any(targets %in% od)) {
  cat("\n Making", targets, "with", dada_cpus, "cores...\n")
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = dada_cpus,
       retries = 1,
       elapsed = 3600*6, #6 hours
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = targets
  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n Long ASV consensus sequences are up-to-date.\n")