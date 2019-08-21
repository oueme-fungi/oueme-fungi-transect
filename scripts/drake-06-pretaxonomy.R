if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  outputs <- unique(unlist(snakemake@output))
} else {
  load("drake.Rdata")
  outputs <- "data/clusters/ITS2.fasta.gz"
}

target <- get_target(default = "big_fasta_ITS2")
outputs <- setdiff(outputs, paste0(".", target))

library(magrittr)
library(backports)
setup_log("pretaxonomy")

#### pre-taxonomy ####
# single-threaded targets after dada2, before taxonomy.
if (target %in% od #|| !all(file.exists(outputs))
    ) {
  cat("\n Making pre-taxonomy targets (loop)...\n")
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus(),
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       cache_log_file = TRUE,
       targets = target

  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n Pre-taxonomy targets are up-to-date. \n")
