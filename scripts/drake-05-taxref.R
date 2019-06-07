if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

targets <- c("unite_db1", "rdp_db1", "silva_db1")
library(magrittr)
library(backports)
setup_log("taxonomy-refs")

#### pre-taxonomy ####
# single-threaded targets after dada2, before taxonomy.
if (any(targets %in% od)) {
  cat("\n Making pre-taxonomy targets (loop)...\n")
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus(),
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       cache_log_file = TRUE,
       targets = pretaxon_targets
  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n Pre-taxonomy targets are up-to-date. \n")
