if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

library(magrittr)
library(backports)
setup_log("preITSx")

#### Pre-ITSx targets ####
# These are computationally easy, but some take a lot of memory, and would be
# inefficient to send to SLURM workers.  It's better to just do them locally.
preitsx_targets <- stringr::str_subset(od, "^split_fasta_")
if (length(preitsx_targets)) {
  cat("\nMaking targets to prepare for ITSx...\n")
  tictoc::tic()
  dconfig <- drake::drake_config(
    plan,
    parallelism = "loop",
    jobs_preprocess = local_cpus(),
    retries = 2,
    elapsed = 3600, # 1 hour
    keep_going = FALSE,
    caching = "worker",
    cache_log_file = TRUE,
    targets = preitsx_targets
  )
  
  dod <- drake::outdated(dconfig)
  
  drake::make(config = dconfig)
  
  tictoc::toc()
  if (any(dod %in% drake::failed()))) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n All pre-itsx targets are up-to-date.\n")
