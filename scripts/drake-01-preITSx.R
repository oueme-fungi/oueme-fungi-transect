if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

library(magrittr)
library(backports)
library(futile.logger)
setup_log("preITSx")

#### Pre-ITSx targets ####
# These are computationally easy, but some take a lot of memory, and would be
# inefficient to send to SLURM workers.  It's better to just do them locally.
preitsx_targets <- c(stringr::str_subset(od, "^split_fasta_"),
                     stringr::str_subset(od, "^join_derep_map"))
if (length(preitsx_targets)) {
  flog.info("\nMaking targets to prepare for ITSx...")
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
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else flog.info("\n All pre-itsx targets are up-to-date.")
