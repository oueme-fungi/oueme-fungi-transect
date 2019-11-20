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
options(clustermq.scheduler = "multicore")

preitsx_cpus <- 1
preitsx_jobs <- max(local_cpus() %/% preitsx_cpus, 1)
preitsx_cpus <- max(local_cpus() %/% preitsx_jobs, 1)

#### Pre-ITSx targets ####
# These are computationally easy, but some take a lot of memory.
# run them locally with 2 cores (worth of memory) each
preitsx_targets <- c(stringr::str_subset(od, "^split_derep_"),
                     stringr::str_subset(od, "^derep_submap_"))
if (length(preitsx_targets)) {
  flog.info("\nMaking targets to prepare for ITSx...")
  tictoc::tic()
  dconfig <- drake::drake_config(
    plan,
    parallelism = "clustermq",
    jobs_preprocess = local_cpus(),
    jobs = preitsx_jobs,
    retries = 2,
    elapsed = 3600, # 1 hour
    keep_going = FALSE,
    caching = "worker",
    cache_log_file = TRUE,
    targets = preitsx_targets,
    memory_strategy = "preclean"
  )
  
  dod <- drake::outdated(dconfig)
  
  drake::make(config = dconfig)
  
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else flog.info("\n All pre-itsx targets are up-to-date.")
