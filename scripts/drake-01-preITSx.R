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

ncpus <- 1
njobs <- max(local_cpus() %/% ncpus, 1)
ncpus <- max(local_cpus() %/% njobs, 1)

#### Pre-ITSx targets ####
# These are computationally easy, but some take a lot of memory.
# run them locally with 2 cores (worth of memory) each
preitsx_targets <- c(stringr::str_subset(od, "^split_derep_"),
                     stringr::str_subset(od, "^derep_submap_"),
                     stringr::str_subset(od, "^qstats_(raw|demux)"))
if (length(preitsx_targets)) {
  flog.info("Making targets to prepare for ITSx...")
  tictoc::tic()
  dconfig <- drake::drake_config(
    plan,
    parallelism = "clustermq",
    jobs_preprocess = local_cpus(),
    jobs = njobs,
    retries = 2,
    elapsed = 3600, # 1 hour
    keep_going = FALSE,
    caching = "worker",
    cache_log_file = TRUE,
    targets = preitsx_targets,
    memory_strategy = "preclean",
    console_log_file = get_logfile("preITSx"),
    cache = cache
  )
  
  dod <- drake::outdated(dconfig)
  
  drake::make(config = dconfig)
  
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
  flog.info("Recalculating outdated targets...")
  od <- drake::outdated(
    drake::drake_config(
      plan,
      jobs_preprocess = local_cpus(),
      cache = cache
    )
  )
  flog.info("Finished.")
} else flog.info("All pre-itsx targets are up-to-date.")

if (exists("snakemake")) {
  writeLines(od, snakemake@output$flag)
}
