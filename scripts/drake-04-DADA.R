if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

targets <- purrr::keep(od, startsWith, "nochim")
targets <- subset_outdated(targets, dconfig)
library(magrittr)
library(backports)
library(futile.logger)
library(clustermq)
options(clustermq.scheduler = "multicore")
setup_log("dada")

#### DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)
jobs <- max(local_cpus() %/% 8, targets)
dada_cpus <- local_cpus() %/% jobs

if (length(targets) > 0) {
  flog.info("Making %d dada targets with %d jobs of %d cores...", length(targets), jobs, dada_cpus)
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "clustermq",
       jobs_preprocess = local_cpus(),
       jobs = jobs,
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
} else flog.info("DADA targets are up-to-date.")
