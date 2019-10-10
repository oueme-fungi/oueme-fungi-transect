if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

targets <- purrr::keep(plan$target, startsWith, "taxon_")
targets <- subset_outdated(targets, dconfig)

library(magrittr)
library(backports)
library(futile.logger)
setup_log("taxonomy")
library(clustermq)
options(clustermq.scheduler = "multicore")

#### Taxonomy targets from DADA2 pipeline ####
# All of the taxonomy assignment algorithms are internally parallel,
# so these need to be sent to nodes with multiple cores (and incidentally a lot
# of memory)

if (length(targets) > 0) {
  jobs <- min(max(local_cpus() %/% 8, 1), length(targets))
  dada_cpus <- max(local_cpus() %/% jobs, 1)
  parallelism <- if(jobs > 1) "clustermq" else "loop"

  flog.info("Building %d taxonomy targets with %d jobs of %d cores...", length(targets), jobs, dada_cpus)
  tictoc::tic()
  drake::make(
    plan,
    parallelism = parallelism,
    jobs_preprocess = local_cpus(),
    jobs = jobs,
    retries = 1,
    elapsed = 3600*6, #6 hours
    keep_going = FALSE,
    caching = "worker",
    cache_log_file = TRUE,
    targets = targets
  )
  tictoc::toc()
#  flog.info("Determining outdated targets...")
#  tictoc::tic()
#  dod <- exp_try(drake::outdated(dconfig), 30, 300)
#  tictoc::toc()
#  flog.info("Making targets...")
#  tictoc::tic()
#  drake::make(config = dconfig)
#  tictoc::toc()
  if (any(targets %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else flog.info("Target is up-to-date.", target)
