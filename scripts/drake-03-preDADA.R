if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

library(magrittr)
library(backports)
library(futile.logger)
setup_log("predada")
options(clustermq.scheduler = "multicore")
library(disk.frame)

#### pre-DADA2 ####
# single-threaded targets after itsx
# for local runs, ITSx targets will also run here.
dada_cpus <- local_cpus()
predada_targets <- c(
  purrr::keep(od, startsWith, "filter_")
)
if (length(predada_targets)) {
  flog.info("Making pre-dada targets (loop)...")
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "clustermq",
       jobs_preprocess = local_cpus(),
       jobs = local_cpus(),
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       cache_log_file = TRUE,
       targets = predada_targets,
       caching = "worker",
       memory_strategy = "preclean"
  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else flog.info("Pre-DADA2 targets are up-to-date.")
