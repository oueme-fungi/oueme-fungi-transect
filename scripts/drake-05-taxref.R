if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

target <- get_target(default = "refdb_unite_ITS")
target <- sub("refdb", "", target)
targets <- plan$target[endsWith(plan$target, target)]

library(magrittr)
library(backports)
library(futile.logger)
setup_log("taxonomy-refs")

#### pre-taxonomy ####
# single-threaded targets after dada2, before taxonomy.
if (any(targets %in% od)) {
  flog.info("Making pre-taxonomy targets for %s (loop)...", target)
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus(),
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       cache_log_file = TRUE,
       targets = pretaxon_targets,
       console_log_file = get_logfile("taonomy-refs")
  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else flog.info("Pre-taxonomy targets are up-to-date.")
