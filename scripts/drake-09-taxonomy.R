if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

target <- get_target(default = "taxon_ITS2_unite_ITS2")
target <- plan$target[startsWith(plan$target, target)]

library(magrittr)
library(backports)
library(futile.logger)
setup_log("taxonomy")

#### Taxonomy targets from DADA2 pipeline ####
# All of the taxonomy assignment algorithms are internally parallel,
# so these need to be sent to nodes with multiple cores (and incidentally a lot
# of memory)

dada_cpus <- local_cpus()
if (any(target %in% od)) {
  flog.info("Making %s with %d cores...", target, dada_cpus)
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = dada_cpus,
       retries = 1,
       elapsed = 3600*6, #6 hours
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = target
  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else flog.info("Target is up-to-date.", target)