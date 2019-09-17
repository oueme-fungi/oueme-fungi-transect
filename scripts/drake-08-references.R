if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

dada_cpus <- local_cpus()

target <- get_target(default = "refdb_unite_ITS")
target <- sub("refdb", "", target)
target <- plan$target[endsWith(plan$target, target)]

library(magrittr)
library(backports)
library(futile.logger)
setup_log("pretaxonomy")

#### pre-taxonomy ####
# training the IDTAXA classifier is the big part of this;
# that is internally parallel
if (any(target %in% od) #|| !all(file.exists(outputs))
    ) {
  flog.info("Making taxonomy reference targets (loop)...")
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus(),
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       cache_log_file = TRUE,
       targets = target

  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else flog.info("Pre-taxonomy targets are up-to-date.")
