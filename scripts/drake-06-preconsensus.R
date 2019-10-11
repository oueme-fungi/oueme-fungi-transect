if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  outputs <- unique(unlist(snakemake@output))
} else {
  load("drake.Rdata")
  outputs <- "data/clusters/ITS2.fasta.gz"
}

targets <- purrr::keep(plan$target, startsWith, "big_fasta")
targets <- subset_outdated(targets, dconfig)
outputs <- setdiff(outputs, ".pretaxonomy")

library(magrittr)
library(backports)
library(futile.logger)
setup_log("pretaxonomy")
library(clustermq)
options(clustermq.scheduler = "multicore")

#### pre-taxonomy ####
# single-threaded targets after dada2, before taxonomy.
if (length(targets) > 0) {
  flog.info("Making pre-taxonomy targets with %d jobs of one core each...")
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "clustermq",
       jobs_preprocess = local_cpus(),
       jobs = local_cpus(),
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       cache_log_file = TRUE,
       targets = targets

  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else flog.info("Pre-taxonomy targets are up-to-date.")

for (f in outputs) Sys.setFileTime(f, Sys.time())
