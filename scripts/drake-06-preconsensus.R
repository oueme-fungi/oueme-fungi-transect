if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  od <- readLines(snakemake@input$flag)
  outputs <- unique(unlist(snakemake@output))
} else {
  load("drake.Rdata")
  outputs <- "data/clusters/ITS2.fasta.gz"
}

targets <- purrr::keep(od, startsWith, "big_fasta")

library(magrittr)
library(backports)
library(futile.logger)
setup_log("pretaxonomy")
options(clustermq.scheduler = "multicore")
ncpus <- 1

#### pre-taxonomy ####
# single-threaded targets after dada2, before consensus.
if (length(targets) > 0) {
  flog.info("Making pre-consensus targets with %d jobs of one core each...", local_cpus())
  tictoc::tic()
  cache <- drake::drake_cache(cache_dir)
  dconfig <- drake::drake_config(plan,
       parallelism = "clustermq",
       jobs_preprocess = local_cpus(),
       jobs = local_cpus(),
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       cache_log_file = TRUE,
       targets = targets,
       console_log_file = get_logfile("pretaxonomy"),
       memory_strategy = "preclean",
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
} else flog.info("Pre-consensus targets are up-to-date.")

if (exists("snakemake")) {
  writeLines(od, snakemake@output$flag)
}
