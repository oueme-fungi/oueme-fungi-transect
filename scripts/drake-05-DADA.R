if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  od <- readLines(snakemake@input$flag)
} else {
  load("drake.Rdata")
}

targets <- c(
  purrr::keep(od, startsWith, "chimeras_"),
  purrr::keep(od, startsWith, "lsux_illumina_")
)
library(magrittr)
library(backports)
library(futile.logger)
options(clustermq.scheduler = "multicore")
setup_log("dada")

#### DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)
jobs <- 1
ncpus <- local_cpus()

if (length(targets) > 0) {
  flog.info("Making %d dada targets with %d jobs of %d cores...", length(targets), jobs, ncpus)
  tictoc::tic()
  cache <- drake::drake_cache(cache_dir)
  dconfig <- drake::drake_config(plan,
       parallelism = if (jobs > 1) "clustermq" else "loop",
       jobs_preprocess = local_cpus(),
       jobs = jobs,
       retries = 1,
       elapsed = 3600*6, #6 hours
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = targets,
       console_log_file = get_logfile("dada"),
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
} else flog.info("DADA targets are up-to-date.")

if (exists("snakemake")) {
  writeLines(od, snakemake@output$flag)
}
