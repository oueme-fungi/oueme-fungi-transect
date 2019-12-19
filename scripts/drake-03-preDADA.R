if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  od <- readLines(snakemake@input$flag)
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
ncpus <- local_cpus()
predada_targets <- c(
  purrr::keep(od, startsWith, "derep2_"),
  purrr::keep(od, startsWith, "qstats")
)
if (length(predada_targets)) {
  flog.info("Making pre-dada targets (%s with %d jobs)", options("clustermq.scheduler"), ncpus)
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
       memory_strategy = "preclean",
       console_log_file = get_logfile("predada")
  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
  flog.info("Recalculating outdated targets...")
  od <- drake::outdated(drake::drake_config(plan, jobs_preprocess = local_cpus()))
  flog.info("Finished.")
} else flog.info("Pre-DADA2 targets are up-to-date.")

if (exists("snakemake")) {
  writeLines(od, snakemake@output$flag)
}
