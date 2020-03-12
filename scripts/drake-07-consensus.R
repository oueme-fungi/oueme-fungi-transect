if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  od <- readLines(snakemake@input$flag)
} else {
  load("data/plan/drake.Rdata")

}

library(magrittr)
library(backports)
library(futile.logger)
setup_log("consensus")
options(clustermq.scheduler = "multicore")

#### Taxonomy targets from DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)
targets <- c("taxon_table",
             plan$target %>% purrr::keep(startsWith, "aln_decipher_"))

targets <- intersect(targets, od)

ncpus <- max(local_cpus() %/% 2L, 1L)
jobs <- max(1, local_cpus() %/% ncpus)

if (length(targets) > 0) {
  jobs <- min(max(local_cpus() %/% 8, 1), length(targets))
  ncpus <- max(local_cpus() %/% jobs, 1)
  parallelism <- if (jobs > 1) "clustermq" else "loop"
  flog.info("Building %d consensus and taxonomy targets with %d jobs of %d cores...",
            length(targets), jobs, ncpus)
  flog.info("Making targets %s...", paste(targets, collapse = "; "))
  cache <- drake::drake_cache(cache_dir)
  tictoc::tic()
  drake::make(
    plan,
    parallelism = parallelism,
    jobs_preprocess = local_cpus(),
    jobs = jobs,
    retries = 1,
    elapsed = 3600*48, #48 hours
    keep_going = FALSE,
    caching = "worker",
    cache_log_file = TRUE,
    targets = targets,
    memory_strategy = "preclean",
    console_log_file = get_logfile("consensus"),
    cache = cache
  )
  tictoc::toc()
  if (any(targets %in% drake::failed())) {
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
} else flog.info("Consensus targets are up-to-date.")

if (exists("snakemake")) {
  writeLines(od, snakemake@output$flag)
}
