if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  od <- readLines(snakemake@input$flag)
} else {
  load("drake.Rdata")
}


library(magrittr)
library(backports)
library(drake)
library(futile.logger)
setup_log("raxml")


targets <- purrr::keep(od, startsWith, "raxml_")

#### RaxML targets ####
# Raxml is parallel,
# so these need to be sent to nodes with multiple cores (and incidentally a lot
# of memory)

ncpus <- local_cpus()
if (length(targets)) {
  flog.info(
    "Building RAxML targets %s with %d cores...",
    paste(targets, collapse = ","),
    ncpus
  )
  tictoc::tic()
  cache <- drake::drake_cache(cache_dir)
  drake::make(
    plan,
    parallelism = "loop",
    jobs_preprocess = local_cpus(),
    jobs = 1,
    retries = 0,
    elapsed = 3600*72, #3 days
    keep_going = TRUE,
    caching = "worker",
    cache_log_file = TRUE,
    targets = targets,
    console_log_file = get_logfile("raxml"),
    prework = "options(expressions = 10000)",
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
} else flog.info("RAxML targets are up-to-date.")

if (exists("snakemake")) {
  writeLines(od, snakemake@output$flag)
}
