if (interactive()) {
  r.dir <- here::here("scripts")
} else if (exists("snakemake")) {
  r.dir <- snakemake@config$rdir
} else {
  r.dir <- Sys.getenv("RDIR")
}

source(file.path(r.dir, "drake.R"))

#### pre-DADA2 ####
# single-threaded targets after itsx
# for local runs, ITSx targets will also run here.
predada_targets <- c(str_subset(od, "^derep2_"))
if (length(predada_targets)) {
  cat("\n Making pre-dada targets (loop)...\n")
  tictoc::tic()
  make(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus,
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = predada_targets
  )
  tictoc::toc()
  if (length(failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n Pre-DADA2 targets are up-to-date. \n")