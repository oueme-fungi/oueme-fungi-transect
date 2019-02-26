if (interactive()) {
  r.dir <- here::here("scripts")
} else if (exists("snakemake")) {
  r.dir <- snakemake@config$rdir
} else {
  r.dir <- Sys.getenv("RDIR")
}

source(file.path(r.dir, "drake.R"))

#### Pre-ITSx targets ####
# These are computationally easy, but some take a lot of memory, and would be
# inefficient to send to SLURM workers.  It's better to just do them locally.
preitsx_targets <- str_subset(od, "^split_fasta_")
if (length(preitsx_targets)) {
  cat("\nMaking targets to prepare for ITSx...\n")
  tictoc::tic()
  make(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus,
       retries = 2,
       elapsed = 3600, # 1 hour
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = preitsx_targets
  )
  tictoc::toc()
  if (length(failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n All pre-itsx targets are up-to-date.\n")
