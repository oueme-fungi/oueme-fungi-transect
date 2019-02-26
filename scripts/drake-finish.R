if (interactive()) {
  r.dir <- here::here("scripts")
} else if (exists("snakemake")) {
  r.dir <- snakemake@config$rdir
} else {
  r.dir <- Sys.getenv("RDIR")
}

source(file.path(r.dir, "drake.R"))



#### Finish ####
# For now the later steps are not very intensive, so they can be done
# using the resources of the master computer.
if (length(od)) {
  cat("\n Making all remaining targets (loop)...\n")
  tictoc::tic()
  future::plan(strategy = "multiprocess")
  make(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus,
       retries = 2,
       elapsed = 600, # 10 minutes
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE
  )
  tictoc::toc()
  if (length(failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
}
cat("\nAll targets are up-to-date.\n")