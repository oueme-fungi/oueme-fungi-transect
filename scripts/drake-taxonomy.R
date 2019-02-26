if (interactive()) {
  r.dir <- here::here("scripts")
} else if (exists("snakemake")) {
  r.dir <- snakemake@config$rdir
} else {
  r.dir <- Sys.getenv("RDIR")
}

source(file.path(r.dir, "drake.R"))



#### Taxonomy targets from DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory) or, if local, use all the cores of
# the local machine

taxon_target <- target
dada_cores <- local_cores
if (target %in% od) {
  cat("\n Making taxonomy targets (local with", dada_cores, "cores)...\n")
  tictoc::tic()
  make(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus,
       retries = 1,
       elapsed = 3600*6, #6 hours
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = dada_target
  )
  tictoc::toc()
  if (length(failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n Taxonomy targets are up-to-date.\n")