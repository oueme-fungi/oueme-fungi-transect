if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")

}

library(magrittr)
library(backports)
library(futile.logger)
setup_log("consensus")
options(clustermq.scheduler = "multicore")

#### Taxonomy targets from DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)
targets <- c("cmaln_long", "guidetree_32S", "taxon_table",
             plan$target %>% purrr::keep(startsWith, "aln_decipher_"))

targets <- subset_outdated(targets, dconfig)
if (!file.exists(config$cmaln_long)) targets <- c(targets, "cmaln_long")
if (!file.exists(config$guide_tree)) targets <- c(targets, "guidetree_32S")

ncpus <- max(local_cpus() %/% 2L, 1L)
jobs <- max(1, local_cpus() %/% ncpus)

if (length(targets) > 0) {
  jobs <- min(max(local_cpus() %/% 8, 1), length(targets))
  ncpus <- max(local_cpus() %/% jobs, 1)
  parallelism <- if(jobs > 1) "clustermq" else "loop"
  flog.info("Building %d consensus and taxonomy targets with %d jobs of %d cores...",
            length(targets), jobs, ncpus)
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = parallelism,
       jobs_preprocess = local_cpus(),
       jobs = jobs,
       retries = 1,
       elapsed = 3600*48, #48 hours
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = targets,
       console_log_file = get_logfile("consensus")
  )
  tictoc::toc()
  flog.info("Finding outdated targets...")
  tictoc::tic()
  od <- subset_outdated(targets, dconfig)
  tictoc::toc()
  flog.info("Found %d outdated targets.", length(od))
  flog.info("Making targets %s...", paste(od, collapse = "; "))
  tictoc::tic()
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(od %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
  if (!file.exists(config$cmaln_long)) {
    flog.info("Creating %s.", config$cmaln_long)
    tictoc::tic()
    drake_build(cmaln_long, dconfig)
    tictoc::toc()
  }
  if (!file.exists(config$guide_tree)) {
    flog.info("Creating	%s.", config$guide_tree)
    tictoc::tic()
    drake_build(guidetree_32S, dconfig)
    tictoc::toc()
  }

} else {
  flog.info("Consensus targets are up-to-date.")
}

Sys.setFileTime(config$cmaln_long, Sys.time())
Sys.setFileTime(config$guide_tree, Sys.time())
