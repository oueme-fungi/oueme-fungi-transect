if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  aln_file_long <- snakemake@output$cmaln_long
  guide_tree_file <- snakemake@output$guide_tree
} else {
  load("drake.Rdata")
  cmaln_file_long <- file.path(locarna_dir, "long_cmalign.aln")
  guide_tree_file <- file.path(locarna_dir, "32S_guide.tree")

}

library(magrittr)
library(backports)
library(futile.logger)
setup_log("consensus")

#### Taxonomy targets from DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)
targets <- c("cmaln_long", "guidetree_32S")

dada_cpus <- local_cpus()

if (any(targets %in% od) || !all(file.exists(aln_file_LSU, aln_file_32S))) {
  flog.info("Configuring drake  %s with %d core(s)...", targets, dada_cpus)
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = dada_cpus,
       retries = 1,
       elapsed = 3600*6, #6 hours
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = targets
  )
  tictoc::toc()
  flog.info("Finding outdated targets...")
  tictoc::tic()
  od <- subset_outdated(targets, dconfig)
  tictoc::toc()
  flog.info("Found %d outdated targets.", length(od))
  flog.info("Making targets...")
  tictoc::tic()
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(od %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
  if (!file.exists(cmaln_file_long)) {
    flog.info("Creating %s.", cmaln_file_long)
    tictoc::tic()
    drake_build(cmaln_long, dconfig)
    tictoc::toc()
  }
  if (!file.exists(guide_tree_file)) {
    flog.info("Creating	%s.", guide_tree_file)
    tictoc::tic()
    drake_build(guidetree_32S, dconfig)
    tictoc::toc()
  }


} else {
  flog.info("Long ASV consensus sequences are up-to-date.")
}

Sys.setFileTime(cmaln_file_long, Sys.time())
Sys.setFileTime(guide_tree_file, Sys.time())
