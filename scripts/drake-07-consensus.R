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

#### Taxonomy targets from DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)
targets <- c("cmaln_long", "guidetree_32S", "realign_long", "big_reconstructed_pb_500")

dada_cpus <- local_cpus()

if (any(targets %in% od) || !all(file.exists(cmaln_file_long, guide_tree_file, mlocarna_result_file))) {
  flog.info("Configuring drake to build consensus targets with %d core(s)...", dada_cpus)
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = dada_cpus,
       retries = 1,
       elapsed = 3600*24, #24 hours
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
  flog.info("Making targets %s...", paste(od, collapse = "; "))
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
  if (!file.exists(mlocarna_aln_file)) {
    flog.info("Creating	%s.", mlocarna_alnfile)
    tictoc::tic()
    drake_build(realign_long, dconfig)
    tictoc::toc()
  }


} else {
  flog.info("Consensus targets are up-to-date.")
}

Sys.setFileTime(cmaln_file_long, Sys.time())
Sys.setFileTime(guide_tree_file, Sys.time())
