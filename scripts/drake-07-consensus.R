if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  aln_file_LSU <- snakemake@output$aln_LSU
  aln_file_32S <- snakemake@output$aln_32S
} else {
  load("drake.Rdata")
  aln_file_LSU <- file.path(pasta.dir, "lsu_ASVs.fasta")
  aln_file_32S <- file.path(pasta.dir, "32S_ASVs.aln")

}

library(magrittr)
library(backports)
library(futile.logger)
setup_log("consensus")

#### Taxonomy targets from DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)
targets <- c("allseqs", "write_aln_LSU", "write_aln_32S")

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
  if (!file.exists(aln_file_LSU)) {
    flog.info("Creating %s.", aln_file_LSU)
    tictoc::tic()
    drake_build(write_aln_LSU, dconfig)
    tictoc::toc()
  }
  if (!file.exists(aln_file_32S)) {
    flog.info("Creating	%s.", aln_file_32S)
    tictoc::tic()
    drake_build(write_aln_32S, dconfig)
    tictoc::toc()
  }


} else {
  flog.info("Long ASV consensus sequences are up-to-date.")
}

Sys.setFileTime(aln_file_LSU, Sys.time())
Sys.setFileTime(aln_file_32S, Sys.time())
