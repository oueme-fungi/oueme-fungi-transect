if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

library(magrittr)
library(backports)
setup_log("pretaxonomy")

#### pre-DADA2 ####
# single-threaded targets after itsx
# for local runs, ITSx targets will also run here.
pretaxon_targets <- stringr::str_subset(od, "^big_seq_table_")
if (length(pretaxon_targets)) {
  cat("\n Making pre-taxonomy targets (loop)...\n")
  tictoc::tic()
  drake::make(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus(),
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       cache_log_file = TRUE,
       targets = pretaxon_targets
  )
  tictoc::toc()
  if (length(drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n Pre-taxonomy targets are up-to-date. \n")