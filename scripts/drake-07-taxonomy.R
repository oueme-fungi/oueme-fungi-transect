if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

target <- get_target(default = "guilds_table_ITS2_unite")
library(magrittr)
library(backports)
setup_log("taxonomy")

#### Taxonomy targets from DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory)

dada_cpus <- local_cpus()
if (target %in% od) {
  cat("\n Making", target, " with", dada_cpus, "cores...\n")
  tictoc::tic()
  drake::make(plan,
       parallelism = "loop",
       jobs_preprocess = dada_cpus,
       retries = 1,
       elapsed = 3600*6, #6 hours
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = target
  )
  tictoc::toc()
  if (length(drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n Target", target, "is up-to-date.\n")