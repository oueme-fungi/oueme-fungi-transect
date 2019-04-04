if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

library(magrittr)
library(backports)
setup_log("finish")

#### Finish ####
# For now the later steps are not very intensive, so they can be done
# using the resources of the master computer.
if (length(od)) {
  cat("\n Making all remaining targets (loop)...\n")
  tictoc::tic()
  future::plan(strategy = "multiprocess")
  drake::make(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus(),
       retries = 2,
       elapsed = 600, # 10 minutes
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE
  )
  tictoc::toc()
  if (length(drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
}
cat("\nAll targets are up-to-date.\n")