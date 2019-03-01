if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}

library(magrittr)
library(backports)
setup_log("ITSx")

#### ITSx ####
# itsx (actually hmmer) does have an internal parallel option, but it isn't
# very efficient at using all the cores.
# Instead, we divide the work into a large number of 
# shards and submit them all as seperate jobs on SLURM.
# failing that, do all the shards locally on the cores we have.
itsx_targets <- stringr::str_subset(od, "^itsx_shard")
# hmmer can use multiple processes per job; it tends to become I/O bound after about 4.
itsx_cpus <- 4L
itsx_cpus <- min(itsx_cpus, max_cpus())
itsx_jobs <- bigsplit

if (length(itsx_targets)) {
  # send two jobs per target to help with load balancing
  # this will tend to be one Ion Torrent job (long, first) and one PacBio job (short, second)
  itsx_jobs <- max(ceiling(length(itsx_targets) / 2), itsx_jobs)
  # in any case never send more jobs than there are targets
  itsx_jobs <- min(itsx_jobs, length(itsx_targets))
  
  # do it as SLURM jobs if possible and
  # - we don't have enough local cores to run 1 job locally
  #   -or-
  # - it would use less than twice as much total CPU time, compared to running locally
  if (is_slurm()
      && ((local_cpus() < itsx_cpus)
          || (length(itsx_targets) / itsx_jobs * (itsx_jobs * itsx_cpus + local_cpus()) <
              2 * length(itsx_targets) * itsx_cpus))) {
    itsx_parallelism <- "clustermq"
    options(clustermq.scheduler = "slurm",
            clustermq.template = "slurm_clustermq.tmpl")
    itsx_template <- list(log_file = glue("logs/itsx-%A_%a.log"),
                          ncpus = itsx_cpus,
                          memory = 7*1024,
                          timeout = 1800) # the master can take a while to send everything.
    cat("\n Making itsx_shard (SLURM with ", itsx_jobs, "worker(s))...\n")
  } else {
    itsx_cpus <- 1
    itsx_jobs <- local_cpus() %/% itsx_cpus
    itsx_cpus <- local_cpus() %/% itsx_jobs
    itsx_parallelism <- if (itsx_jobs > 1) "clustermq" else "loop"
    options(clustermq.scheduler = "multicore")
    itsx_template = list()
    cat("\n Making itsx_shard (local with ", itsx_jobs, "worker(s))...\n")
  }
  tictoc::tic()
  drake::make(plan,
       parallelism = itsx_parallelism,
       template = itsx_template,
       jobs = itsx_jobs,
       elapsed = 3600*6, #6 hours
       jobs_preprocess = local_cpus(),
       caching = "worker",
       cache_log_file = TRUE,
       targets = itsx_targets
  )
  tictoc::toc()
  if (length(drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n ITSx targets are up-to-date. \n")
