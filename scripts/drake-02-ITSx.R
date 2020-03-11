if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  od <- readLines(snakemake@input$flag)
} else {
  load("drake.Rdata")
}

library(magrittr)
library(backports)
library(futile.logger)
setup_log("ITSx")

#### ITSx ####
# itsx (actually hmmer) does have an internal parallel option, but it isn't
# very efficient at using all the cores.
# Instead, we divide the work into a large number of 
# shards and submit them all as seperate jobs on SLURM.
# failing that, do all the shards locally on the cores we have.

# Note: not actually using ITSx anymore
itsx_targets <- stringr::str_subset(od, "^lsux_") %>%
  purrr::discard(startsWith, "lsux_illumina")

# hmmer and cmalign can use multiple processes per job;
# it tends to become I/O bound after about 4.
ncpus <- 4L
ncpus <- min(ncpus, max_cpus())
njobs <- config$bigsplit

if (length(itsx_targets)) {
  # send two jobs per target to help with load balancing
  # this will tend to be one Ion Torrent job (long, first) and one PacBio job (short, second)
  njobs <- max(ceiling(length(itsx_targets) / 2), njobs)
  # in any case never send more jobs than there are targets
  njobs <- min(njobs, length(itsx_targets))
  
  # do it as SLURM jobs if possible and
  # - we don't have enough local cores to run 1 job locally
  #   -or-
  # - it would use less than twice as much total CPU time, compared to running locally
  if (is_slurm()
      && local_cpus() < max_cpus()
      && ((local_cpus() < ncpus)
          || (length(itsx_targets) / njobs * (njobs * ncpus + local_cpus()) <
              2 * length(itsx_targets) * ncpus))) {
    itsx_parallelism <- "clustermq"
    options(clustermq.scheduler = "slurm",
            clustermq.template = "slurm_clustermq.tmpl")
    itsx_template <- list(log_file = glue::glue("logs/itsx-%A_%a.log"),
                          ncpus = ncpus,
                          memory = 7*1024,
                          timeout = 1800) # the master can take a while to send everything.
    flog.info(" Making ITSx and LSUx shards (SLURM with %d worker(s))...", njobs)
  } else {
    njobs <- max(local_cpus() %/% ncpus, 1)
    ncpus <- max(local_cpus() %/% njobs, 1)
    itsx_parallelism <- if (njobs > 1) "clustermq" else "loop"
    options(clustermq.scheduler = "multicore")
    itsx_template = list()
    flog.info("Making ITSx and LSUx shards (local with %d worker(s) of %d cpus each)...", njobs, ncpus)
  }
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = itsx_parallelism,
       template = itsx_template,
       jobs = njobs,
       elapsed = 3600*6, #6 hours
       jobs_preprocess = local_cpus(),
       caching = "worker",
       cache_log_file = TRUE,
       targets = itsx_targets,
       console_log_file = get_logfile("ITSx"),
       cache = cache
  )
  
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
  flog.info("Recalculating outdated targets...")
  od <- drake::outdated(
    drake::drake_config(
      plan,
      jobs_preprocess = local_cpus(),
      cache = cache
    )
  )
  flog.info("Finished.")
} else flog.info("ITSx targets are up-to-date.")

if (exists("snakemake")) {
  writeLines(od, snakemake@output$flag)
}
