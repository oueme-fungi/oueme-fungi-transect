if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
  od <- readLines(snakemake@input$flag)
} else {
  load("data/plan/drake.Rdata")
}

library(magrittr)
library(backports)
library(futile.logger)
setup_log("finish")

#### Finish ####
# For now the later steps are not very intensive, so they can be done
# using the resources of the master computer.
targets <- od

jobs <- local_cpus()
ncpus <- 1
options(clustermq.scheduler = "multicore")

if (jobs > 1 && length(targets) > 1) parallelism <- "clustermq" else parallelism <- "loop"

cache <- drake::drake_cache(cache_dir)
if (length(targets)) {
  flog.info(
    "Making all remaining targets: [%s] (%s with %d job(s) of %d core(s))...",
    paste(targets, collapse = ","),
    parallelism,
    jobs,
    ncpus
  )
  tictoc::tic()
  dconfig <- drake::drake_config(plan,
       parallelism = parallelism,
       jobs_preprocess = local_cpus(),
       jobs = jobs,
       retries = 2,
       elapsed = 600, # 10 minutes
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       console_log_file = get_logfile("finish"),
       targets = targets,
       cache = cache
  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
}
flog.info("All targets are up-to-date.")


cache2 <- drake::drake_cache(".drake")
if (is.null(cache2)) cache2 <- drake::new_cache(".drake")
exports <- intersect(
  c(
    "allseqs",
    "taxon_table",
    dplyr::filter(plan, step == "taxon")$target,
    "raxml_decipher_unconst_long",
    dplyr::filter(plan, step == "big_seq_table")$target,
    dplyr::filter(plan, step == "err")$target,
    dplyr::filter(plan, step == "chimeras")$target,
    dplyr::filter(plan, step == "allchimeras")$target,
    purrr::keep(plan$target, startsWith, "file_meta_"),
    purrr::keep(plan$target, startsWith, "illumina_group_"),
    "qstats"
  ),
  cache$list()
)
cache2$import(
  from = cache,
  list = exports
)
