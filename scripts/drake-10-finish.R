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
setup_log("finish")

#### Finish ####
# For now the later steps are not very intensive, so they can be done
# using the resources of the master computer.
if (length(od)) {
  flog.info("Making all remaining targets (loop)...")
  tictoc::tic()
  future::plan(strategy = "multiprocess")
  dconfig <- drake::drake_config(plan,
       parallelism = "loop",
       jobs_preprocess = local_cpus(),
       retries = 2,
       elapsed = 600, # 10 minutes
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       console_log_file = get_logfile("finish")
  )
  dod <- drake::outdated(dconfig)
  drake::make(config = dconfig)
  tictoc::toc()
  if (any(dod %in% drake::failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
}
flog.info("All targets are up-to-date.")

cache1 <- drake::drake_cache()

cache2 <- drake::drake_cache(".light")
if (is.null(cache2)) cache2 <- drake::new_cache(".light")
cache2$import(
  from = cache1,
  list = c(
    "allseqs",
    "taxon_table",
    dplyr::filter(plan, step == "taxon")$target,
    "raxml_decipher_LSU",
    "raxml_decipher_long",
    "raxml_epa_full",
    "raxml_decipher_unconst_long",
    "raxml_epa_decipher_unconst_full",
    dplyr::filter(plan, step == "big_seq_table")$target,
    dplyr::filter(plan, step == "err")$target,
    dplyr::filter(plan, step == "chimeras")$target,
    dplyr::filter(plan, step == "allchimeras")$target
  )
)
