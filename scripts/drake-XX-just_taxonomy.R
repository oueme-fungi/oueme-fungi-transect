## This script is kind of a cheat, which

if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("data/plan/drake.Rdata")
}

library(magrittr)
library(backports)
library(futile.logger)
library(drake)

ncpus <- local_cpus()
njobs <- 1
cache <- drake::drake_cache(cache_dir)
skip_allseqs <- TRUE
tax_targets <- purrr::keep(plan$target, startsWith, "taxon_")
flog.info("Making ONLY taxonomy targets, without taking earlier dependencies into account.")
make(
  plan,
  targets = tax_targets,
  parallelism = "loop",
  jobs_preprocess = local_cpus(),
  jobs = njobs,
  retries = 2,
  elapsed = 3600, # 1 hour
  keep_going = FALSE,
  caching = "worker",
  cache_log_file = TRUE,
  memory_strategy = "preclean",
  console_log_file = get_logfile("just_taxonomy"),
  cache = cache
)

cache2 <- drake::drake_cache(".drake")
if (is.null(cache2)) cache2 <- drake::new_cache(".drake")
exports <- intersect(tax_targets, cache$list())
cache2$import(
  from = cache,
  list = exports
)
