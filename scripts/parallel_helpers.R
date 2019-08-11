#### parallel setup ####

# are we running slurm?
is_slurm <- function() nchar(Sys.which("sbatch")) > 0
is_local <- function() !is_slurm()

# how many cpus do we have on the local machine?
# if we're not running on the cluster, leave one cpu free.
local_cpus <- function() {
  if (exists("snakemake")) {
    snakemake@threads
  } else if (is_slurm()) {
    out <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
    assertthat::assert_that(assertthat::is.count(out))
    out
  } else {
    max(parallel::detectCores() - 1, 1)
  }
}

# how many cpus can we get on one node?
max_cpus <- function() {
  if (is_slurm()) {
    slurminfo <- system("sinfo -O cpus", intern = TRUE)[-1]
    slurminfo <- as.integer(slurminfo)
    return(max(slurminfo))
  } else {
    return(local_cpus())
  }
}

#### If we were passed a target, what is it?
get_target <- function(default) {
  if (interactive()) {
    return(default)
  } else if (exists("snakemake")) {
    return(sub("^.", "", snakemake@output[1]))
  } else {
    return(sub("^.", "", Sys.getenv("TARGETLIST")))
  }
}

#### what is our logfile called?
get_logfile <- function(default) {
  if (interactive()) {
    return(file.path("logs", paste0(default, ".log")))
  } else if (exists("snakemake")) {
    return(file(snakemake@log[[1]], open = "at"))
  } else {
    return(file(file.path(Sys.getenv("LOGDIR"), paste0(default, ".log"))))
  }
}

#### Start logging ####
setup_log <- function(default) {
  logfile <- get_logfile(default)
  sink(logfile, split = TRUE)
  sink(logfile, type = "message")
}

#### combine drake targets with their names ####
drake_combine <- function(...) UseMethod("drake_combine")  

drake_combine.default <- function(...) {
   list <- list(...)
   names(list) <- as.character(match.call(expand.dots = FALSE)$...)
   list
}

#### return the argument's name as a string ####
name_to_string <- function(x) deparse(substitute(x))

#### Get the name of the current Conda environment ####
get_conda_env <- function() {
   Sys.getenv("CONDA_PREFIX")
}
