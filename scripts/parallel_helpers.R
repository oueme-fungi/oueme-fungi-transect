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
  flog.appender(appender.tee(logfile))
}

#### combine drake targets with their names ####
drake_combine <- function(...) UseMethod("drake_combine")  

drake_combine.default <- function(...) {
   list <- list(...)
   names(list) <- as.character(match.call(expand.dots = FALSE)$...)
   list
}

#### return the argument's name(s) as (a) string(s) ####
name_to_string <- function(x) deparse(substitute(x))
names_to_strings <- function(...) {
  s <- rlang::ensyms(...)
  purrr::map_chr(s, rlang::as_string)
}

#### undo conversion of values to symbols ####
# transform = combine() in drake will convert columns which are not intended to
# be symbols into symbols, e.g. the value 3L is converted to `3L` and "yes" is
# converted to `"yes"`.  This function converts them back.
symbols_to_values <- function(...) {
  s <- rlang::ensyms(...)
  s <- purrr::map_chr(s, rlang::as_string)
  s <- purrr::map(s, ~parse(text = .))
  purrr::map_chr(s, eval)
}

#### Get the name of the current Conda environment ####
get_conda_env <- function() {
   Sys.getenv("CONDA_PREFIX")
}


#### Try several times with an exponential backoff and logging.
exp_try <- function(x, t, max) {
  if (t >= max) return(x)
  x <- rlang::enquo(x)
  futile.logger::ftry(rlang::eval_tidy(x),
                      error = function(e) {
                        Sys.sleep(t)
                        exp_try(rlang::eval_tidy(x), t*2, max)
                      })
}


is_outdated <- function(task, dconfig) {
  assertthat::assert_that(assertthat::is.string(task))
  tryCatch(
    any(drake::deps_profile(task, dconfig, character_only = TRUE)$changed),
    error = function(e) TRUE
  )
}

which_outdated <- function(tasks, dconfig) {
  assertthat::assert_that(is.character(tasks))
  which(purrr::map_lgl(tasks, is_outdated, dconfig = dconfig))
}

#### Count the number of tasks which still need to be done from a list (without calculating all outdated tasks)
n_outdated <- function(tasks, dconfig) {
  length(which_outdated(tasks, dconfig))
}

subset_outdated <- function(tasks, dconfig) {
  tasks[which_outdated(tasks, dconfig)]
}
