# utility functions, mostly related to drake
# author Brendan Furneaux

# make relative symlinks from all files in one directory tree to another
# directory tree
mirror_dir <- function(from, to) {
  files <- list.files(from, recursive = TRUE, include.dirs = FALSE)
  dirs <- unique(dirname(files))
  dir.create(
    file.path(to, dirs),
    recursive = TRUE,
    showWarnings = FALSE
  )
  newfiles <- files[!file.exists(file.path(to, files))]
  oldfiles <- setdiff(files, newfiles)
  wrongfiles <- oldfiles[normalizePath(file.path(from, oldfiles)) !=
                           normalizePath(file.path(to, oldfiles))]
  todofiles <- c(newfiles, wrongfiles)
  if (length(todofiles) > 0) {
    reldirs <- vapply(
      dirs,
      function(dir)
        system2(
          "realpath",
          c(
            "--relative-to",
            shQuote(file.path(to, dir)),
            shQuote(from)
          ),
          stdout = TRUE
        ),
      ""
    )[dirname(todofiles)]

    unlink(file.path(to, wrongfiles))
    file.symlink(
      file.path(reldirs, todofiles),
      file.path(to, todofiles)
    )
  }
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

#### combine FST targets into a disk.frame using symlinks ####
combine_diskframe <- function(..., dir = drake::drake_tempfile(), cache = drake::drake_cache()) {
  if (is.character(cache)) cache <- drake::drake_cache(cache)
  dir.create(dir, recursive = TRUE)
  shards <- rlang::ensyms(...)
  for (i in seq_along(shards)) {
    shardfile <- cache$file_return_key(rlang::expr_text(shards[[i]]))
    file.link(shardfile, file.path(dir, paste0(i, ".fst")))
  }
  disk.frame::disk.frame(dir)
}

combine_dynamic_diskframe <- function(dd, dir = drake::drake_tempfile(), cache = drake::drake_cache()) {
  if (is.character(cache)) cache <- drake::drake_cache(cache)
  dir.create(dir, recursive = TRUE)
  for (i in seq_along(dd)) {
    shardfile <- file.path(cache$path_return, dd[i])
    stopifnot(file.exists(shardfile))
    file.link(shardfile, file.path(dir, paste0(i, ".fst")))
  }
  disk.frame::disk.frame(dir)
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

#### Works like readd, but takes multiple arguments and makes a list

gett_list <- function(..., cache = drake::drake_cache()) {
  if (is.character(cache)) cache <- drake::drake_cache(cache)
  out <- lapply(
    as.character(match.call(expand.dots = FALSE)$...),
    cache$get
  )
  names(out) <- as.character(match.call(expand.dots = FALSE)$...)
  out
}

readd_list <- function(..., cache = drake::drake_cache()) {
  if (is.character(cache)) cache <- drake::drake_cache(cache)
  out <- lapply(
    as.character(match.call(expand.dots = FALSE)$...),
    readd,
    character_only = TRUE,
    cache = cache
  )
  names(out) <- as.character(match.call(expand.dots = FALSE)$...)
  out
}

readd_c <- function(...) {
  call <- match.call()
  call[[1]] <- readd_list
  do.call(c, eval(call))
}

readd_bind_rows <- function(...) {
  call <- match.call()
  call[[1]] <- readd_list
  dplyr::bind_rows(readd_list(...))
}
