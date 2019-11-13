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
