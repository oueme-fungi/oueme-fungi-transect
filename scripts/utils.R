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
  )[dirname(files)]
  file.symlink(
    file.path(reldirs, files),
    file.path(to, files)
  )
}