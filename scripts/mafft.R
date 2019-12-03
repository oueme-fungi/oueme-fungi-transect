mafft_add <- function(x, y, add = c("add", "addfragments"),
                      method = c("localpair", "globalpair", "lastpair", "10merpair", "6merpair"),
                      maxiterate = c("0", "2", "1000"), threads = -1L,
                      quiet = FALSE, exec = Sys.which("mafft")) {
  assertthat::assert_that(
    methods::is(x, "DNAStringSet"),
    methods::is(y, "DNAStringSet"),
    rlang::is_scalar_integerish(threads),
    threads >= -1L,
    assertthat::is.flag(quiet)
  )
  x_file <- tempfile(pattern = "mafft", fileext = ".fasta")
  Biostrings::writeXStringSet(x, x_file)
  on.exit(unlink(x_file))
  y_file <- tempfile(pattern = "mafft", fileext = ".fasta")
  Biostrings::writeXStringSet(y, y_file)
  on.exit(unlink(y_file))
  
  add <- match.arg(add)
  method <- match.arg(method)
  maxiterate <- as.character(maxiterate)
  maxiterate <- match.arg(maxiterate)
  
  args <- c(
    paste0("--", add), y_file,
    paste0("--", method),
    "--maxiterate", maxiterate,
    "--thread", threads
  )
  if (isTRUE(quiet)) args <- c(args, "--quiet")
  
  out_file <- tempfile(pattern = "mafft", fileext = ".fasta")
  args <- c(args, x_file, ">", out_file)
  message(paste(exec, paste(args, collapse = " ")))
  system2(exec, args)
  on.exit(unlink(out_file))
  Biostrings::readDNAStringSet(out_file)
}
