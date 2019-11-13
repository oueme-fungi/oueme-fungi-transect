# Function(s) to streamline calling (m)locarna from R

#' Realign the nonconserved parts of an Infernal alignment with mLocarna
#'
#' @param alignment (\code{character} scalar) path to the alignment, in 
#'        locarna's extended clustal format
#' @param target_dir (\code{character} scalar) path to a directory to put the 
#'        output in.  Locarna will issue a warning if the directory already
#'        exists.
#' @param guide_tree (\code{character} scalar) path to a strictly bifurcating
#'        guide tree in Newick format.  All the sequence names in
#'        \code{alignment} should exist as tip labels in \code{guide_tree},
#'        but additional tips are allowed.
#' @param probabilistic (\code{logical} scalar) if \code{TRUE}, run Locarna in
#'        probabalistic mode.  This is much slower, but should yield more
#'        reliable results.
#' @param extended_pf (\code{logical} scalar) if \code{TRUE}, use extended
#'        precision calculations in probabalistic mode.
#' @param stockholm (\code{logical} scalar) if \code{TRUE}, write results in
#'        stockholm format in addition to clustal format.
#' @param consensus_structure (\code{character} scalar) how to calculate
#'        consensus structure to include in stockholm and screen output.
#' @param skip_pp (\code{logical} scalar) if \code{TRUE}, skip calculating
#'        pair probabilities if they already exist in the output directory.
#' @param cache_dir (\code{character} scalar) path to a directory for caching
#'        intermediat calculations.
#' @param verbose (\code{logical} scalar) if \code{TRUE}, print extra output.
#' @param quiet (\code{logical} scalar) if \code{TRUE}, don't print output to
#'        screen.
#' @param cpus (\code{integer} scalar) number of CPUs to use for parallel
#'        computation.  The progressive alignment stage of mlocarna is not
#'        effectively parallelized, but the initial all-to-all stages of the
#'        alignment (pair probabilities, pairwise alignments if no guide tree
#'        is given, pairwise probabilities in probabilistic mode, etc.) can be
#'        run in parallel.
#'
#' @return md5sum of the output
#' @export
mlocarna_realign <- function(alignment,
                             target_dir,
                             guide_tree = NULL,
                             probabilistic = FALSE,
                             extended_pf = FALSE,
                             stockholm = TRUE,
                             consensus_structure = c("none", "alifold", "mea"),
                             skip_pp = FALSE,
                             cache_dir = NULL,
                             only_dps = FALSE,
                             verbose = FALSE,
                             quiet = FALSE,
                             cpus = 1L,
                             pw_aligner = NULL,
                             pw_aligner_options = NULL,
                             pw_aligner_p = NULL,
                             pw_aligner_p_options = NULL) {
  assertthat::assert_that(
    assertthat::is.string(alignment),
    assertthat::is.readable(alignment),
    assertthat::is.string(target_dir),
    assertthat::is.flag(probabilistic),
    assertthat::is.flag(extended_pf),
    assertthat::is.flag(stockholm),
    is.character(consensus_structure),
    assertthat::is.flag(skip_pp),
    assertthat::is.flag(only_dps),
    assertthat::is.flag(verbose),
    assertthat::is.flag(quiet),
    !(verbose && quiet),
    assertthat::is.count(cpus)
  )
  args <- c("--realign", alignment, "--tgtdir", target_dir)
  
  if (!is.null(guide_tree)) {
    assertthat::assert_that(
      assertthat::is.string(guide_tree),
      assertthat::is.readable(guide_tree)
    )
    args <- c(args, "--treefile", guide_tree)
  }
  
  if (isTRUE(probabilistic)) args <- c(args, "--probabilistic")
  if (isTRUE(extended_pf)) {
    if (!isTRUE(probabilistic)) stop("'extended_pf' requires 'probabilistic'.")
    args <- c(args, "--extended-pf")
  }
  
  if (isTRUE(stockholm)) args <- c(args, "--stockholm")
  
  consensus_structure <- match.arg(consensus_structure)
  args <- c(args, "--consensus-structure", consensus_structure)
  
  if (isTRUE(skip_pp)) args <- c(args, "--skip-pp")
  if (isTRUE(only_dps)) args <- c(args, "--only-dps")
  if (isTRUE(verbose)) args <- c(args, "--verbose")
  if (isTRUE(quiet)) args <- c(args, "--quiet")
  
  if (!missing(cache_dir)) {
    assertthat::assert_that(
      assertthat::is.string(cache_dir)
    )
    args <- c(args, "--dp-cache", cache_dir)
  }
  if (!is.null(pw_aligner)) {
    assertthat::assert_that(
      assertthat::is.string(pw_aligner)
    )
    args <- c(args, "--pw-aligner", pw_aligner)
  }
  
  if (!is.null(pw_aligner_options)) {
    assertthat::assert_that(
      assertthat::is.string(pw_aligner_options)
    )
    args <- c(args, "--pw-aligner-options", paste0("\"", pw_aligner_options, "\""))
  }
  
  if (!is.null(pw_aligner_p)) {
    assertthat::assert_that(
      assertthat::is.string(pw_aligner_p)
    )
    args <- c(args, "--pw-aligner-p", pw_aligner_p)
  }
  
  if (!is.null(pw_aligner_p_options)) {
    assertthat::assert_that(
      assertthat::is.string(pw_aligner_p_options)
    )
    args <- c(args, "--pw-aligner-p-options", paste0("\"", pw_aligner_p_options, "\""))
  }
  
  if (!missing(cpus)) args <- c(args, "--cpus", cpus)
  cat("mlocarna", args, sep = " ")
  system2("mlocarna", args = args)
  outfile <- file.path(target_dir, "results", "result.aln")
  # return md5 of the output file
  if (!only_dps) return(tools::md5sum(outfile))
  # if only calculating initial files, calculate md5 of all of them and then
  # md5 of the list
  seqnames <- readLines(alignment)[-1]
  seqnames <- sub(" .*", "", seqnames)
  seqnames <- seqnames[!grepl("^#", seqnames)]
  seqnames <- seqnames[nchar(seqnames) > 0]
  seqnames <- unique(seqnames)
  ppdir <- file.path(target_dir, "input", seqnames)
  md5s <- tools::md5sum(ppdir)
  return(digest::digest(md5s, "md5"))
}

seqhash_safe <- function(seq) {
  unq <- unique(seq)
  h <- tzara::seqhash(unq, algo = "xxhash32")
  n <- 8
  while (length(unique(h)) != length(unq) && n <= 16) {
    h <- tzara::seqhash(unq, algo = "xxhash64", len = n)
    n <- n + 1
  }
  while (length(unique(h)) != length(unq) && n <= 128) {
    h <- tzara::seqhash(unq, algo = "sha512", len = n)
    n <- n + 1
  }
  h[match(unq, seq)]
}
