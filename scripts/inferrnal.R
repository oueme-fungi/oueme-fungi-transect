

#' Search for Covariance Models (CM) in a set of sequences.
#' 
#' This function calls "\code{cmsearch}" from Infernal.  Infernal must be installed.
#' Many parameters are not included (yet!), and the function is focused on
#' retrieving the hits table and, optionally, producing an alignment.
#'
#' @param cm (character of length 1) Path to the covariance model (.cm) file.
#'  The covariance model must include calibration data from running 
#'  "\code{cmcalibrate}".
#' @param seq (filename, character vector, or 
#'   \code{\link[Biostrings]{XStringSet}}) Sequences to search with the CM.
#'   If a filename, the file can be of any type supported by Infernal.
#' @param glocal (logical of length 1) Whether to run the search in glocal mode
#'   (global with respect to the model, local with respect to the sequence).
#'   When \code{TRUE}, the search is faster, but will fail to find matches with
#'   only partially homologous structure.
#' @param alignment (filename) A file to save the aligned hits to.  If given,
#'   the alignment is saved in Stockholm format with annotations for secondary
#'   structure, posterior probablility, etc.
#' @param cpu (integer of length 1) The number of threads to use in the search.
#'
#' @return a \code{\link[tibble]{tibble}} with columns:
#'   \itemize{
#'     \item{target_name}{ (character) the name of the target sequence}
#'     \item{taget_accession}{(character) the target's accession number}
#'     \item{query_name}{(character) the name of the query CM}
#'     \item{query_accession}{(character) the query CM's accession number}
#'     \item{mdl}{(character) the model type ("cm" or "hmm")}
#'     \item{mdl_from}{(integer) the start location of the hit in the model}
#'     \item{mdl_to}{(integer) the end location of the hit in the model}
#'     \item{seq_from}{(integer) the start location of the hit in the sequence}
#'     \item{seq_to}{(integer) the end location of the hit in the sequence}
#'     \item{strand}{(character) the strand the hit was found on ("+" or "-")}
#'     \item{trunc)}{(character) whether the hit is truncated, and where ("no", "5'", "3'", "5'&3'", or "-" for hmm hits).}
#'     \item{pass}{(integer) which algorithm pass the hit was found on.}
#'     \item{gc}{(numeric) GC content of the hit}
#'     \item{bias}{(numeric) biased composition correction.  See the Infernal documentation.}
#'     \item{score}{(numeric) bit-score of the hit, including the biased
#'     composition correction.}
#'     \item{E_value}{(numeric) Expectation value for the hit.}}
#'     \item{inc}{(character) "!" if the sequence meets the inclusion threshold,
#'     "?" if it only meets the reporting threshold.}
#' @export
cmsearch <- function(cm,
                     seq,
                     glocal = TRUE,
                     alignment,
                     cpu) {
  assertthat::assert_that(assertthat::is.string(cm),
                          file.exists(cm),
                          assertthat::is.flag(glocal))
  tablefile <- tempfile("cmsearch", fileext = ".dat")
  on.exit(unlink(tablefile))
  args <- c("--tblout", tablefile, "--toponly", "--noali")
  if (isTRUE(glocal)) args <- c(args, "-g")
  
  if (!missing(cpu)) {
    assertthat::assert_that(assertthat::is.count(cpu))
    args <- c(args, "--cpu", cpu)
  }
  
  if (!missing(alignment)) {
    assertthat::assert_that(assertthat::is.string(alignment))
    d <- dirname(alignment)
    if (nchar(d) > 0 && !dir.exists(d)) dir.create(d, recursive = TRUE)
    args <- c(args, "--A", alignment)
  }
  
  args <- c(args, cm)
  
  seqfile <- NULL
  if (assertthat::is.string(seq) && file.exists(seq)) {
    seqfile <- seq
  } else {
    seqfile <- tempfile("seq", fileext = ".fasta")
    if (is.character(seq)) {
      seq <- Biostrings::BStringSet(seq)
      abc <- Biostrings::uniqueLetters(seq)
      if (all(abc %in% Biostrings::DNA_ALPHABET)) {
        seq <- Biostrings::DNAStringSet(seq)
        seq <- Biostrings::RNAStringSet(seq)
      } else if (all(abc %in% Biostrings::RNA_ALPHABET)) {
        seq <- Biostrings::RNAStringSet(seq)
      } else stop("Sequence alphabet should be DNA or RNA for CMalign.")
    }
    if (methods::is(seq, "XStringSet")) {
      Biostrings::writeXStringSet(seq, seqfile)
      on.exit(unlink(seqfile))
    } else  {
      stop("'seq' should be a filename, XStringSet, or character vector.")
    }
  }
  args <- c(args, seqfile)
  
  system2("cmsearch", args)
  
  readr::read_table2(tablefile,
                     col_names = c("target_name", "target_accession",
                                   "query_name", "query_accession",
                                   "mdl", "mdl_from", "mdl_to",
                                   "seq_from", "seq_to", "strand",
                                   "trunc", "pass", "gc", "bias",
                                   "score", "E_value", "inc",
                                   "description"),
                     col_types = "ccccciiiicciddddcc",
                     comment = "#")
}

cmalign <- function(cmfile, seq, glocal = TRUE, cpu) {
  assertthat::assert_that(assertthat::is.readable(cmfile),
                          assertthat::is.flag(glocal))
  args <- "cmalign"
  if (isTRUE(glocal)) args <- c(args, "-g")
  if (!missing(cpu)) {
    assertthat::assert_that(assertthat::is.count(cpu))
    args <- c(args, "--cpu", cpu)
  }
  
  seqfile <- NULL
  if (assertthat::is.string(seq) && file.exists(seq)) {
    seqfile <- seq
  } else {
    seqfile <- tempfile("seq", fileext = ".fasta")
    if (is.character(seq)) {
      seq <- Biostrings::BStringSet(seq)
      abc <- Biostrings::uniqueLetters(seq)
      if (all(abc %in% Biostrings::DNA_ALPHABET)) {
        seq <- Biostrings::DNAStringSet(seq)
        seq <- Biostrings::RNAStringSet(seq)
      } else if (all(abc %in% Biostrings::RNA_ALPHABET)) {
        seq <- Biostrings::RNAStringSet(seq)
      } else stop("Sequence alphabet should be DNA or RNA for CMalign.")
    }
    if (methods::is(seq, "XStringSet")) {
      Biostrings::writeXStringSet(seq, seqfile)
      on.exit(unlink(seqfile))
    } else  {
      stop("'seq' should be a filename, XStringSet, or character vector.")
    }
  }
  args <- c(args, cmfile, seqfile)
  args <- paste(args, collapse = " ")
  alnpipe <- pipe(args)
  read_stockholm_msa(alnpipe)
}
