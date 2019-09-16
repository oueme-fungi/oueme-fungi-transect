

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
      } else if (all(abc %in% Biostrings::RNA_ALPHABET)) {
        seq <- Biostrings::RNAStringSet(seq)
      }
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

merge_5_8S <- function(itsx_result, csearch_result) {
  csearch_result <- dplyr::select(csearch_result, 
                                  seq = target_name,
                                  seq_from, seq_to) %>%
    dplyr::mutate(region = "5_8S", seq = as.integer(seq))
  
  # CMsearch is better at finding 5.8S than ITSx is at finding anything,
  # so csearch_result probably has hits in sequences that did not pass ITSx.
  # We don't have a comparable alternative to find LSU (it is so long that the
  # CM is much slower) so for now ignore those.
  out <- dplyr::full_join(itsx_result, csearch_result,
                          by = c("seq", "region"))
  
  # remove chimeras
  out <- out %>%
    dplyr::filter(is.na(comment) | 
                    stringr::str_detect(comment, "Chimer", negate = TRUE)) %>%
    dplyr::group_by(seq) %>%
    dplyr::filter(sum(region == "5_8S") <= 1) %>%
    dplyr::ungroup()
  
  # find the location difference between the ITSx HMM and the Rfam CM
  to_shift <- dplyr::filter(out, region == "5_8S") %$%
    as.integer(median(end - seq_to, na.rm = TRUE))
  from_shift <- dplyr::filter(out, region == "5_8S") %$%
    as.integer(median(start - seq_from, na.rm = TRUE))
  
  # When ITSx didn't find 5.8S, use the adjusted CM location
  out <-
    dplyr::mutate(out,
                  start = dplyr::coalesce(start, pmax(seq_from + from_shift, 1L)),
                  end = dplyr::coalesce(end, seq_to + to_shift))
  
  # Adjust ITS1 end and ITS2 start to match
  # This involves switching to "wide" format
  # Also remove comments about missing 5.8S if 5.8S has been found.
  out <- dplyr::full_join(dplyr::select(out, seq, length, comment, region, start) %>%
                            tidyr::spread(key = "region", value = "start"),
                          dplyr::select(out, seq, length, comment, region, end) %>%
                            tidyr::spread(key = "region", value = "end"),
                          by = c("seq", "length", "comment"),
                          suffix = c("_start", "_end")) %>%
    dplyr::mutate(ITS1_end = dplyr::coalesce(`5_8S_start` - 1L, ITS1_end),
                  ITS1_start = dplyr::coalesce(ITS1_start, 1L),
                  ITS2_start = dplyr::coalesce(`5_8S_end` + 1L, ITS2_start),
                  comment = ifelse(is.na(`5_8S_start`) | is.na(`5_8S_end`),
                                   comment,
                                   stringr::str_replace(comment,
                                                        "[^!]*5.8S[^!]*! *",
                                                        "")) %>%
                    dplyr::na_if(""))
  # Switch back to "long" format
  out <-
    dplyr::full_join(dplyr::select(out, seq, length, comment, ends_with("_start")) %>%
                       dplyr::rename_all(stringr::str_replace, "_start", "") %>%
                       tidyr::gather(key = "region", value = "start", -(1:3)),
                     dplyr::select(out, seq, length, comment, ends_with("_end")) %>%
                       dplyr::rename_all(stringr::str_replace, "_end", "") %>%
                       tidyr::gather(key = "region", value = "end", -(1:3)),
                     by = c("seq", "length", "comment", "region"))
  
  out
}

read_stockholm_rf <- function(stockholm) {
  assertthat::assert_that((assertthat::is.string(stockholm) &&
                             file.exists(stockholm)) || 
                            methods::is(stockholm, "connection"))
  f <- function(x, pos, acc) {
    stringr::str_match(x, "#=GC +RF +(.+)")[,2] %>%
      na.omit() %>%
      paste(collapse = "")
  }
  readr::read_lines_chunked(stockholm,
                            readr::AccumulateCallback$new(f, acc = ""))
}


parse_stockholm_msa_chunk <- function(x, pos, acc) {
  x <- stringr::str_match(x, "^([^#][^ ]*) +([^ ]+)$")[,2:3]
  x <- x[complete.cases(x),]
  for (i in 1:nrow(x)) {
    if (x[i,1] %in% names(acc)) {
      acc[[x[i,1]]] <- paste0(acc[[x[i,1]]], x[i,2])
    } else {
      acc[[x[i,1]]] <- x[i,2]
    }
  }
  acc
}

read_stockholm_msa <- function(stockholm) {
  assertthat::assert_that((assertthat::is.string(stockholm) &&
                             file.exists(stockholm)) || 
                            methods::is(stockholm, "connection"))
  
  seqs <- 
    readr::read_lines_chunked(stockholm,
                              readr::AccumulateCallback$new(parse_stockholm_msa_chunk, acc = list()))
  Biostrings::RNAMultipleAlignment(unlist(seqs))
}

map_position <- function(alignment, x) {
  if (length(x) == 1 && is.na(x)) return(rep(NA_integer_, nrow(alignment)))
  assertthat::assert_that(assertthat::is.count(x),
                          x >= 1,
                          methods::is(alignment, "MultipleAlignment"),
                          x <= ncol(alignment))
  if (x == 1L) return(rep(1L, nrow(alignment)))
  
  trimaln <- substr(alignment, 1L, x - 1L)
  gapcounts <- stringr::str_count(trimaln, "[.-]")
  pmax(x - gapcounts, 1L)
}

extract_rf_region <- function(rf, n, names) {
  assertthat::assert_that(assertthat::is.string(rf),
                          is.character(n),
                          all(nchar(n) == 1))
  if (!missing(names)) {
    assertthat::assert_that(is.character(names),
                            length(names) == length(n))
  }
  out <- stringr::str_locate(rf, paste0(n, ".*", n))
  if (!missing(names)) {
    rownames(out) <- names
  }
  out
}

LSUx <- function(stockholm, include_incomplete = FALSE) {
  assertthat::assert_that(assertthat::is.string(stockholm),
                          assertthat::is.readable(stockholm))
  rf <- read_stockholm_rf(stockholm)
  aln <- read_stockholm_msa(stockholm)
  
  limits <- extract_rf_region(rf, c(1:9, LETTERS[1:10]),
                              c("5_8S", paste0("LSU", 1:18)))
  
  out <- tibble::tibble(seq_name = names(aln@unmasked))
  for (region in rownames(limits)) {
    for (boundary in colnames(limits)) {
      out[[paste0(region, "_", boundary)]] <- map_position(aln, limits[region, boundary])
    }
  }
  out[["ITS2_start"]] <- out[["5_8S_end"]] + 1L
  out[["ITS2_end"]] <- out[["LSU1_start"]] - 1L
  for (i in 2L:18L) {
    prename <- paste0("LSU", i - 1, "_end")
    postname <- paste0("LSU", i, "_start")
    if (include_incomplete || 
        any(!is.na(out[[prename]]) & !is.na(out[[postname]]))) {
    out[[paste0("V", i, "_start")]] <- out[[prename]] + 1L
    out[[paste0("V", i, "_end")]] <- out[[postname]] - 1L
    }
  }
  out <- purrr::discard(out, ~all(is.na(.)))
  outhead <- out["seq_name"]
  outvals <- dplyr::select(out, -seq_name)
  outvals <- outvals[order(apply(outvals, 2, median, na.rm = TRUE))]
  dplyr::bind_cols(outhead, outvals)
}

gather_regions <- function(pos) {
  starts <- dplyr::select(pos, -tidyselect::ends_with("_end"))
  starts <- tidyr::gather(starts, key = "region", value = "start", tidyselect::ends_with("_start"))
  starts <- dplyr::mutate_at(starts, "region", stringr::str_replace, "_start$", "")
  
  ends <- dplyr::select(pos, -tidyselect::ends_with("_start"))
  ends <- tidyr::gather(ends, key = "region", value = "end", tidyselect::ends_with("_end"))
  ends <- dplyr::mutate_at(ends, "region", stringr::str_replace, "_end$", "")
  
  hvars <- names(pos)
  hvars <- purrr::discard(hvars, endsWith, "_start")
  hvars <- purrr::discard(hvars, endsWith, "_end")
  joinvars <- c(hvars, "region")
  out <- dplyr::full_join(starts, ends, by = joinvars)
  dplyr::arrange(out, !!!rlang::parse_exprs(hvars), start)
}

spread_regions <- function(pos) {
  hvars <- setdiff(names(pos), c("region", "start", "end"))
  
  starts <- dplyr::select(pos, -end)
  starts <- dplyr::mutate_at(starts, "region", paste0, "_start")
  starts <- tidyr::spread(starts, key = "region", value = "start")
  
  ends <- dplyr::select(pos, -start)
  ends <- dplyr::mutate_at(ends, "region", paste0, "_end")
  ends <- tidyr::spread(ends, key = "region", value = "end")
  
  out <- dplyr::full_join(starts, ends, by = hvars)
  outhead <- out[hvars]
  outvals <- dplyr::select(out, -!!hvars)
  outvals <- outvals[order(apply(outvals, 2, median))]
  dplyr::bind_cols(outhead, outvals)
}
