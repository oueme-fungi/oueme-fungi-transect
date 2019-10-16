

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

  gc <- stringr::str_match(x, "#=GC +([^ ]+) +(.+)")[,2:3]
  gc <- gc[complete.cases(gc),]
  for (i in seq_len(nrow(gc))) {
    attr(acc, gc[i,1]) <- paste0(attr(acc, gc[i,1]), gc[i,2])
  }

  x <- stringr::str_match(x, "^(\\d+\\|)?([^#][^ ]*) +([^ ]+)$")[,3:4]
  x <- x[complete.cases(x),]
  for (i in seq_len(nrow(x))) {
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
  out <- attributes(seqs)
  out[["alignment"]] <- Biostrings::RNAMultipleAlignment(unlist(seqs))
  out
}

read_stockholm_msa_rf <- function(stockholm) {
  
}

map_position <- function(alignment, x) {
  if (length(x) == 1 && is.na(x)) return(rep(NA_integer_, nrow(alignment)))
  assertthat::assert_that(assertthat::is.count(x),
                          x >= 1,
                          methods::is(alignment, "MultipleAlignment"),
                          x <= ncol(alignment))
  if (x == 1L) return(rep(1L, nrow(alignment)))
  
  trimaln <- substr(alignment, 1L, x)
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


gap_free_width <- function(x, gapchars = ".-") {
  if (methods::is(x, "MultipleAlignment")) x <- x@unmasked
  assertthat::assert_that(methods::is(x, "XStringSet"))
  Biostrings::width(x) - c(Biostrings::letterFrequency(x, gapchars))
}

# extract_LSU <- function(stockholm, include_incomplete = FALSE) {
#   assertthat::assert_that(assertthat::is.string(stockholm),
#                           assertthat::is.readable(stockholm))
#   
#   aln <- read_stockholm_msa(stockholm)
#   rf <- aln$RF
#   aln <- aln$alignment

extract_LSU <- function(aln, rf, include_incomplete = FALSE) {
  limits <- extract_rf_region(rf, c(1:9, LETTERS[1:10]),
                              c("5_8S", paste0("LSU", 1:18)))
  
  outhead <- tibble::tibble(seq_name = names(aln@unmasked),
                        length = gap_free_width(aln@unmasked))
  
  out <- tibble::tibble(.rows = nrow(outhead))
  
  for (region in rownames(limits)) {
    startcol <- paste0(region, "_start")
    endcol <- paste0(region, "_end")
    
    out[[startcol]] <- map_position(aln, limits[region, "start"])
    out[[endcol]] <- map_position(aln, limits[region, "end"])
    
    out_of_range <- out[[startcol]] >= outhead$length & out[[endcol]] >= outhead$length
    out[[startcol]] <- ifelse(out_of_range, NA_integer_, out[[startcol]])
    out[[endcol]] <- ifelse(out_of_range, NA_integer_, out[[endcol]])
  }
  
  out[["ITS2_start"]] <- out[["5_8S_end"]] + 1L
  out[["ITS2_end"]] <- out[["LSU1_start"]] - 1L
  for (i in 2L:18L) {
    prename <- paste0("LSU", i - 1, "_end")
    postname <- paste0("LSU", i, "_start")
    if (isTRUE(include_incomplete) || 
        any(!is.na(out[[prename]]) & !is.na(out[[postname]]))) {
      zerosize <- out[[postname]] - out[[prename]] == 1L
      out[[paste0("V", i, "_start")]] <- ifelse(zerosize, NA_integer_,
                                                out[[prename]] + 1L)
      out[[paste0("V", i, "_end")]] <- ifelse(zerosize, NA_integer_,
                                              out[[postname]] - 1L)
    }
  }
  out <- purrr::discard(out, ~all(is.na(.)))
  out <- out[order(apply(out, 2, median, na.rm = TRUE))]
  dplyr::bind_cols(outhead, out)
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


LSUx <- function(seq, cm_5.8S, cm_32S, glocal = TRUE, ITS1 = FALSE, cpu) {
  assertthat::assert_that(assertthat::is.readable(cm_5.8S),
                          assertthat::is.readable(cm_32S),
                          assertthat::is.flag(glocal),
                          assertthat::is.flag(ITS1))
  
  futile.logger::flog.info("Beginning CM search.", name = "LSUx")
  cms <- cmsearch(cm = cm_5.8S, seq = seq, glocal = glocal, cpu = cpu)
  
  # remove multiple hits
  cms <- dplyr::group_by(cms, target_name)
  cms <- dplyr::filter(cms, all(inc == "?") | inc == "!")
  cms <- dplyr::filter(cms, dplyr::n() == 1)
  cms <- dplyr::ungroup(cms)
  futile.logger::flog.info("%d sequences contained a single 5.8S hit.", nrow(cms),
                           name = "LSUx")
  seq_32S <- IRanges::narrow(seq[cms$target_name], start = cms$seq_from)
  
  futile.logger::flog.info("Beginning CM alignment.", name = "LSUx")
  aln <- cmalign(cmfile = cm_32S, seq = seq_32S, glocal = glocal, cpu = cpu)
  
  futile.logger::flog.info("Extracting LSU regions.", name = "LSUx")
  pos <- extract_LSU(aln = aln$alignment, rf = aln$RF)
  pos <- dplyr::mutate_at(pos, "seq_name", stringr::str_replace, "^[^|]*\\|", "")
  pos <- dplyr::mutate_if(pos, is.integer, add, cms$seq_from - 1L)
  if (isTRUE(ITS1)) {
    no_ITS1 <- pos$`5_8S_start` == 1 | is.na(pos$`5_8S_start`)
    if (all(no_ITS1)) {
      futile.logger::flog.warn("ITS1 annotation was requested, but no bases before 5.8S were found.",
                name = "LSUx")
    } else {
      pos$ITS1_start <- ifelse(no_ITS1, NA_integer_, 1L)
      pos$ITS1_end <- ifelse(no_ITS1, NA_integer_, pos$`5_8S_start` - 1L)
      pos <- dplyr::select(pos, seq_name, length, ITS1_start, ITS1_end,
                           dplyr::everything())
    }
  }
  pos
}

write_clustalw_ss <- function(aln, sec_str, file, ref = sec_str, seq_names = names(aln), write_ref = TRUE) {
  assertthat::assert_that(methods::is(aln, "XStringSet"))
  assertthat::assert_that(is.character(seq_names),
                          length(seq_names) == length(aln),
                          assertthat::is.string(file),
                          assertthat::is.string(sec_str),
                          assertthat::is.string(ref),
                          assertthat::is.flag(write_ref),
                          all(nchar(sec_str) == Biostrings::width(aln)),
                          nchar(sec_str) == nchar(ref))
  
  aln <- Biostrings::RNAStringSet(aln)
  sec_str <- chartr("{[<>]},:_-", "((())).x..", sec_str)
  con <- file(file, "wt")
  on.exit(close(con))
  writeLines("CLUSTALW", con)
  start <- 1
  width <- unique(Biostrings::width(aln))
  namewidth <- max(nchar(seq_names))
  seq_names <- stringr::str_pad(seq_names, namewidth, "right")
  str_name <- stringr::str_pad("#S", namewidth, "right")
  ref_name <- stringr::str_pad("#A1", namewidth, "right")
  while (start <= width) {
    end <- min(start + 59, width)
    writeLines("", con)
    writeLines(paste(seq_names, substr(aln, start, end), end), con)
    writeLines(paste(str_name, substr(sec_str, start, end), end), con)
    if (!missing(ref) || write_ref)
    writeLines(paste(ref_name, substr(ref, start, end), end), con)
    start <- end + 1
  }
}

parse_clustal_ss_chunk <- function(x, pos, acc) {
  ss <- stringr::str_match(x, "#S *([()\\[\\]{}<>.,:_x-]+) ?\\d*$")
  ss <- ss[complete.cases(ss),]
  for (i in seq_len(nrow(ss))) {
    attr(acc, "SS_cons") <- paste0(attr(acc, "SS_cons"), ss[i, 2])
  }
  
  x <- stringr::str_match(x, "^([^# ][^ ]*) +([^ ]+) *\\d*$")[,2:3]
  x <- x[complete.cases(x),]
  for (i in seq_len(nrow(x))) {
    if (x[i,1] %in% names(acc)) {
      acc[[x[i,1]]] <- paste0(acc[[x[i,1]]], x[i,2])
    } else {
      acc[[x[i,1]]] <- x[i,2]
    }
  }
  acc
}

read_clustalw_ss <- function(clustal) {
  assertthat::assert_that((assertthat::is.string(clustal) &&
                             file.exists(clustal)) || 
                            methods::is(clustal, "connection"))
  
  seqs <- 
    readr::read_lines_chunked(
      clustal,
      readr::AccumulateCallback$new(parse_clustal_ss_chunk, acc = list())
    )
  
  #names(seqs) <- stringr::str_replace(names(seqs), "^\\d+\\|", "")
  out <- attributes(seqs)
  out[["alignment"]] <- Biostrings::RNAMultipleAlignment(unlist(seqs))
  out
}

remove_nonconsensus_nongaps <- function(aln, gapfrac = 1) {
  # rle gives a list with "value" and "length" for each run
  # turn it into a tibble, and calculate the begin and end of each run
  # TRUE = nonconsensus position, FALSE = consensus position
  rle(seqinr::s2c(aln$SS_cons) == ".") %>%
    unclass() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      end = cumsum(lengths),
      start = end - lengths + 1
    ) %>%
    # take only the consensus values
    dplyr::filter(!.data$values) %>%
    purrr::pmap(
      function(start, end, ...)
        substr(aln$alignment@unmasked, start, end)
    ) %>%
    # paste them together
    do.call(paste0, .) %>%
    magrittr::set_names(aln$names) %>%
    Biostrings::RNAMultipleAlignment() %>%
    # mask all-gap columns
    Biostrings::maskGaps(gapfrac, 1) %>%
    # converting to XStringSet removes masked columns
    methods::as("RNAStringSet")
}

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
#' @return the path to the resulting alignment
#' @export
mlocarna_realign <- function(alignment,
                             target_dir,
                             guide_tree = NULL,
                             probabilistic = FALSE,
                             extended_pf = FALSE,
                             skip_pp = FALSE,
                             cache_dir = NULL,
                             verbose = FALSE,
                             quiet = FALSE,
                             cpus = 1L) {
  assertthat::assert_that(
    assertthat::is.string(alignment),
    assertthat::is.readable(alignment),
    assertthat::is.string(target_dir),
    assertthat::is.flag(probabilistic),
    assertthat::is.flag(extended_pf),
    assertthat::is.flag(skip_pp),
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
  
  if (isTRUE(skip_pp)) args <- c(args, "--skip-pp")
  if (isTRUE(verbose)) args <- c(args, "--verbose")
  if (isTRUE(quiet)) args <- c(args, "--quiet")
  
  if (!missing(cache_dir)) {
    assertthat::assert_that(
      assertthat::is.string(cache_dir)
    )
    args <- c(args, "--dp-cache", cache_dir)
  }
  
  if (!missing(cpus)) args <- c(args, "--cpus", cpus)
  cat("mlocarna", args, sep = " ")
  system2("mlocarna", args = args)
  file.path(target_dir, "results", "result.aln")
}
