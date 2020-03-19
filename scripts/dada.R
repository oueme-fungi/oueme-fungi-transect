#' Combine consensus sequences and ASVs for different regions
#'
#' @param conseqs (list of \code{\link[tibble]{tibble}}) ASVs for one region 
#'   (in the column named by \code{conseq_key}) and the corresponding consensus
#'   sequences for a linked region, as well as the number of reads for the ASV 
#'   (in column \code{nread}).
#' @param seq_tables (named list of character matrix, as returned by 
#'   \code{\link[dada2]{makeSequenceTable}}) ASV matrices for the same (and
#'   optionally additional) regions as given in \code{conseqs}.  The list names
#'   should be the same as the column names for the regions in \code{conseq},
#'   with the optional addition of a prefix.
#' @param conseq_key (character of length 1) The region whose ASVs are used to
#'   group the other regions in \code{conseqs}.
#' @param seq_table_prefix (character of length one) A regular expression which
#'   be removed from the names of \code{seq_tables}.  It does not actually
#'   have to be a prefix.
#' @param label_order (character) regions in order of priority for naming the
#'   ASV.
#'
#' @return A tibble with columns for each of the regions in \code{conseqs} and
#'  \code{seq_tables}, where each row corresponds to sequences which match the 
#'  same ASV for the region given by \code{conseq_key}.  ASVs which are not
#'  represented in the consensus sequences are given on their own rows, where all
#'  other sequences are \code{NA}.  A column named \code{hash} is also created,
#'  giving the hash of the highest priority (as determined by \code{label_order})
#'  sequence which is present.
#' @export
make_allseq_table <- function(conseqs, seq_tables,
                              conseq_key = "ITS2",
                              seq_table_prefix = "big_seq_table_") {
  conseqs <- purrr::reduce(conseqs,
                           dplyr::full_join,
                           by = c(conseq_key, "nread"))
  
  names(seq_tables) <- stringr::str_replace(names(seq_tables),
                                            seq_table_prefix,
                                            "")
  
  # make sure all the names are present in the consensus table, so that the
  # join will work
  for (n in names(seq_tables)) {
    if (!n %in% names(conseqs)) seq_tables[[n]] <- NULL
  }
  
  seq_tables <- purrr::compact(seq_tables)
  
  seq_tables <- purrr::imap(seq_tables,
                            ~ tibble::tibble(x = colnames(.x)) %>%
                              set_colnames(.y))
  
  conseqs <-
    purrr::reduce(seq_tables,
                  dplyr::full_join,
                  .init = conseqs)
  
  dplyr::mutate(conseqs,
                hash = tzara::seqhash(.data[[conseq_key]]) %>%
                    unname())
}

# write the sequences from a dada2-style sequence table to a fasta file for
# clustering in [uv]search. Sequences are identified by hash of the sequence,
# number of reads, and the sample.
write_big_fasta <- function(big_seq_table, filename) {
  if (!dir.exists(dirname(filename))) dir.create(dirname(filename))
  tibble::as_tibble(big_seq_table, rownames = "filename") %>%
    tidyr::gather(key = "seq", value = "size", -1) %>%
    dplyr::filter(size >= 1) %>%
    tidyr::extract(
      col = "filename",
      into = c("seq_run", "plate", "well", "region", "dir"),
      regex = "([a-zA-Z]+[-_]\\d+)_(\\d+)_([A-H]1?[0-9])([fr]?)_([:alnum:]+).+"
    ) %>%
    dplyr::left_join(
      dplyr::group_by(., seq_run, plate, well) %>%
        dplyr::summarize(total = sum(size)),
      by = c("seq_run", "plate", "well")) %>%
    dplyr::group_by(seq_run, plate, well, seq, total) %>%
    dplyr::summarize(size = sum(size)) %>%
    
    dplyr::mutate(
      f = size/total,
      hash = tzara::seqhash(chartr("T", "U", seq)),
      header = glue::glue("{hash};size={size};sample={seq_run}_{plate}{well};") %>%
        as.character()
    ) %>%
    dplyr::arrange(desc(f)) %$%
    Biostrings::DNAStringSet(magrittr::set_names(seq, header)) %T>%
    Biostrings::writeXStringSet(filepath = filename, compress = "gzip")
}

# Use a temp file to dereplicate ShortReadQ objects using dada
derepShortReadQ <- function(reads, n = 1e+06, verbose = FALSE, qualityType = "Auto") {
  UseMethod("derepShortReadQ")
}

derepShortReadQ.list <- function(reads, n = 1e+06, verbose = FALSE, qualityType = "Auto") {
  nonnulls <- which(purrr::map_lgl(reads, is.null))
  assertthat::assert_that(
    all(purrr::map_lgl(reads[nonnulls], methods::is, "ShortReadQ"))
  )
  dir <- tempdir()
  fnames <- tempfile(as.character(seq_along(reads[nonnulls])), dir, ".fastq.gz")
  on.exit(unlink(fnames))
  for (i in seq_along(nonnulls)) {
    ShortRead::writeFastq(reads[[nonnulls[i]]], fnames[i], compress = TRUE, qualityType = qualityType)
  }
  derep <- dada2::derepFastq(fnames, n = n, verbose = verbose, qualityType = qualityType)
  out <- vector("list", length(reads))
  out[nonnulls] <- derep
  out
}

derepShortReadQ.ShortReadQ <- function(reads, n = 1e+06, verbose = FALSE, qualityType = "Auto") {
  if (length(reads) == 0) return(NULL)
  fname <- tempfile("reads", fileext = ".fastq.gz")
  on.exit(unlink(fname))
  ShortRead::writeFastq(reads, fname, compress = TRUE, qualityType = qualityType)
  dada2::derepFastq(fname, n = n, verbose = verbose, qualityType = qualityType)
}

filterReads <- function(reads, maxLen = Inf, minLen = 0,
                        maxEE = Inf) {
  assertthat::assert_that(methods::is(reads, "ShortReadQ"))
  reads <- reads[ShortRead::width(reads) <= maxLen]
  reads <- reads[ShortRead::width(reads) >= minLen]
  reads <- reads[!grepl("N", as.character(reads@sread))]
  ee <- rowSums(10 ^ (-1 * (methods::as(reads@quality, "matrix") / 10)),
                na.rm = TRUE)
  reads <- reads[ee <= maxEE]
  reads
}

filterReadPairs <- function(reads1, reads2, trimR = 0, truncQ = 0, maxLen = Inf, minLen = 1,
                        maxEE = Inf) {
  assertthat::assert_that(
    methods::is(reads1, "ShortReadQ"),
    methods::is(reads2, "ShortReadQ"),
    is.numeric(trimR),
    length(trimR) %in% 1L:2L,
    is.numeric(truncQ),
    length(truncQ) %in% 1L:2L,
    is.numeric(maxLen),
    length(maxLen) %in% 1L:2L,
    is.numeric(minLen),
    length(minLen) %in% 1L:2L,
    is.numeric(maxEE),
    length(maxEE) %in% 1L:2L
  )
  
  if (!missing(trimR)) {
    if (length(trimR == 1)) trimR <- c(trimR, trimR)
    reads1 <- ShortRead::narrow(reads1, end = pmax(ShortRead::width(reads1) - trimR[1], 0))
    reads2 <- ShortRead::narrow(reads2, end = pmax(ShortRead::width(reads2) - trimR[1], 0))
  }
  
  if (!missing(truncQ)) {
    if (length(truncQ == 1)) truncQ <- c(truncQ, truncQ)
    w1 <- apply(methods::as(reads1@quality, "matrix") <= truncQ[1],
                MARGIN = 1, match, x = TRUE)
    reads1 <- ShortRead::narrow(reads1, end = w1 - 1)
    w2 <- apply(methods::as(reads2@quality, "matrix") <= truncQ[2],
                MARGIN = 1, match, x = TRUE)
    reads2 <- ShortRead::narrow(reads2, end = w2 - 1)
  }
  
  if (length(maxLen == 1)) maxLen <- c(maxLen, maxLen)
  shortenough <- ShortRead::width(reads1) <= maxLen[1] &
    ShortRead::width(reads2) <= maxLen[2]
  reads1 <- reads1[shortenough]
  reads2 <- reads2[shortenough]
  if (length(reads1) == 0) return(list(R1 = reads1, R2 = reads2))
  
  if (length(minLen == 1)) minLen <- c(minLen, minLen)
  longenough <- ShortRead::width(reads1) >= minLen[1] &
    ShortRead::width(reads2) >= minLen[1]
  reads1 <- reads1[longenough]
  reads2 <- reads2[longenough]
  if (length(reads1) == 0) return(list(R1 = reads1, R2 = reads2))
  
  noN <- !grepl("N", as.character(reads1@sread)) &
    !grepl("N", as.character(reads2@sread))
  reads1 <- reads1[noN]
  reads2 <- reads2[noN]
  if (length(reads1) == 0) return(list(R1 = reads1, R2 = reads2))
  
  if (length(maxEE == 1)) maxEE <- c(maxEE, maxEE)
  ee1 <- rowSums(10 ^ (-1 * (methods::as(reads1@quality, "matrix") / 10)),
                na.rm = TRUE)
  ee2 <- rowSums(10 ^ (-1 * (methods::as(reads2@quality, "matrix") / 10)),
                 na.rm = TRUE)
  reads1 <- reads1[ee1 <= maxEE[1] & ee2 <= maxEE[2]]
  reads2 <- reads2[ee1 <= maxEE[1] & ee2 <= maxEE[2]]
  list(R1 = reads1, R2 = reads2)
}

extract_and_derep <- function(positions, trim_file, region, region_start, region_end,
                              max_length, min_length, max_ee) {
  if (nrow(positions) == 0) return(NULL)
  pos <- dplyr::group_by(positions, trim_file)
  filekey <- dplyr::select(positions, "trim_file", "seq") %>%
    unique()
  regions <- 
    tzara::extract_region(
      seq = dplyr::group_keys(pos)$trim_file,
      region = region_start,
      region2 = region_end,
      positions = dplyr::group_split(pos)
    )
  qstats_region <- q_stats(
    regions,
    step = "lsux",
    file = plyr::mapvalues(
      as.character(regions@id),
      filekey$seq,
      filekey$trim_file,
      warn_missing = FALSE
    ),
    region = region
  )
  filter <- filterReads(
    regions,
    maxLen = max_length,
    minLen = min_length,
    maxEE = max_ee
  )
  qstats_filter <- q_stats(
    filter,
    step = "filter",
    file = plyr::mapvalues(
      as.character(filter@id),
      filekey$seq,
      filekey$trim_file,
      warn_missing = FALSE
    ),
    region = region
  )
  qstats <- dplyr::bind_rows(qstats_region, qstats_filter)
  if (nrow(qstats) == 0) {
    qstats <- tibble::tibble(
      file = dplyr::group_keys(pos)$trim_file,
      length = NA_integer_,
      minq = NA_integer_,
      eexp = NA_real_,
      erate = NA_real_,
      p.noerr = NA_real_
    )
  }
  if (length(filter) == 0) return(structure(list(), qstats = qstats))
  out <- 
    derepShortReadQ(
    reads = filter,
    n = 1e4,
    qualityType = "FastqQuality",
    verbose = TRUE
  )
  
  out[["names"]] <- as.character(filter@id)
  attr(out, "qstats") <- qstats
  out
}

filter_and_derep_pairs <- function(trim_file_1, trim_file_2, trimR, truncQ, max_length,
                                   min_length, max_ee, ID, ...) {
  reads1 <- ShortRead::readFastq(trim_file_1)
  reads2 <- ShortRead::readFastq(trim_file_2)
  
  filter_reads <- filterReadPairs(reads1, reads2, trimR, truncQ, max_length, min_length, max_ee)
  reads1 <- filter_reads$R1
  reads2 <- filter_reads$R2
  
  qstats <- 
    dplyr::bind_rows(
      q_stats(reads1, step = "filter", file = trim_file_1, read = "R1", ...),
      q_stats(reads2, step = "filter", file = trim_file_2, read = "R2", ...)
    )
  if (nrow(qstats) == 0) {
    qstats <- tibble::tibble(
      file = c(trim_file_1, trim_file_2),
      read = c("R1", "R2"),
      length = NA_integer_,
      minq = NA_integer_,
      eexp = NA_real_,
      erate = NA_real_,
      p.noerr = NA_real_
    )
  }
  if (length(reads1) == 0) return(structure(list(), qstats = qstats))
  out <- list()
  out$R1 <- derepShortReadQ(
      reads = reads1,
      n = 1e4,
      qualityType = "FastqQuality",
      verbose = TRUE
    )
  out$R1[["names"]] <- as.character(reads1@id)
  out$R2 <- derepShortReadQ(
    reads = reads2,
    n = 1e4,
    qualityType = "FastqQuality",
    verbose = TRUE
  )
  out$R2[["names"]] <- as.character(reads2@id)
  attr(out, "qstats") <- qstats
  out
}

# call dada, but be tolerant of NULL inputs and empty list inputs
robust_dada <- function(derep, ...) {
  UseMethod("robust_dada")
}

robust_dada.derep <- function(derep, ...) {
  dada2::dada(derep, ...)
}

robust_dada.list <- function(derep, ...) {
  nonnulls <- which(!vapply(derep, is.null, TRUE) & vapply(derep, length, 1L) > 0)
  assertthat::assert_that(
    all(vapply(derep[nonnulls], methods::is, TRUE, "derep"))
  )
  dada <- dada2::dada(derep[nonnulls], ...)
  if (methods::is(dada, "dada")) {
    dada <- list(dada)
  }
  out <- vector("list", length(derep))
  out[nonnulls] <- dada
  names(out) <- names(derep)
  out
}

robust_dada.character <- function(derep, ...) {
  existing <- which(!vapply(derep, file.exists), TRUE)
  dada <- dada2::dada(derep[existing], ...)
  out <- vector("list", length(derep))
  out[existing] <- dada
  out
}

multidada <- function(dereplist, dadalist, region, ..., keyvars = NULL) {
  dadamap <- tzara::dadamap(dereplist, dadalist, region = region, ...)
  if (nrow(dadamap) == 0) return(dadamap)
  regions <- unique(dadamap$region)
  dadamap %>%
    dplyr::select(-name, -derep.idx, -derep.seq, -dada.idx) %>%
    dplyr::mutate_at("dada.seq", chartr, old = "T", new = "U") %>%
    tidyr::spread(key = "region", value = "dada.seq") %>%
    dplyr::group_by_at(regions) %>%
    dplyr::summarize(nread = dplyr::n()) %>%
    dplyr::ungroup()
}
