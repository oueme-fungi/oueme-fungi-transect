#' Combine consensus sequences and ASVs for different regions
#'
#' @param conseqs (list of \code{\link[tibble]{tibble}}) ASVs for one region 
#'   (in the column named by \code{conseq_key}) and the corresponding consensus
#'   sequences for a linked region, as well as the number of reads for the ASV 
#'   (in column \code{nreads}).
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
                           by = c(conseq_key, "nreads"))
  
  names(seq_tables) <- stringr::str_replace(names(seq_tables),
                                            seq_table_prefix,
                                            "")
  
  # make sure all the names are present in the consensus table, so that the
  # join will work
  for (n in names(seq_tables)) {
    if (!n %in% names(conseqs)) conseqs[[n]] <- NA_character_
  }
  
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
    tidyr::extract(col = "filename", into = c("tech", "run", "plate", "well", "region", "dir"), regex = "([a-z]+)_(\\d+)_(\\d+)-([A-H]1?[0-9])([fr]?)-([:alnum:]+).+") %>%
    dplyr::left_join(
      dplyr::group_by(., tech, run, plate, well) %>%
        dplyr::summarize(total = sum(size)),
      by = c("tech", "run", "plate", "well")) %>%
    dplyr::group_by(tech, run, plate, well, seq, total) %>%
    dplyr::summarize(size = sum(size)) %>%
    
    dplyr::mutate(f = size/total,
                  hash = tzara::seqhash(seq),
                  header = glue::glue("{hash};size={size};sample={tech}_{run}_{plate}{well};")) %>%
    dplyr::arrange(desc(f)) %$%
    Biostrings::DNAStringSet(magrittr::set_names(seq, header)) %T>%
    Biostrings::writeXStringSet(filepath = filename, compress = "gzip")
}

derepShortReadQ <- function(reads, n = 1e+06, verbode = FALSE, qualityType = "Auto") {
  UseMethod("derep")
}

derepShortReadQ.list <- function(reads, n = 1e+06, verbose = FALSE, qualityType = "Auto") {
  assertthat::assert_that(
    all(purrr::map_lgl(reads, methods::is, "ShortReadQ"))
  )
  dir <- tempdir()
  fnames <- tempfile(seq_along(reads), dir, ".fastq.gz")
  for (i in seq_along(reads)) {
    ShortRead::writeFastq(reads[[i]], fnames[i], compress = TRUE, qualityType = qualityType)
  }
  dada2::derepFastq(fnames, n = n, verbose = verbose, qualityType = qualityType)
}

derepShortReadQ.ShortReadQ <- function(reads, n = 1e+06, verbose = FALSE, qualityType = "Auto") {
  fname <- tempfile(1, fileext = ".fastq.gz")
  ShortRead::writeFastq(reads, fname, compress = TRUE, qualityType = qualityType)
  dada2::derepFastq(fname, n = n, verbose = verbose, qualityType = qualityType)
}

filterReads <- function(reads, maxLen = Inf, minLen = 0,
                        maxEE = Inf) {
  assertthat::assert_that(methods::is(reads, "ShortReadQ"))
  reads <- reads[ShortRead::width(reads) <= maxLen]
  reads <- reads[ShortRead::width(reads) >= minLen]
  ee <- rowSums(10 ^ (-1 * (methods::as(reads@quality, "matrix") / 10)),
                na.rm = TRUE)
  reads <- reads[ee <= maxEE]
  reads
}
