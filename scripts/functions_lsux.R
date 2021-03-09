# functions used in conjunction with LSUx
# most of this has been moved to the LSUx package.
# author Brendan Furneaux

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

trim_LSU_intron <- function(aln) {
  consensus <- DECIPHER::ConsensusSequence(
    aln,
    threshold = 0.5,
    ambiguity = FALSE
  )
  consensus <- as.character(consensus)
  site <- stringi::stri_locate_first_fixed(consensus, "CAAAUUUGGGUAUAG")[1, 'end']
  IRanges::narrow(aln, start = 1, end = site)
}

# add a column called "out_col" to a region table by concatenating the sequences
# from the columns "regions" in order.  If out_col is already a column in the
# table, its current value is used as a backup in case one of the "regions" is
# NA.  If "key_col" is given, then it is required that the sequence in "key_col"
# is a subsequence of "out_col", or "out_col" will be NA.
region_concat <- function(table, out_col, regions, key_col = NULL) {
  if (!out_col %in% names(table)) return(table)
  table[[out_col]] <- dplyr::coalesce(
    do.call(stringr::str_c, table[,regions]),
    table[[out_col]]
    )
  if (!is.null(key_col)) {
    table[[out_col]] <- ifelse(
      stringi::stri_detect_fixed(table[[out_col]], table[[key_col]]),
      table[[out_col]],
      NA_character_
    )
  }
  table
}

# as region_concat, but the final sequence is "seeded" with the key column,
# and extended in each direction only until it reaches an NA sequence
stepwise_region_concat <- function(table, out_col, regions, key_col) {
  # regions occuring before the key, in reverse order
  preregions <- rev(regions[!dplyr::cumany(regions == key_col)])
  # regions occurring after the key, in forward order
  postregions <- rev(rev(regions)[!dplyr::cumany(rev(regions) == key_col)])
  # temporary column; make sure its name is unique
  pre_col <- dplyr::last(make.names(c(names(table), "pre"), unique = TRUE))
  # initialize with the key column
  table[[pre_col]] <- table[[key_col]]
  table[[out_col]] <- table[[pre_col]]
  # extend the 5' end of the key column
  # pre_col will be NA after the first time an NA region is encountered.
  for (r in preregions) {
    table[[pre_col]] <- stringr::str_c(table[[r]], table[[pre_col]])
    table[[out_col]] <- dplyr::coalesce(table[[pre_col]], table[[out_col]])
  }
  # extend the 3' end of the key column
  # pre_col will be NA after the first time an NA region is encountered.
  table[[pre_col]] <- table[[out_col]]
  for (r in postregions) {
    table[[pre_col]] <- stringr::str_c(table[[pre_col]], table[[r]])
    table[[out_col]] <- dplyr::coalesce(table[[pre_col]], table[[out_col]])
  }
  table[[pre_col]] <- NULL
  table
}
