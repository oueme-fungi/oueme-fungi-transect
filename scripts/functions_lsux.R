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
