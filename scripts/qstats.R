q_stats <- function(sreadq, ...) {
  UseMethod("q_stats")
}


q_stats.ShortReadQ <- function(sreadq, ...){
  q <- as(sreadq@quality, "matrix")
  minq <- matrixStats::rowMins(q, na.rm = TRUE)
  q <- 10^(-q/10)
  out <- tibble::tibble(
    ...,
    length = ShortRead::width(sreadq),
    minq = minq,
    eexp = Matrix::rowSums(q, na.rm = TRUE),
    erate = eexp/length,
    p.noerr = exp(Matrix::rowSums(log1p(-q), na.rm = TRUE))
  )
  if (nrow(out) == 0) {
    out <- tibble::tibble(
      ...,
      length = NA_integer_,
      minq = NA_integer_,
      eexp = NA_real_,
      erate = NA_real_,
      p.noerr = NA_real_
    )
  }
  out
}

q_stats.character <- function(sreadq, ..., qualityType = "FastqQuality") {
  assertthat::assert_that(
    assertthat::is.string(sreadq),
    file.exists(sreadq)
  )
  fqs <- ShortRead::FastqStreamer(sreadq, n = 10000)
  on.exit(close(fqs))
  out <- list()
  while (length(fq <- ShortRead::yield(fqs, qualityType = qualityType))) {
    out <- c(out, list(
      q_stats.ShortReadQ(fq, file = sreadq, ...)
    ))
  }
  out <- dplyr::bind_rows(out)
  
  if (nrow(out) == 0) {
    out <- tibble::tibble(
      ...,
      length = NA_integer_,
      minq = NA_integer_,
      eexp = NA_real_,
      erate = NA_real_,
      p.noerr = NA_real_
    )
  }
  out
}

parse_qstat <- function(d) {
  mutate_at(d, "file", basename) %>%
    tidyr::extract(col = "file", into = c("seq_run", "plate", "well"),
                   regex = "([pi][bs]_\\d{3})_(\\d{3})-?([A-H]1?\\d)?[rf]?[.].*")
}
