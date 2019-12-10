q_stats <- function(sreadq, ...) {
  UseMethod("q_stats")
}


q_stats.ShortReadQ <- function(sreadq, ...){
  q <- as(sreadq@quality, "matrix")
  minq <- matrixStats::rowMins(q, na.rm = TRUE)
  q <- 10^(-q/10)
  tibble::tibble(...,
                 length = ShortRead::width(sreadq),
                 minq = minq,
                 eexp = Matrix::rowSums(q, na.rm = TRUE),
                 erate = eexp/length,
                 p.noerr = exp(Matrix::rowSums(log1p(-q), na.rm = TRUE)))
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
  dplyr::bind_rows(out)
}
