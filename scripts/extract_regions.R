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
  q_stats.ShortReadQ(ShortRead::readFastq(sreadq, qualityType = qualityType), file = sreadq, ...)
}
