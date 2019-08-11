q_stats <- function(file) {
  fq <- ShortRead::readFastq(file)
  q <- as(fq@quality, "matrix")
  minq <- matrixStats::rowMins(q, na.rm = TRUE)
  q <- 10^(-q/10)
  tibble::tibble(file = basename(file),
                 length = ShortRead::width(fq),
                 minq = minq,
                 eexp = Matrix::rowSums(q, na.rm = TRUE),
                 erate = eexp/length,
                 p.noerr = exp(Matrix::rowSums(log1p(-q), na.rm = TRUE)))
}
