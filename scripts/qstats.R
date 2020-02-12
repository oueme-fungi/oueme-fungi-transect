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
  # for floating points, we need to make bins so that there aren't
  # too many unique values.  Do transforms to preserve relevant distinctions
  # first
  out <- dplyr::mutate(
    out,
    eexp = round(log(eexp), 2),
    eexp = exp(eexp),
    erate = round(log(erate) - log1p(-erate), 2),
    erate = exp(erate)/(1 + exp(erate)),
    p.noerr = round(log(p.noerr) - log1p(-p.noerr), 2),
    p.noerr = exp(p.noerr)/(1 + exp(p.noerr)),
    n = NA_real_
  )
  out <- tidyr::pivot_longer(
    out,
    cols = c("length", "minq", "eexp", "erate", "p.noerr", "n"),
    names_to = "stat",
    values_to = "value"
  )
  out <- dplyr::group_by_all(out)
  out <- dplyr::summarize(out, nreads = dplyr::n())
  out <- dplyr::ungroup(out)
  
  if (nrow(out) == 0) {
    out <- tibble::tibble(
      ...,
      stat = c("length", "minq", "eexp", "erate", "p.noerr"),
      value = NA_real_,
      nreads = 0
    )
  }
  out
}

q_stats.character <- function(sreadq, ..., qualityType = "FastqQuality") {
  assertthat::assert_that(
    assertthat::is.string(sreadq),
    file.exists(sreadq)
  )
  infile <- sreadq
  if (endsWith(sreadq, ".bam")) infile <- pipe(paste("samtools", "fastq", shQuote(sreadq)))
  fqs <- ShortRead::FastqStreamer(infile, n = 100000)
  on.exit(close(fqs))
  out <- tibble::tibble()
  while (length(fq <- ShortRead::yield(fqs, qualityType = qualityType))) {
    out <- dplyr::bind_rows(out, q_stats.ShortReadQ(fq, file = sreadq, ...))
    out <- dplyr::group_by_at(out, dplyr::vars(-nreads))
    out <- dplyr::summarize(out, nreads = sum(nreads))
  }
  out <- dplyr::bind_rows(out)
  
  if (nrow(out) == 0) {
    out <- tibble::tibble(
      ...,
      length = NA_integer_,
      minq = NA_integer_,
      eexp = NA_real_,
      erate = NA_real_,
      p.noerr = NA_real_,
      reads = 0
    )
  }
  out
}

parse_qstat <- function(d) {
  dplyr::mutate_at(d, "file", basename) %>%
    tidyr::extract(col = "file", into = c("seq_run", "plate", "well", "read"),
                   regex = "([piS][bsH][-_]\\d{3,4})[-_](\\d{3}|OT\\d)-?([A-H]1?\\d)?(?:_S\\d_L\\d{3})?(_R[12])?[rf]?(?:_\\d{3})?[.].*")
}
