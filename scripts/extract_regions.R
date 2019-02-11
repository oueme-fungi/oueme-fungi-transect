
extract_region <- function(infile, outfile, region, positions) {
  assert_that(is.string(infile),
             file.exists(infile),
             is.string(outfile),
             is.string(region))
  
  #create the output directory if needed
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  assert_that(dir.exists(dirname(outfile)))
  
  # if the region is "full", then we don't need to cut anything.
  if (region %in% c("full", "long", "short")) {
    return(file.copy(from = infile, 
              to = outfile,
              overwrite = TRUE))
  }
  
  # make sure the file exists even if we don't have anything to write.
  if (file.exists(outfile)) file.remove(outfile)
  ShortRead::writeFastq(ShortRead::ShortReadQ(), outfile)
  
  p <- dplyr::filter(positions, region == !!region,
              !is.na(start),
              !is.na(end))
  fastq <- ShortRead::readFastq(infile)
  fastq <- fastq[p$idx] %>%
    ShortRead::narrow(start = p$start, end = p$end)
  ShortRead::writeFastq(fastq, outfile, mode = "a")
}

q_stats <- function(file) {
  fq <- ShortRead::readFastq(file)
  q <- as(fq@quality, "matrix")
  minq <- matrixStats::rowMins(q, na.rm = TRUE)
  q <- 10^(-q/10)
  tibble(file = basename(file),
         length = ShortRead::width(fq),
         minq = minq,
         eexp = rowSums(q, na.rm = TRUE),
         erate = eexp/length,
         p.noerr = exp(rowSums(log1p(-q), na.rm = TRUE)))
}
