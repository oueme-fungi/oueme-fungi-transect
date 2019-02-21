
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
  
  if (region == "ITS") {
    positions <- gather(positions, key = "border", value = "loc", start, end) %>%
      dplyr::filter((border == "start" & region == "ITS1") |
               (border == "end" & region == "ITS2")) %>%
      mutate(region = "ITS") %>%
      spread(key = "border", value = "loc")
  }
  
  p <- dplyr::filter(positions, region == !!region,
              !is.na(start),
              start > 0,
              !is.na(end),
              end > 0,
              end <= readr::parse_number(length),
              end > start)
  fastq <- ShortRead::readFastq(infile)
  if (nrow(p)) {
    fastq <- fastq[p$idx] %>%
      ShortRead::narrow(start = p$start, end = p$end)
    ShortRead::writeFastq(fastq, outfile, mode = "a")
  }
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
