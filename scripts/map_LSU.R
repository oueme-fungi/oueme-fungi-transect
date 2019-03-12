# hash a set of sequences
seqhash <- function(seq, len = 8) UseMethod("seqhash")
seqhash.character <- function(seq, len = 8) {
  h <- purrr::map_chr(seq, digest::digest)
  stringr::str_sub(h, end = len)
}
seqhash.XStringSet <- function(seq, len = 8) {
  seqhash.character(as.character(seq))
}


# Take a list of \code{ShortReadQ} and return the reads, ids, and expected error
# in \code{tibble} form
raw_reads <- function(..., filenames, max_ee = Inf) {
  purrr::map2_dfr(list(...), filenames,
          function(x, name) {
            if (!methods::is(x, "ShortReadQ")) return(tibble::tibble(
              Seq.Run = character(),
              Region = character(),
              seq.id = character(),
              seq = character(),
              ee = numeric()))
            name <-
              stringr::str_match(name,
                        "regions_([:alpha:]{2}_\\d{3})(\\d{3})([A-H]1?\\d)([rf]?)_([:alnum:]+)")
            name <- c(name)[-1]
            
            tibble::tibble(
              Seq.Run = name[1],
              Plate = name[2],
              Well = name[3],
              Direction = name[4],
              Region = name[5],
              seq.id = as.character(x@id),
              seq = as.character(x@sread),
              ee = rowSums(10^-(as(x@quality, "matrix")/10), na.rm = TRUE)) %>%
              dplyr::filter(ee <= max_ee) %>%
              dplyr::select(-ee)
          })
}


combine_bigmaps <- function(dadamap, rawdata) {
  purrr::map_dfr(dadamap, ~tibble::tibble(file = names(.), data = .)) %>%
    tidyr::extract(
      col = "file",
      into = c("Seq.Run", "Plate", "Well", "Direction", "Region"),
      regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.qfilt\\.fastq\\.gz") %>%
    tidyr::unnest(data) %>%
    dplyr::full_join(rawdata) %>%
    dplyr::group_by(seq.id) %>%
    dplyr::filter(any(!is.na(asv.idx))) %>%
    dplyr::mutate(seq = dplyr::coalesce(asv.seq, derep.seq, seq)) %>%
    dplyr::select(-derep.seq, -derep.idx, -asv.seq, -asv.idx) %>%
    tidyr::spread(key = Region, value = seq)
}

calculate_consensus <- function(seq, names, ncpus = 1) {
  seq <- rlang::set_names(seq, names)
  seq <- stats::na.omit(seq)
  if (length(seq) < 3) return(NA_character_)
  cat("Calculating consensus of", length(seq), "sequences...\n")
  tictoc::tic("total")
  on.exit(tictoc::toc())
  seq <- Biostrings::DNAStringSet(seq) %>%
    Biostrings::RNAStringSet()
  
  cat(" Aligning...\n")
  tictoc::tic("  alignment")
  aln <- DECIPHER::AlignSeqs(seq, processors = ncpus, verbose = FALSE)
  tictoc::toc()
  
  cat(" Removing outliers...\n")
  tictoc::tic("  outliers")
  outliers <- odseq::odseq(Biostrings::RNAMultipleAlignment(aln))
  cat("  -removed", sum(outliers), "/", length(outliers),
      "sequences as outliers.\n")
  aln <- aln[!outliers]
  tictoc::toc()
  
  cat(" Masking gaps...\n")
  tictoc::tic("  masking")
  aln <- aln %>%
    Biostrings::RNAMultipleAlignment() %>%
    Biostrings::maskGaps(min.fraction = 0.5, min.block.width = 1) %>%
    as("RNAStringSet")
  tictoc::toc()
  
  cat(" Calculating consensus...\n")
  tictoc::tic("  consensus")
  on.exit(tictoc::toc(), add = TRUE)
  DECIPHER::ConsensusSequence(aln,
                              threshold = 0.5,
                              ambiguity = TRUE,
                              ignoreNonBases = TRUE,
                              includeTerminalGaps = FALSE)
}
