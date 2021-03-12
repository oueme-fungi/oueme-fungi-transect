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
      stat = "n",
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
  if (endsWith(sreadq, ".sam.genozip")) infile <- pipe(paste("genocat", shQuote(sreadq), "| samtools fastq" ))
  fqs <- ShortRead::FastqStreamer(infile, n = 100000)
  on.exit(close(fqs))
  out <- tibble::tibble()
  while (length(fq <- ShortRead::yield(fqs, qualityType = qualityType))) {
    out <- dplyr::bind_rows(out, q_stats.ShortReadQ(fq, file = sreadq, ...))
    out <- dplyr::group_by_at(out, dplyr::vars(-nreads))
    out <- dplyr::summarize(out, nreads = sum(nreads))
  }
  out <- dplyr::ungroup(out)

  if (nrow(out) == 0) {
    out <- tibble::tibble(
      ...,
      stat = "n",
      value = NA_real_,
      nreads = 0
    )
  }
  out
}

parse_qstat <- function(d) {
  dplyr::mutate_at(d, "file", basename) %>%
    tidyr::extract(col = "file", into = c("seq_run", "plate", "well", "read"),
                   regex = "([piS][bsH][-_]\\d{3,4})[-_](\\d{3}|OT\\d)-?([A-H]1?\\d)?(?:_S\\d_L\\d{3})?(_R[12])?[rf]?(?:_\\d{3})?(?:[.].+)?")
}

compile_bioinf_table <- function(rawcounts, demuxcounts, filtercounts_full,
                                 regioncounts, filtercounts_ITS2, region_table,
                                 asv_table, allseqs, datasets) {
  bind_rows(
    enframe(rawcounts, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "Raw"
      ),
    enframe(demuxcounts, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "Trim"
      ),
    enframe(filtercounts_full, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "Filter (full)"
      ),
    enframe(regioncounts, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "LSUx"
      ),
    enframe(filtercounts_ITS2, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "Filter (ITS2)"
      ),
    filter(region_table, region == "ITS2") %>%
      select(seq_run, reads, step = region, ASVs),
    pivot_longer(asv_table, -1, names_to = "seq_run", values_to = "reads") %>%
      filter(reads > 0) %>%
      mutate_at("seq", chartr, old = "T", new = "U") %>%
      inner_join(
        select(allseqs, seq = ITS2, long, short, ITS, LSU) %>%
          pivot_longer(-1, names_to = "region", values_to = "consensus") %>%
          filter(!is.na(consensus)) %>%
          select(-consensus) %>%
          unique(),
        by = "seq"
      ) %>%
      group_by(seq_run, region) %>%
      summarize(reads = sum(reads, na.rm = TRUE), ASVs = n()) %>%
      select(seq_run, reads, step = region, ASVs)#,
    # pivot_longer(otu_table, -1, names_to = "seq_run", values_to = "reads") %>%
    #   filter(reads > 0) %>%
    #   mutate_at("seq", chartr, old = "T", new = "U") %>%
    #   inner_join(
    #     select(readd(allseqs, cache = cache), seq = ITS2, long, short, ITS, LSU) %>%
    #       pivot_longer(-1, names_to = "region", values_to = "consensus") %>%
    #       filter(!is.na(consensus)) %>%
    #       select(-consensus) %>%
    #       unique(),
    #     by = "seq"
    #   ) %>%
    #   group_by(seq_run, region) %>%
    #   summarize(reads = sum(reads, na.rm = TRUE), OTUs = n()) %>%
    #   select(seq_run, reads, step = region, OTUs)
  ) %>%
    left_join(select(datasets, seq_run, tech, amplicon), by = "seq_run") %>%
    select(tech, amplicon, step, reads, ASVs) %>%
    pivot_longer(c("reads", "ASVs"), names_to = "type", values_to = "count") %>%
    filter(!is.na(count)) %>%
    mutate(
      step = factor(step, levels = c("Raw", "Trim", "Filter (full)", "LSUx", "Filter (ITS2)", "ITS2",
                                     "short", "ITS", "LSU", "long")),
      tech = factor(tech, levels = c("PacBio", "Ion Torrent", "Illumina")),
      amplicon = factor(stringr::str_to_title(amplicon), levels = c("Long", "Short")),
      type = factor(type, c("ASVs", "reads"))
    ) %>%
    arrange(tech, amplicon, type) %>%
    # filter(ifelse(amplicon == "Short", !step  %in% c("ITS", "LSU", "long"), step != "short")) %>%
    pivot_wider(
      id_cols = "step",
      names_from = c("tech", "amplicon", "type"),
      values_from = c("count")
    ) %>%
    arrange(step) %>%
    column_to_rownames("step")
}
