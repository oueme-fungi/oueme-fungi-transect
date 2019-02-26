dadamap <- function(derep, asv) {
  purrr::map2(derep, asv,
       function(derep, asv) {
         m <- tibble::tibble(seq.id = names(derep$map),
                             derep.idx = derep$map,
                             derep.seq = names(derep$uniques)[derep.idx])
         
         dplyr::left_join(m,
                          tibble::tibble(
                            asv.idx = asv$map,
                            derep.idx = seq_along(asv.idx),
                            asv.seq = asv$sequence[asv.idx]))
       })
}

join_seqs <- function(seq.tabs) {
  purrr::map(seq.tabs, as.data.frame) %>%
    purrr::map(tibble::rownames_to_column, "file") %>%
    dplyr::bind_rows() %>%
    dplyr::mutate_all(replace_na, 0L) %>%
    tibble::column_to_rownames("file") %>%
    as.matrix
}

taxonomy <- function(seq.table, reference, multithread = FALSE) {
  assertthat::assert_that(file.exists(reference),
                          assertthat::is.readable(reference))
  
  tax <- seq.table %>%
    colnames %>%
    dada2::assignTaxonomy(reference, multithread = multithread) %>%
    tibble::as_tibble(rownames = "seq") %>%
    # remove taxon rank prefixed from Unite reference
    dplyr::mutate_at(vars(-seq), str_replace, "^[kpcofgs]__", "") %>%
    dplyr::mutate(
      # add Species to RDP reference
      Species = if ("Species" %in% names(.)) Species else NA_character_,
      Species = ifelse(is.na(Genus) | is.na(Species),
                       NA_character_,
                       paste(Genus, Species))) %>%
    dplyr::mutate(Taxonomy = paste(Kingdom, Phylum, Class, Order,
                            Family, Genus, Species,
                            sep = ";") %>%
             stringr::str_replace_all(fixed(";NA"), "")) %>%
    dplyr::left_join(tibble(seq = colnames(seq.table),
                     nreads = Matrix::colSums(seq.table)),
              by = "seq")
}

its_join <- function(bigmaps) {
  combined_map <-
    purrr::map_dfr(bigmaps, ~tibble::tibble(file = names(.), data = .)) %>%
    tidyr::extract(
      col = "file",
      into = c("Seq.Run", "Plate", "Well", "Direction", "Region"),
      regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.qfilt\\.fastq\\.gz") %>%
    tidyr::unnest(data) %>%
    dplyr::group_by(seq.id) %>%
    dplyr::filter(any(!is.na(asv.idx))) %>%
    dplyr::mutate(seq = dplyr::coalesce(asv.seq, derep.seq)) %>%
    dplyr::select(-derep.seq, -derep.idx, -asv.seq, -asv.idx) %>%
    tidyr::spread(key = Region, value = seq)
  
  group_map <- combined_map %>%
    dplyr::group_by(ITS, ITS1, ITS2, long, LSU) %>%
    dplyr::summarize(reads = n()) %>%
    dplyr::mutate(group = NA_integer_,
                  chimera = FALSE) %>%
    dplyr::filter(!(is.na(ITS1) & is.na(ITS2) & is.na(LSU)))
  g <- 0
  ungrouped <- Inf
  regions <- c("ITS", "ITS1", "ITS2", "LSU") %>% rlang::set_names()
  # while there are ITS1 and ITS2 sequences that have not been
  # assigned to a group
  while (any(!(is.na(group_map$ITS1)
               & is.na(group_map$ITS2))
             & is.na(group_map$group)
             & !group_map$chimera)) {
    g <- g + 1
    # find the most abundant sequence which has not been assigned
    toassign <- which(is.na(group_map$group) & !group_map$chimera)
    core <- toassign[which.max(group_map$reads[toassign])]
    group_map$group[core] <- g
    # while we are still assigning sequences to this group
    while (ungrouped >
           (ungrouped <- sum(is.na(group_map$group))
            - sum(group_map$chimera))) {
      # get the sequences already assigned to this group
      ingroup <- purrr::map(regions,
                     ~ unique(stats::na.omit(group_map[[.]][group_map$group == g])))
      candidates <- purrr::map(regions,
                        ~group_map[[.]] %in% ingroup[[.]]) %>%
        purrr::reduce(or) %>%
        which
      cand <- purrr::map(regions,
                  ~unique(stats::na.omit(group_map[[.]][candidates])))
      chimeras <- purrr::map2(
        cand, ingroup,
        function(cand, ingroup) {
          if (length(ingroup) * length(cand) == 0) return(character(0))
          tidyr::crossing(cand = cand, ingroup = ingroup) %>%
            dplyr::mutate(
              dist = purrr::pmap_dbl(.,
                                     ~ Biostrings::stringDist(c(.x, .y)))) %>%
            dplyr::group_by(cand) %>%
            dplyr::summarize(dist = min(dist)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(dist > 5) %$%
            cand
        })
      for (r in regions) {
        group_map$chimera[group_map[[r]] %in% chimeras[[r]]] <- TRUE
      }
      for (r in regions) {
        group_map$group[group_map[[r]] %in% ingroup[[r]]
                        & !group_map$chimera] <- g
      }
    }
  }
  return(group_map)
}