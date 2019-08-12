join_seqs <- function(seq.tabs) {
  purrr::map(seq.tabs, as.data.frame) %>%
    purrr::map(tibble::rownames_to_column, "file") %>%
    dplyr::bind_rows() %>%
    dplyr::mutate_all(replace_na, 0L) %>%
    tibble::column_to_rownames("file") %>%
    as.matrix
}

taxonomy_dada <- function(seq.table, reference, multithread = FALSE) {
  assertthat::assert_that(file.exists(reference),
                          assertthat::is.readable(reference))
  # we can take a community matrix (in which case the sequences are the column
  # names) or the sequences
  if (is.matrix(seq.table)) {
    nreads <- Matrix::colSums(seq.table)
    seq <- colnames(seq.table)
  } else {
    nreads <- 1
    seq <- as.character(seq.table)
  }
  tax <- dada2::assignTaxonomy(seqs = seq, 
                               refFasta = reference,
                               multithread = multithread) %>%
    tibble::as_tibble(rownames = "seq") %>%
    # remove taxon rank prefixed from Unite reference
    dplyr::mutate_at(dplyr::vars(-seq), stringr::str_replace, "^[kpcofgs]__", "") %>%
    dplyr::mutate(
      # add Species to RDP reference
      Species = if ("Species" %in% names(.)) Species else NA_character_,
      Species = ifelse(is.na(Genus) | is.na(Species),
                       NA_character_,
                       paste(Genus, Species)),
      # Put the whole classification in one field for FUNGuild
      Taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species,
                       sep = ";") %>%
        stringr::str_replace_all(stringr::fixed(";NA"), ""),
      # Use the finest classification available as a short "name" for trees
      Name = dplyr::coalesce(Species, Genus, Family, Order, Class,
                             Phylum, Kingdom, "unknown") %>%
        stringr::str_replace_all("\\s", "_")) %>%
    # If there are duplicate names, number them.
    dplyr::group_by(Name) %>%
    dplyr::mutate(newname = if (dplyr::n() > 1) {
      paste(Name, seq_along(Name), sep = "_")
    } else {
      Name
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Name) %>%
    dplyr::rename(Name = newname) %>%
    
    # add in the reads if we had them.
    dplyr::left_join(tibble::tibble(seq = seq,
                                    nreads = nreads),
                     by = "seq")
}

its_join <- function(combined_map,
                     regions = c("ITS", "ITS1", "ITS2", "LSU", "long"),
                     joinregions = c("ITS1", "ITS2"),
                     maxdist = c(ITS = 20,
                                 ITS1 = 10,
                                 ITS2 = 10,
                                 LSU = 20,
                                 long = 40),
                     verbose = FALSE) {
  if (is.null(names(regions))) names(regions) <- regions
  if (is.null(names(joinregions))) names(joinregions) <- joinregions
  
  group_map <- combined_map %>%
    dplyr::group_by_at(regions) %>%
    dplyr::summarize(reads = n()) %>%
    dplyr::mutate(group = NA_integer_,
                  chimera = FALSE) %>%
    dplyr::filter(!(is.na(ITS1) & is.na(ITS2) & is.na(LSU))) %>%
    dplyr::arrange(dplyr::desc(reads)) %>%
    dplyr::ungroup()
  g <- 1
  for (i in 1:nrow(group_map)) {
    if (verbose) cat("sequence", i, "(", group_map$reads[i], "reads )\n")
    # to start with, make the next sequence its own group
    group_map$group[i] <- g
    # find all the candidates for merging
    cand <-
      purrr::map_dfr(
        joinregions,
        ~ dplyr::filter(group_map,
                        !chimera,
                        row_number() < i,
                        group_map[[.]] == group_map[[.]][i])) %>%
      unique()
    # if there are no matches, then let it continue as its own group.
    if (nrow(cand) == 0) {
      if (verbose) cat(" starting new group", g, "\n")
      g <- g + 1
    } else {
      # there are matches.  We need to make pairwise comparisons to ensure
      # there aren't problems.
      cand_groups <- unique(c(cand$group, g))
      cand_pairs <- tidyr::crossing(a = cand_groups,
                                    b = cand_groups,
                                    region = regions) %>%
        dplyr::filter(b > a) %>%
        dplyr::mutate(dmax =
                        purrr::pmap_dbl(
                          .,
                          function(a, b, region) {
                            aseq <- group_map[[region]][group_map$group == a] %>%
                              na.omit %>%
                              unique()
                            bseq <- group_map[[region]][group_map$group == b] %>%
                              na.omit %>%
                              unique()
                            if (length(aseq) == 0) return(NA_real_)
                            if (length(bseq) == 0) return(NA_real_)
                            fwd <- adist(aseq, bseq, partial = TRUE)
                            rev <- adist(bseq, aseq, partial = TRUE)
                            
                            return(max(pmin(fwd, t(rev))))
                          }),
                      dlimit = maxdist[region]) %>%
        dplyr::filter(stats::complete.cases(.))
      
      if (any(cand_pairs$dmax > cand_pairs$dlimit)) {
        if (verbose) {
          cat(" marking as chimera:\n")
          print(cand_pairs)
        }
        group_map$chimera[i] <- TRUE
        group_map$group[i] <- NA_integer_
      } else {
        if (verbose) {
          cat(" joining groups", paste(cand_groups, collapse = ", "), "\n")
          cand_pairs %>%
            dplyr::group_by(region) %>%
            dplyr::summarize(dmax = max(dmax),
                             dlimit = max(dlimit)) %>%
            print
        }
        group_map$group[group_map$group %in% cand_groups] <- min(cand_groups)
        cand_groups <- setdiff(cand_groups, min(cand_groups))
        cand_groups <- sort(cand_groups, decreasing = TRUE)
        for (j in cand_groups) {
          group_map %<>% mutate(
            group = ifelse(group > j,
                           group - 1,
                           group))
        }
      }
    }
  }
  return(group_map)
}

write_big_fasta <- function(big_seq_table, filename) {
  if (!dir.exists(dirname(filename))) dir.create(dirname(filename))
  tibble::as_tibble(big_seq_table, rownames = "filename") %>%
    tidyr::gather(key = "seq", value = "size", -1) %>%
    dplyr::filter(size >= 1) %>%
    tidyr::extract(col = "filename", into = c("tech", "run", "plate", "well", "region", "dir"), regex = "([a-z]+)_(\\d+)_(\\d+)-([A-H]1?[0-9])([fr]?)-([:alnum:]+).+") %>%
    dplyr::left_join(
      dplyr::group_by(., tech, run, plate, well) %>%
        dplyr::summarize(total = sum(size)),
      by = c("tech", "run", "plate", "well")) %>%
    dplyr::group_by(tech, run, plate, well, seq, total) %>%
    dplyr::summarize(size = sum(size)) %>%
    
    dplyr::mutate(f = size/total,
                  hash = tzara::seqhash(seq),
                  header = glue::glue("{hash};size={size};sample={tech}_{run}_{plate}{well};")) %>%
    dplyr::arrange(desc(f)) %$%
    Biostrings::DNAStringSet(magrittr::set_names(seq, header)) %T>%
    Biostrings::writeXStringSet(filepath = filename, compress = "gzip")
}
