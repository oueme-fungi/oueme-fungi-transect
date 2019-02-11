dadamap <- function(derep, asv) {
  map2(derep, asv,
       function(derep, asv) {
         m <- tibble(seq.idx = seq_along(derep$map),
                     derep.idx = derep$map,
                     derep.seq = names(derep$uniques)[derep.idx])
         m %<>%
           left_join(
             tibble(
               asv.idx = asv$map,
               derep.idx = seq_along(asv.idx),
               asv.seq = asv$sequence[asv.idx]))
         m
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
  assert_that(file.exists(reference),
              is.readable(reference))
  
  tax <- seq.table %>%
    colnames %>%
    assignTaxonomy(reference, multithread = multithread) %>%
    as_tibble(rownames = "seq") %>%
    # remove taxon rank prefixed from Unite reference
    mutate_at(vars(-seq), str_replace, "^[kpcofgs]__", "") %>%
    mutate(
      # add Species to RDP reference
      Species = if ("Species" %in% names(.)) Species else NA_character_,
      Species = ifelse(is.na(Genus) | is.na(Species),
                       NA_character_,
                       paste(Genus, Species)))
  
  seq.table %>%
    t %>%
    as_tibble(rownames = "seq") %>%
    left_join(tax, by = "seq") %>%
    mutate(Taxonomy = paste(Kingdom, Phylum, Class, Order,
                            Family, Genus, Species,
                            sep = ";") %>%
             str_replace_all(fixed(";NA"), ""))
}
