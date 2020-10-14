rank_factor <- function(r,
                        ranks = c("rootrank", "domain", "kingdom", "phylum",
                                  "class", "order", "family", "genus",
                                  "species"),
                        abbrev = FALSE) {
  if (abbrev) {
    factor(r, levels = substr(ranks, 1, 1), labels = ranks, ordered = TRUE)
  } else {
    factor(r, levels = ranks, ordered = TRUE)
  }
}

# combine a named list of taxonomy tables
combine_taxon_tables <- function(tables, allseqs) {
  tibble::enframe(tables) %>%
    tidyr::extract(
      col = name,
      into = c("region", "reference", "ref_region", "method"),
      regex = "taxon_(32S|5_8S|ITS[12]?|LSU|long|short)_(unite|warcup|rdp_train)_(ITS[12]?|LSU)_(idtaxa|dada2|sintax)",
      remove = FALSE
    ) %>%
    tidyr::unnest("value") %>%
    dplyr::filter(
      !is.na(rank),
      !is.na(label),
      confidence >= 0.5,
      !startsWith(taxon, "unidentified"),
      !endsWith(taxon, "incertae_sedis"),
      reference != "warcup" | rank != "kingdom"
    ) %>%
    dplyr::mutate_at("rank", rank_factor) %>%
    dplyr::group_by(label, rank) %>%
    dplyr::filter((!"ITS" %in% region) | region != "short") %>%
    dplyr::mutate(
      n_tot = dplyr::n(),
      n_diff = dplyr::n_distinct(taxon, na.rm = TRUE),
      n_method = dplyr::n_distinct(method, na.rm = TRUE),
      n_reference = dplyr::n_distinct(reference, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      dplyr::select(allseqs, label = "hash", n_reads = "nread") %>%
        dplyr::group_by(label) %>%
        dplyr::summarize(n_reads = sum(n_reads)) %>%
        dplyr::ungroup(),
      by = "label"
    )
}

#### taxon_labels ####
# make labels summarizing the taxonomy of each sequence
make_taxon_labels <- function(t) {
    dplyr::group_by(t, label, rank, n_reads) %>%
    dplyr::summarize(
      taxon =
        table(taxon) %>%
        paste0(names(.), collapse = "/") %>%
        gsub(pattern = "(.+/.+)", replacement = "<\\1>") %>%
        gsub(pattern = "(mycota|mycetes|ales|aceae)", replacement = "") %>%
        gsub(pattern = "incertae_sedis", replacement = "i_s") %>%
        gsub(pattern = "Fungi\\b", replacement = "F") %>%
        gsub(pattern = "Basidio\\b", replacement = "B") %>%
        gsub(pattern = "Asco\\b", replacement = "A") %>%
        gsub(pattern = "Chytridio\\b", replacement = "Chy") %>%
        gsub(pattern = "Zygo\\b", replacement = "Z")
    ) %>%
    dplyr::group_by(label, n_reads) %>%
    dplyr::arrange(rank) %>%
    dplyr::summarize(tip_label = paste(label[1],
                                       format(n_reads[1], width = 5),
                                       paste0(taxon, collapse = "-")))
}

#### relabel_tree ####
# replaces tree tip labels from old with labels from new
relabel_tree <- function(tree, old, new, chimeras = character(0)) {
  tree <- ape::drop.tip(tree, intersect(chimeras, tree$tip.label))
  tree$tip.label <-
    plyr::mapvalues(tree$tip.label, old, paste0('"', new, '"'), warn_missing = FALSE)
  tree
}
