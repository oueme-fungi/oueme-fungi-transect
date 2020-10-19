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
      taxon != "NA",
      confidence >= 0.5,
      !startsWith(taxon, "unidentified"),
      !endsWith(taxon, "incertae_sedis"),
      reference != "warcup" | rank != "kingdom"
    ) %>%
    dplyr::mutate_at("rank", rank_factor) %>%
    dplyr::left_join(
      dplyr::select(allseqs, label = "hash", n_reads = "nread") %>%
        dplyr::group_by(label) %>%
        dplyr::summarize(n_reads = sum(n_reads)) %>%
        dplyr::ungroup(),
      by = "label"
    )
}

#### Find ASVs with consistent kingdom-level assignments
# i.e., at least min_n assignments with greater than min_confidence
# confidence, and also no conflicting assignments at any confidence level
# they should also be present in the tree, and not one of the labels given
# by "ignore"
identify_taxa <- function(taxon_table, tree, rank, confidence, n,
                          ignore = character()) {
  taxon_table %>%
    # look at kingdom assignments from the long amplicon dataset
    # don't include tulasnella
    dplyr::filter(
      rank == !!rank,
      label %in% tree$tip.label,
      !label %in% ignore
    ) %>%
    dplyr::group_by(label) %>%
    dplyr::filter(
      # at least 3 "confident" assignments
      sum(confidence >= !!confidence) >= !!n,
      # no conflicting assignments
      dplyr::n_distinct(taxon) == 1
    ) %>%
    dplyr::select(label, taxon) %>%
    dplyr::ungroup() %>%
    unique()
}

# warcup doesn't jave kingdom annotations, because it only includes fungi.
# however, if there is an annotation, then call it a fungus.
add_warcup_kingdom <- function(taxon_table) {
  taxon_table %>%
    dplyr::filter(
      reference == "warcup",
      rank == "phylum"
    ) %>%
    dplyr::mutate(
      rank = "kingdom",
      taxon = "Fungi"
    ) %>%
    dplyr::bind_rows(
      taxon_table,
      .
    )
}

## identifies a clade which is probably bikonta based on the confident
# kingdom assignments, and returns all the tips in it
extract_bikonta <- function(kingdoms, tree) {
  bikonta_list <- c("Alveolata", "Rhizaria", "Stramenopila", "Viridiplantae")
  bikonta <- dplyr::filter(kingdoms, taxon %in% bikonta_list)
  nonbikonta <- dplyr::filter(kingdoms, !taxon %in% bikonta_list)
  biconta_mrca <- ape::getMRCA(tree, bikonta$label)
  # error if the root is currently inside bikonta
  stopifnot(biconta_mrca != phangorn::getRoot(tree))
  bikonta <- tree$tip.label[unlist(phangorn::Descendants(tree, biconta_mrca))]
  # error if bikonta is not monophyletic with respect to the other
  # confident assignments
  stopifnot(!any(nonbikonta$label %in% bikonta))
  bikonta
}

# combine a set of taxon identifications and a set of sample read counts
select_taxon_reads <- function(taxa, reads, ...,
                               method = first(taxa$assigned$method)) {
  taxa$assigned%>%
    select("method", "label", "rank", "taxon", "region") %>%
    unique() %>%
    pivot_wider(names_from = "rank", values_from = "taxon") %>%
    right_join(filter(reads, ...), by = "label") %>%
    mutate_at("kingdom", na_if, "NA") %>%
    mutate(kingdom = ifelse(endsWith(phylum, "mycota"), "Fungi", kingdom),
           method = !!method)
}

# this function would be helpful in FUNGuildR
widen_taxonomy <- function(taxa) {
    dplyr::left_join(
      tidyr::spread(taxa, key = rank, value = taxon),
      dplyr::group_by(taxa, label) %>%
        dplyr::filter(!is.na(taxon)) %>%
        dplyr::summarize(Taxonomy = paste(taxon, collapse = ";")),
      by = "label"
    )
}
