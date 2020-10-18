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

# calculate last common ancestor consensus
# i.e., if multipe assignments disagree, take the last common ancestor
# that does agree.
# this is the same thing as strict consensus at each rank.
lca_consensus <- function(
  taxa, ranks = NULL,
  method = if (utils::hasName(taxa, "method")) "LCA" else NULL
) {
  method <- phylotax:::check_method(taxa, method)
  taxa <- dplyr::select(taxa, "label", dplyr::one_of(names(method)), "rank",
                        "taxon")
  taxa <- phylotax:::check_ranks(taxa, ranks)
  taxa <- phylotax:::count_assignments(taxa)
  tip_taxa <- dplyr::group_by_at(taxa, c("label", names(method))) %>%
    dplyr::arrange(.data$rank) %>%
    dplyr::filter(dplyr::cumall(.data$n_diff == 1)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::one_of(c(names(method), "n_diff", "n_tot"))) %>%
    unique()
  for (n in names(method)) {
    tip_taxa[[n]] <- unname(method[n])
  }
  taxa <- dplyr::select(taxa, -"n_diff", -"n_tot")
  structure(
    list(
      node_taxa = NULL,
      tree = NULL,
      tip_taxa = tip_taxa,
      retained = dplyr::semi_join(taxa, tip_taxa, by = c("label", "rank")),
      rejected = dplyr::anti_join(taxa, tip_taxa, by = c("label", "rank")),
      missing = dplyr::filter(taxa, FALSE)
    ),
    class = "phylotax"
  )
}

combotax <- function(phylotax, lca = NULL, method = if (utils::hasName(phylotax$tip_taxa, "method")) "PHYLOTAX" else NULL) {
  method <- phylotax:::check_method(phylotax$tip_taxa, method)
  assertthat::assert_that(methods::is(phylotax, "phylotax"))
  if (is.null(lca)) {
    lca <- lca_consensus(phylotax$missing, method = method)
  } else {
    for (n in names(method))
      lca$tip_taxa[[n]] <- unname(method[n])
  }
  phylotax$tip_taxa <- dplyr::bind_rows(
    phylotax$tip_taxa,
    dplyr::anti_join(lca$tip_taxa, phylotax$tip_taxa, by = "label")
  )
  phylotax$rejected <- unique(dplyr::bind_rows(phylotax$rejected, lca$rejected))
  phylotax$retained <- unique(dplyr::bind_rows(phylotax$retainsd, lca$retained))
  phylotax$missing <- purrr::reduce(
    list(phylotax$missing, phylotax$rejected, phylotax$retained),
    dplyr::anti_join, by = names(method)
  )
  phylotax
}

# combine a set of taxon identifications and a set of sample read counts
select_taxon_reads <- function(taxa, reads, ..., method = first(taxa$tip_taxa$method)) {
  taxa$tip_taxa %>%
    select("method", "label", "rank", "taxon", "region") %>%
    unique() %>%
    pivot_wider(names_from = "rank", values_from = "taxon") %>%
    right_join(filter(reads, ...), by = "label") %>%
    mutate_at("kingdom", na_if, "NA") %>%
    mutate(kingdom = ifelse(endsWith(phylum, "mycota"), "Fungi", kingdom),
           method = !!method)
}
