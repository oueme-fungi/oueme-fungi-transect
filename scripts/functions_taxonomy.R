# functions to deal with taxonomic annotations
# many of the functions originally here have been moved the phylotax package.
# author Brendan Furneaux

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
most_common <- function(x) {
  ux <- unique(x)
  ux[which.max(table(match(x, ux)))]
}

# combine a set of taxon identifications and a set of sample read counts
select_taxon_reads <- function(taxa, reads, ...,
                               method = most_common(taxa$assigned$method),
                               region = most_common(taxa$assigned$region),
                               reference = most_common(taxa$assigned$reference)) {
  taxa$assigned %>%
    select("method", "label", "rank", "taxon", "region") %>%
    unique() %>%
    pivot_wider(names_from = "rank", values_from = "taxon") %>%
    right_join(filter(reads, ...), by = "label") %>%
    mutate_at("kingdom", na_if, "NA") %>%
    mutate(kingdom = ifelse(endsWith(phylum, "mycota"), "Fungi", kingdom),
           method = !!method, region = !!region, reference = !!reference)
}

# fraction of reads for each ASV in each sequencing run
make_taxon_reads_proto <- function(proto_physeq, datasets) {
  # Get the number of reads per sequencing run for each ASV
  phyloseq::merge_samples(proto_physeq, group = "seq_run") %>%
    phyloseq::otu_table() %>%
    t() %>%
    as("matrix") %>%
    as_tibble(rownames = "label") %>%
    # Normalize reads to a fraction of the sequencing run.
    mutate_if(is.numeric, ~./sum(.)) %>%
    # Pivot so that each ASV/sequencing run combination is its own line, and remove 0's
    pivot_longer(-1, names_to = "seq_run", values_to = "reads") %>%
    filter(reads > 0) %>%
    # Add the amplicon and technology used for each sequencing run
    left_join(select(datasets, seq_run, amplicon, tech), by = "seq_run") %>%
    # put the amplicon in title case for display
    mutate_at("amplicon", stringr::str_to_title)
}

compile_primary_taxon_reads <- function(proto, taxon_table_primary) {
  # Make all valid combinations of reference, region, and method
  expand_grid(
    proto,
    tibble(
      reference = c("unite", "unite", "warcup", "warcup", "rdp_train"),
      region = c("ITS", "short", "ITS", "short", "LSU")
    ),
    method = c("dada2", "idtaxa", "sintax")
  ) %>%
    filter(
      ifelse(amplicon == "Short", region == "short", TRUE),
      ifelse(amplicon == "Long", region %in% c("ITS", "LSU"), TRUE)
    ) %>%
    # Add the identifications from the taxon table.
    # This will leave NAs where a method did not make an identification at a rank,
    # Including full NA rows where no identification was made at all.
    left_join(
      taxon_table_primary %>%
        mutate_at("taxon", na_if, "NA") %>%
        select(label, method, reference, rank, taxon, region) %>%
        pivot_wider(
          names_from = "rank",
          values_from = "taxon"
        ) %>%
        mutate(
          family = ifelse(is.na(genus), family, coalesce(family, "incertae_sedis")),
          order = ifelse(is.na(family), order, coalesce(order, "incertae_sedis")),
          class = ifelse(is.na(order), class, coalesce(class, "incertae_sedis")),
          phylum = ifelse(!is.na(class) & class == "Leotiomycetes", "Ascomycota", phylum)
        ),
      by = c("label", "reference", "region", "method")
    )
}

format_taxon_reads <- function(taxon_reads, funguild_db) {
  mutate(
    taxon_reads,
    family = ifelse(is.na(genus), family, coalesce(family, "incertae_sedis")),
    order = ifelse(is.na(family), order, coalesce(order, "incertae_sedis")),
    class = ifelse(is.na(order), class, coalesce(class, "incertae_sedis")),
    # Not sure what happened here...
    phylum = ifelse(!is.na(class) & class == "Leotiomycetes", "Ascomycota", phylum)
  ) %>%
    mutate_at(
      "method",
      factor,
      levels = c("PHYLOTAX", "combined", "Consensus", "dada2", "sintax", "idtaxa"),
      labels = c("PHYLOTAX", "PHYLOTAX", "Consensus", "RDPC", "SINTAX", "IDTAXA")
    ) %>%
    mutate_at("reference", replace_na, "All") %>%
    mutate_at(
      "reference",
      factor,
      levels = c("All", "unite", "warcup", "rdp_train"),
      labels = c("All", "Unite", "Warcup", "RDP")
    ) %>%
    mutate_at("reference", replace_na, "All") %>%
    mutate_at("tech", factor, levels = c("PacBio", "Illumina", "Ion Torrent")) %>%
    rename(Algorithm = method) %>%
    mutate(
      Taxonomy = paste(kingdom, phylum, class, order, family, genus, sep = ";") %>%
        str_replace_all(";NA", "")
    ) %>%
    left_join(
      FUNGuildR::funguild_assign(
        unique(select(., Taxonomy)),
        funguild_db
      ) %>%
        # Recent versions of FUNGuild (2021-2-3) have some taxa with two entries
        # this combines them.
        group_by(Taxonomy) %>%
        summarize(
          guild = str_c(unlist(strsplit(guild, "-")), collapse = "-"),
          confidenceRanking = unique(confidenceRanking)
        ) %>%
        ungroup() %>%
        select(Taxonomy, guild, confidenceRanking),
      by = "Taxonomy"
    ) %>%
    mutate(
      ECM = grepl("Ectomycorrhizal", guild) %>%
        ifelse(paste(confidenceRanking, "ECM"),
               ifelse(is.na(guild), NA, "non-ECM")) %>%
        ifelse(is.na(family), NA, .)) %>%
    mutate_at(
      "ECM",
      factor,
      levels = c(NA_character_, "non-ECM", "Possible ECM",
                 "Probable ECM", "Highly Probable ECM"),
      exclude = NULL,
      ordered = TRUE
    ) %>%
    select(-Taxonomy, -guild, -confidenceRanking)
}

compile_taxon_chart <- function(taxon_reads) {
  pivot_longer(
    taxon_reads,
    names_to = "rank",
    cols = c("kingdom", "phylum", "class", "order", "family", "genus"),
    values_to = "taxon"
  ) %>%
    filter(Algorithm != "Consensus" | amplicon != "Short" | region == "short") %>%
    select(label, amplicon, tech, Algorithm, reference, rank, taxon, reads) %>%
    unique() %>%
    group_by(amplicon, tech, Algorithm, reference, rank, ID = !is.na(taxon)) %>%
    summarize(reads = sum(reads), ASVs = n()) %>%
    group_by(amplicon, tech, Algorithm, reference, rank) %>%
    mutate(ASVs = ASVs / sum(ASVs)) %>%
    ungroup() %>%
    filter(ID) %>%
    mutate_at("rank", factor,
              levels = c("kingdom", "phylum", "class", "order", "family", "genus")
    ) %>%
    rename(Reads = reads) %>%
    mutate_at("reference", as.character) %>%
    mutate_at("reference", factor, levels = c("All", "Unite", "Warcup", "RDP")) %>%
    pivot_longer(cols = c("Reads", "ASVs"), names_to = "type", values_to = "frac")
}

draw_taxon_chart <- function(taxon_chart) {
  dplyr::arrange(taxon_chart, rank) %>%
    ggplot(aes(x = Algorithm, y = frac, ymax = frac, group = rank, fill = rank)) +
    ggnomics::facet_nested(
      tech + amplicon ~ type + reference,
      resect = unit(1, "mm"),
      nest_line = TRUE,
      bleed = FALSE,
      scales = "free_x",
      space = "free_x"
    ) +
    # geom_ribbon(ymin = 0, alpha = 0.5) +
    # geom_line(aes(color = Algorithm)) +
    geom_col(position = "identity") +
    scale_fill_brewer(type = "qual", palette = 2, direction = -1, guide = guide_legend(nrow = 1)) +
    ylab("Fraction of reads assigned") +
    xlab(NULL) +
    theme(axis.text.x = element_text(vjust = 0.5, angle = 90),
          legend.direction = "horizontal",
          legend.position = "bottom",
          strip.background = element_blank())
}



# Plot distribution between different taxa
taxon_plot <- function(.data, rank, ..., y = reads, x = Algorithm,
                       facets = vars(tech, amplicon, reference), cutoff = NULL,
                       datasets) {
  rank <- enquo(rank)
  y <- enquo(y)
  x <- enquo(x)
  facets <- enquo(facets)
  facets <- enquo(facets)
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
  .data <- .data %>%
    group_by(!!x) %>%
    group_by_at(rlang::eval_tidy(facets), .add = TRUE) %>%
    mutate(ASVs = n_distinct(label)) %>%
    ungroup() %>%
    filter(...) %>%
    arrange_at(ranks) %>%
    mutate_at(
      ranks,
      ~ factor(
        .,
        levels = c(NA, "other", discard(unique(.), is.na)),
        exclude = NULL
      )
    ) %>%
    group_by(!!x, !!rank) %>%
    group_by_at(rlang::eval_tidy(facets), .add = TRUE) %>%
    summarize(reads = sum(reads), ASVs = n_distinct(label)/max(ASVs)) %>%
    ungroup()

  if (!is.null(cutoff)) {
    prelevels <- levels(pull(.data, !!rank))
    .data <- group_by(.data, !!rank) %>%
      group_map(

        ~ if (all(pull(.x, !!y) < cutoff)) {
          mutate(.x, !!rank := factor("other", levels = levels(!!rank)), exclude = NULL)
        } else {
          .x
        },
        keep = TRUE
      ) %>%
      bind_rows() %>%
      group_by(!!x, !!rank) %>%
      group_by_at(rlang::eval_tidy(facets), .add = TRUE) %>%
      summarize(reads = sum(reads), ASVs = sum(ASVs)) %>%
      ungroup() %>%
      mutate(!!rank := factor(!!rank, levels = prelevels, exclude = NULL))
  }

  rank_label <- as_label(rank)

  .data <- mutate(.data, !!rank := fct_drop(!!rank))
  vals <- levels(pull(.data, !!rank))
  # if ("other" %in% vals) vals <- c("other", vals) %>% magrittr::extract(!duplicated(.))
  # if (any(is.na(vals))) vals <- c(NA, vals) %>% magrittr::extract(!duplicated(.))


  if (rank_label == str_to_lower(rank_label)) rank_label <- str_to_title(rank_label)
  y_label <- as_label(y)
  if (y_label == str_to_lower(y_label)) y_label <- str_to_title(y_label)
  .data <- .data %>%
    mutate(tech = fct_relabel(tech, ~ paste(., plyr::mapvalues(., datasets$tech, datasets$machine, FALSE))
    ))
  ggplot(.data, aes(x = !!x, y = !!y, fill = !!rank)) +
    geom_bar(position = "stack", stat = "identity", color = "white", size = 0.2) +
    ggnomics::facet_nested(
      cols = rlang::eval_tidy(facets),
      scales = "free_x",
      space = "free_x",
      nest_line = TRUE,
      bleed = FALSE,
      resect = unit(1, "mm")
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.background = element_blank(),
          panel.spacing = unit(3, "pt")) +
    scale_fill_discrete(
      breaks = vals,
      labels = replace_na(as.character(vals), "unidentified"),
      name = rank_label
    ) +
    ylab(paste("Fraction of", y_label))

}

draw_fungal_classes <- function(taxon_reads, datasets) {
  filter(
    taxon_reads,
    kingdom == "Fungi",
    !(Algorithm == "Consensus" & region == "ITS"),
    Algorithm == "PHYLOTAX"
  ) %>%
    group_by(tech, amplicon) %>%
    mutate(reads = reads/sum(reads), ASVs = 1/n()) %>%
    pivot_longer(cols = c("reads", "ASVs"), names_to = "type", values_to = "reads") %>%
    taxon_plot(
      rank = class,
      datasets = datasets,
      cutoff = 0.02,
      y = reads,
      x = type,
      facets = vars(tech, amplicon)
    ) +
    ylab("Fraction of Community") +
    xlab(NULL)
}

draw_short_kingdoms <- function(seq_table, taxon_reads, datasets) {
    colSums(seq_table) %>%
    enframe(name = "ITS2", value = "nreads") %>%
    mutate(length = nchar(ITS2)) %>%
    filter(length <= 140) %>%
    mutate_at("ITS2", ~tzara::seqhash(chartr("T", "U", .))) %>%
    semi_join(taxon_reads, ., by = c("label" = "ITS2")) %>%
    group_by(tech, amplicon, Algorithm, reference) %>%
    mutate(reads = reads/sum(reads)) %>%
    ungroup() %>%
    # group_by(ITS2) %>%
    # filter(!"Short" %in% amplicon) %>%
    taxon_plot(kingdom, y = reads, datasets = datasets)
}

draw_ecm_assignment <- function(taxon_reads, datasets) {
    filter(taxon_reads,
           ifelse(Algorithm == "Consensus", region == "hybrid", TRUE)
    ) %>%
    mutate_at("reference", replace_na, "All") %>%
    taxon_plot(ECM, kingdom == "Fungi", datasets = datasets) +
    theme(legend.position = "bottom", legend.title = element_blank())
}
