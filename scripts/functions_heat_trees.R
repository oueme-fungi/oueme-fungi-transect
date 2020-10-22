# functions for making heat trees
# author Brendan Furneaux

generate_taxdata <- function(taxon_reads, asv_table) {
  ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
  out <- taxon_reads %>%
    mutate(domain = "Root", ASVs = 1) %>%
    select(-region, -seq_run, -label) %>%
    mutate_at(
      ranks,
      replace_na,
      "?"
    ) %>%
    group_by_at(vars(-one_of(ranks, "reads", "ASVs"))) %>%
    mutate(ASVs = ASVs / sum(ASVs))
  for (i in seq_along(ranks)) {
    out <- out %>%
      group_by_at(rev(ranks)[i:length(ranks)]) %>%
      group_map(
        function(x, y, i) {
          d <- group_by_at(x, vars(-one_of(ranks, "reads", "ASVs"))) %>%
            summarize_at(c("reads", "ASVs"), sum)
          d2 <- group_by_at(d, "amplicon") %>% summarize_at(c("reads", "ASVs"), mean)
          reads <- select(d2, "amplicon", "reads") %>% deframe()
          ASVs <- select(d2, "amplicon", "ASVs") %>% deframe()
          if (max(c((d$reads), (d$ASVs))) < 0.01) {
            if (
              (isTRUE(abs(log10(reads["Long"]) - log10(reads["Short"])) < 1) &
               isTRUE(abs(log10(ASVs["Long"]) - log10(ASVs["Short"])) < 1)) |
              max(c((d$reads), (d$ASVs))) < 0.01
            ) {
              mutate_at(x, rev(ranks)[i], ~"*")
            } else if (i > 1) {
              mutate_at(x, rev(ranks)[i - 1], ~"*") %>%
                mutate_at(rev(ranks)[i], ~paste0("(",., ")"))
            } else {
              mutate_at(x, rev(ranks)[i], ~paste0("(",., ")"))
            }
          } else {
            x
          }
        },
        i = i,
        keep = TRUE
      ) %>%
      bind_rows()
  }
  out <- group_by_at(out, vars(-one_of("reads", "ASVs", "ECM"))) %>%
    summarize(reads = sum(reads), ASVs = sum(ASVs), ECM = min(ECM)) %>%
    pivot_wider(
      names_from = c("tech", "amplicon", "reference", "Algorithm"),
      values_from = c("reads", "ASVs"),
      values_fill = list(reads = 0, ASVs = 0),
      names_repair = "universal"
    ) %>%
    taxa::parse_tax_data(
      class_cols = ranks
    ) %>%
    taxa::filter_taxa(
      taxon_names != "incertae_sedis",
      drop_obs = TRUE,
      reassign_obs = FALSE,
      reassign_taxa = TRUE
    )

  read_cols <- keep(names(out$data$tax_data), startsWith, "reads_")
  asv_cols <- keep(names(out$data$tax_data), startsWith, "ASVs_")

  out$data$asv_table <- metacoder::calc_obs_props(out, "tax_data", cols = asv_cols)
  out$data$read_table <- metacoder::calc_obs_props(out, "tax_data", cols = read_cols)
  out$data$tax_asv <- metacoder::calc_taxon_abund(out, "asv_table", cols = asv_cols)
  out$data$tax_read <- metacoder::calc_taxon_abund(out, "read_table", cols = read_cols)

  child_of_unknown <-  unlist(
    out$supertaxa(
      recursive = FALSE,
      value = "taxon_names",
      na = TRUE
    )
  ) %in% c("?", "*")
  unknown_taxa <- out$taxon_names() %in% c("?", "*", "incertae_sedis")
  out <- out$filter_taxa(
    !(child_of_unknown & unknown_taxa) | taxon_names == "Root",
    reassign_obs = FALSE
  )
  dangling_taxa <- map_lgl(
    out$subtaxa(recursive = FALSE)[unlist(out$supertaxa(recursive = FALSE, na = TRUE))],
    ~all(out$taxon_names()[.] %in% c("?", "*") | is.na(out$taxon_names()[.]))
  )
  out <- out$filter_taxa(
    !dangling_taxa | taxon_names == "Root",
    reassign_obs = FALSE
  )

  retaxa_table <- tibble::tibble(
    classification = out$classifications(),
    id = out$taxon_ids()
  ) %>%
    group_by(classification) %>%
    mutate(newid = first(id)) %>%
    filter(id != newid)
  for (n in names(out$data)) {
    if ("taxon_id" %in% names(out$data[[n]]))
      out$data[[n]]$taxon_id <- plyr::mapvalues(
        out$data[[n]]$taxon_id,
        retaxa_table$id,
        retaxa_table$newid,
        FALSE
      )
  }
  for (n in c("asv_table", "read_table", "tax_asv", "tax_read")) {
    out$data[[n]] <- out$data[[n]] %>%
      group_by(taxon_id) %>%
      summarize_all(sum)
  }
  out <- out$filter_taxa(!taxon_ids %in% retaxa_table$id, drop_obs = TRUE,
                         reassign_obs = FALSE, reassign_taxa = FALSE)

  out$data$tax_read$mean = rowMeans(out$data$tax_read[,-1])
  out$data$tax_read$max = do.call(pmax, out$data$tax_read[,-1])
  out$data$tax_asv$mean = rowMeans(out$data$tax_asv[,-1])
  out$data$tax_asv$mean = do.call(pmax, out$data$tax_asv[,-1])

  out$data$diff_read <- metacoder::compare_groups(
    out,
    data = "tax_read",
    cols = c(read_cols, "mean"),
    groups = c(read_cols, "mean"),
    combinations = lapply(read_cols, c, "mean")
  )
  out$data$diff_asv <- metacoder::compare_groups(
    out,
    data = "tax_read",
    cols = c(read_cols, "mean"),
    groups = c(read_cols, "mean"),
    combinations = lapply(read_cols, c, "mean")
  )
  out
}

generate_ecm_taxdata <- function(taxon_reads, asv_table) {
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
  out <- taxon_reads %>%
    filter(!is.na(as.character(ECM)), ECM != "non-ECM") %>%
    mutate(ASVs = 1) %>%
    select(-region, -seq_run, -label) %>%
    mutate_at(
      ranks,
      replace_na,
      "?"
    ) %>%
    group_by_at(vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
    mutate(ASVs = ASVs / sum(ASVs))
  for (i in seq_along(ranks)) {
    out <- out %>%
      group_by_at(rev(ranks)[i:length(ranks)]) %>%
      group_map(
        function(x, y, i) {
          d <- group_by_at(x, vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
            summarize_at(c("reads", "ASVs"), sum)
          d2 <- group_by_at(d, "amplicon") %>% summarize_at(c("reads", "ASVs"), mean)
          reads <- select(d2, "amplicon", "reads") %>% deframe()
          ASVs <- select(d2, "amplicon", "ASVs") %>% deframe()
          if (max(c((d$reads), (d$ASVs))) < 0.01) {
            if (
              isTRUE(abs(log10(reads["Long"]) - log10(reads["Short"])) < 1) &
              isTRUE(abs(log10(ASVs["Long"]) - log10(ASVs["Short"])) < 1)
            ) {
              mutate_at(x, rev(ranks)[i], ~"*")
            } else if (i > 1) {
              mutate_at(x, rev(ranks)[i - 1], ~"*") %>%
                mutate_at(rev(ranks)[i], ~paste0("(",., ")"))
            } else {
              mutate_at(x, rev(ranks)[i], ~paste0("(",., ")"))
            }
          } else {
            x
          }
        },
        i = i,
        keep = TRUE
      ) %>%
      bind_rows()
  }
  out <- group_by_at(out, vars(-one_of("reads", "ASVs", "ECM"))) %>%
    summarize(reads = sum(reads), ASVs = sum(ASVs), ECM = min(ECM)) %>%
    pivot_wider(
      names_from = c("tech", "amplicon", "reference", "Algorithm"),
      values_from = c("reads", "ASVs"),
      values_fill = list(reads = 0, ASVs = 0),
      names_repair = "universal"
    ) %>%
    taxa::parse_tax_data(
      class_cols = ranks
    ) %>%
    taxa::filter_taxa(
      taxon_names != "incertae_sedis",
      drop_obs = TRUE,
      reassign_obs = FALSE,
      reassign_taxa = TRUE
    )

  read_cols <- keep(names(out$data$tax_data), startsWith, "reads_")
  asv_cols <- keep(names(out$data$tax_data), startsWith, "ASVs_")

  out$data$asv_table <- metacoder::calc_obs_props(out, "tax_data", cols = asv_cols)
  out$data$read_table <- metacoder::calc_obs_props(out, "tax_data", cols = read_cols)
  out$data$tax_asv <- metacoder::calc_taxon_abund(out, "asv_table", cols = asv_cols)
  out$data$tax_read <- metacoder::calc_taxon_abund(out, "read_table", cols = read_cols)

  child_of_unknown <-  unlist(
    out$supertaxa(
      recursive = FALSE,
      value = "taxon_names",
      na = TRUE
    )
  ) %in% c("?", "*")
  out <- out$filter_taxa(
    !child_of_unknown | taxon_names == "Fungi",
    reassign_obs = FALSE
  )
  dangling_taxa <- map_lgl(
    out$subtaxa(recursive = FALSE)[unlist(out$supertaxa(recursive = FALSE, na = TRUE))],
    ~all(out$taxon_names()[.] %in% c("?", "*") | is.na(out$taxon_names()[.]))
  )
  out <- out$filter_taxa(
    !dangling_taxa | taxon_names == "Fungi",
    reassign_obs = FALSE
  )

  retaxa_table <- tibble::tibble(
    classification = out$classifications(),
    id = out$taxon_ids()
  ) %>%
    group_by(classification) %>%
    mutate(newid = first(id)) %>%
    filter(id != newid)
  for (n in names(out$data)) {
    if ("taxon_id" %in% names(out$data[[n]]))
      out$data[[n]]$taxon_id <- plyr::mapvalues(
        out$data[[n]]$taxon_id,
        retaxa_table$id,
        retaxa_table$newid,
        FALSE
      )
  }
  for (n in c("asv_table", "read_table", "tax_asv", "tax_read")) {
    out$data[[n]] <- out$data[[n]] %>%
      group_by(taxon_id) %>%
      summarize_all(sum)
  }
  out <- out$filter_taxa(!taxon_ids %in% retaxa_table$id, drop_obs = TRUE,
                         reassign_obs = FALSE, reassign_taxa = FALSE)

  out$data$tax_read$mean = rowMeans(out$data$tax_read[,-1])
  out$data$tax_read$max = do.call(pmax, out$data$tax_read[,-1])
  out$data$tax_asv$mean = rowMeans(out$data$tax_asv[,-1])
  out$data$tax_asv$mean = do.call(pmax, out$data$tax_asv[,-1])

  out$data$diff_read <- metacoder::compare_groups(
    out,
    data = "tax_read",
    cols = c(read_cols, "mean"),
    groups = c(read_cols, "mean"),
    combinations = lapply(read_cols, c, "mean")
  )
  out$data$diff_asv <- metacoder::compare_groups(
    out,
    data = "tax_read",
    cols = c(read_cols, "mean"),
    groups = c(read_cols, "mean"),
    combinations = lapply(read_cols, c, "mean")
  )
  out
}

draw_full_heattree <- function(taxdata, outfile) {
  set.seed(4)
  theme_update(panel.border = element_blank())
  taxdata %>%
    metacoder::heat_tree(
      layout = "davidson-harel",
      # initial_layout = "reingold-tilford",
      node_size = rowMeans(.$data$tax_asv[,-1]),
      node_color = rowMeans(.$data$tax_asv[,-1]),
      node_label = replace_na(taxon_names, "unknown"),
      node_size_range = c(0.002, .035),
      node_size_axis_label = "ASV count",
      node_color_axis_label = "ASV count",
      node_label_size_range = c(0.01, 0.03),
      edge_size = rowMeans(.$data$tax_read[,-1]),
      edge_size_trans = "linear",
      edge_size_axis_label = "Read abundance",
      edge_size_range = c(0.001, 0.03),
      edge_color = rowMeans(.$data$tax_read[,-1]),
      edge_color_axis_label = "Read abundance",
      aspect_ratio = 9/5,
      make_node_legend = FALSE,
      make_edge_legend = FALSE,
      output_file = outfile
    )

}

log_read_ratio <- function(abund1, abund2) {
  list(read_ratio = log10(mean(abund1) / mean(abund2)) %>% ifelse(is.nan(.), 0, .))
}

log_ASV_ratio <- function(abund1, abund2) {
  list(ASV_ratio = log10(mean(abund1) / mean(abund2)) %>% ifelse(is.nan(.), 0, .))
}

draw_amplicons_heattree <- function(taxdata, outfile) {
  taxdata2 <- taxdata
  val_cols = names(taxdata2$data$tax_read) %>%
    keep(str_detect, "reads_.+_.+_.+_Consensus")

  taxdata2$filter_obs(
    "tax_read",
    rowSums(select_at(taxdata2$data$tax_read, val_cols)) > 0,
    drop_taxa = TRUE,
    drop_obs = TRUE,
    supertaxa = TRUE,
    reassign_obs = FALSE) %>%
    invisible()

  taxdata2$data$diff_read <- metacoder::compare_groups(
    taxdata2,
    data = "tax_read",
    cols = val_cols,
    groups = val_cols %>%
      str_replace("reads_.+_(.+)_.+_.+", "\\1"),
    func = log_read_ratio
  )
  taxdata2$data$diff_asv <- metacoder::compare_groups(
    taxdata2,
    data = "tax_asv",
    cols = names(taxdata2$data$tax_asv) %>% keep(startsWith, "ASVs_"),
    groups = names(taxdata2$data$tax_asv) %>%
      keep(startsWith, "ASVs_") %>%
      str_replace("ASVs_.+_(.+)_.+_.+", "\\1"),
    func = log_ASV_ratio
  )

  set.seed(4)
  theme_update(panel.border = element_blank())
  taxdata2 %>%
    metacoder::heat_tree(
      layout = "davidson-harel",
      # initial_layout = "reingold-tilford",

      node_size = rowMeans(.$data$tax_asv[,-1]),
      node_size_range = c(0.002, .035),
      node_size_axis_label = "ASV richness",

      node_color = ASV_ratio,
      node_color_range = metacoder::diverging_palette(),
      node_color_trans = "linear",
      node_color_axis_label = "Log richness ratio",
      node_color_interval = c(-1.5, 1.5),

      edge_size = rowMeans(.$data$tax_read[,-1]),
      edge_size_range = c(0.001, 0.03),
      edge_size_trans = "linear",
      edge_size_axis_label = "Read abundance",

      edge_color = read_ratio,
      edge_color_range = metacoder::diverging_palette(),
      edge_color_trans = "linear",
      edge_color_axis_label = "Log abundance ratio",
      edge_color_interval = c(-3, 3),

      node_label = ifelse(
        pmax(abs(.$data$diff_read$read_ratio) > 0.5,
             abs(.$data$diff_asv$ASV_ratio) > 0.25),
        taxon_names,
        ""
      ),
      # aspect_ratio = 3/2,
      node_label_size_range = c(0.010, 0.03),
      output_file = outfile
      # title = unique(paste(.$data$diff_read$treatment_1, "vs.", .$data$diff_read$treatment_2))
    )
}

draw_ecm_heattree <- function(taxdata_ECM, outfile) {
  taxdata2 <- taxdata_ECM

  val_cols <- names(taxdata2$data$tax_read) %>%
    keep(str_detect, "reads_.+_.+_.+_Consensus")

  taxdata2$filter_obs(
    "tax_read",
    rowSums(select_at(taxdata2$data$tax_read, val_cols)) > 0,
    drop_taxa = TRUE,
    drop_obs = TRUE,
    supertaxa = TRUE,
    reassign_obs = FALSE) %>%
    invisible()

  taxdata2$data$diff_read <- metacoder::compare_groups(
    taxdata2,
    data = "tax_read",
    cols = val_cols,
    groups = val_cols %>%
      str_replace("reads_.+_(.+)_.+_Consensus", "\\1"),
    func = log_read_ratio
  )
  taxdata2$data$diff_asv <- metacoder::compare_groups(
    taxdata2,
    data = "tax_asv",
    cols = names(taxdata2$data$tax_asv) %>% keep(startsWith, "ASVs_"),
    groups = names(taxdata2$data$tax_asv) %>%
      keep(startsWith, "ASVs_") %>%
      str_replace("ASVs_.+_(.+)_.+_.+", "\\1"),
    func = log_ASV_ratio
  )

  taxdata2 <- taxdata2$filter_taxa(!is.nan(ASV_ratio), !is.nan(read_ratio))

  set.seed(3)
  theme_update(panel.border = element_blank())
  taxdata2 %>%
    metacoder::heat_tree(
      layout = "davidson-harel",
      initial_layout = "reingold-tilford",

      node_size = rowMeans(.$data$tax_asv[,-1]),
      node_size_range = c(0.002, .035),
      node_size_axis_label = "ASV richness",

      node_color = ASV_ratio,
      node_color_range = metacoder::diverging_palette(),
      node_color_trans = "linear",
      node_color_axis_label = "Log richness ratio",
      node_color_interval = c(-0.9, 0.9),

      edge_size = rowMeans(.$data$tax_read[,-1]),
      edge_size_range = c(0.001, 0.03),
      edge_size_trans = "linear",
      edge_size_axis_label = "Read abundance",

      edge_color = read_ratio,
      edge_color_range = metacoder::diverging_palette(),
      edge_color_trans = "linear",
      edge_color_axis_label = "Log abundance ratio",
      edge_color_interval = c(-1.8, 1.8),

      node_label = taxon_names,
      node_label_size_range = c(0.02, 0.03),
      # title = unique(paste(.$data$diff_read$treatment_1, "vs.", .$data$diff_read$treatment_2)),
      aspect_ratio = 3/2,
      output_file = outfile
    )
}

draw_buffer_compare_heattree <- function(taxmap, amplicon, outfile) {
  {
    value_cols <- names(taxmap$data$read_count) %>% keep(endsWith, amplicon)
    taxmap$data$read_diff <- metacoder::compare_groups(
      taxmap,
      "read_count",
      cols = value_cols,
      groups = value_cols %>%
        str_match("(2015|2016)_(Xpedition|LifeGuard)") %>% {
          paste(.[,3], .[,2])
        },
      func = log_read_ratio,
      combinations = list(
        c("Xpedition 2016", "Xpedition 2015"),
        c("Xpedition 2016", "LifeGuard 2016"),
        c("Xpedition 2015", "LifeGuard 2016")
      )
    )
    theme_update(panel.border = element_blank())
    taxa::filter_taxa(
      taxmap,
      taxmap$data$read_count %>% select_at(value_cols) %>% apply(MARGIN = 1, FUN = max) > 0.001,
      drop_obs = TRUE,
      reassign_obs = FALSE
    ) %>%
      metacoder::heat_tree_matrix(
        "read_diff",
        node_label = taxon_names,
        node_color = read_ratio,
        node_size = rowMeans(.$data$read_count[-1]),
        node_size_range = c(0.0001, 0.03),
        node_size_axis_label = "Read abund.",
        node_label_size_range = c(0.01, 0.05),
        node_color_range = metacoder::diverging_palette(),
        node_color_axis_label = "Log abund. ratio",
        node_color_interval = c(-3, 3),
        node_color_trans = "linear",
        # edge_size = rowMeans(.$data$read_count[-1]),
        # edge_size_axis_label = "Read abundance",
        # edge_size_range = c(0.001, 0.03),
        # edge_size_trans = "linear",
        aspect_ratio = 1,
        layout = "davidson-harel",
        initial_layout = "reingold-tilford",
        output_file =
      )
  }
}
