# functions for making heat trees
# author Brendan Furneaux

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
