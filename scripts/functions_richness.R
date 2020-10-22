# functions to calculate targets for ASV/OTU richness comparisons

nonzero_asv_counts <- function(seq_table, platemap) {
  tibble::as_tibble(seq_table, rownames  = "file") %>%
  tidyr::extract(
    file,
    c("seq_run", "plate", "well"),
    regex = "([ipS][sbH][-_]\\d{3,4})_([0OT]{2}\\d)_([A-H]1?[0-9])_ITS2"
  ) %>%
  dplyr::semi_join(platemap, by = c("plate", "well")) %>%
  dplyr::select(-plate) %>%
  dplyr::group_by(seq_run, well) %>%
  dplyr::summarize_all(sum) %>%
  tidyr::unite("sample", c(seq_run, well)) %>%
  tidyr::pivot_longer(-sample, names_to = "ASV", values_to = "nread") %>%
  dplyr::filter(nread > 0) %>%
  dplyr::group_by(sample) %>%
  dplyr::filter(sum(nread) >= 100) %>%
  dplyr::ungroup() %>%
  dplyr::select(sample, nread) %>%
  tibble::deframe() %>%
  split(names(.))
}

nonzero_otu_counts <- function(seq_table_file, platemap) {
  read_tsv(
    seq_table_file,
    col_types = cols(.default = col_integer(), `#OTU ID` = col_character())
  ) %>%
    dplyr::rename(OTU = `#OTU ID`) %>%
    tidyr::pivot_longer(-OTU, names_to = "sample", values_to = "nread") %>%
    dplyr::filter(nread > 0) %>%
    tidyr::extract(
      col = sample,
      into = c("seq_run", "plate", "well"),
      regex = "([piS][bsH][-_]\\d{3,4})_(\\d{3})([A-H]1?[0-9])"
    ) %>%
    dplyr::semi_join(readd(platemap), by = c("plate", "well")) %>%
    tidyr::unite("sample", seq_run, well) %>%
    dplyr::group_by(sample, OTU) %>%
    dplyr::summarize_at("nread", sum) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(sum(nread) >= 100) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, nread) %>%
    tibble::deframe() %>%
    split(names(.))
}

reshape_rarefy_data <- function(sample_rarefy) {
  purrr::map_dfr(sample_rarefy, tibble::enframe) %>%
    tidyr::extract(
      name,
      c("seq_run", "well"),
      regex = "([iSp][sHb][-_][0-9]+)_([A-H]1?[0-9])"
    ) %>%
    dplyr::left_join(datasets, by = "seq_run") %>%
    dplyr::mutate(strategy = paste(tech, amplicon)) %>%
    dplyr::mutate_at(
      "strategy",
      factor,
      ordered = TRUE,
      levels = c("Ion Torrent Short", "Illumina Short", "PacBio Short",
                 "PacBio Long")
    ) %>%
    dplyr::select(well, strategy, value) %>%
    dplyr::left_join(., ., by = "well", suffix = c("_x", "_y")) %>%
    dplyr::filter(strategy_x < strategy_y)
}

fit_deming_slope <- function(rarefy_data) {
  dplyr::group_nest(rarefy_data, strategy_x, strategy_y, .key = "deming") %>%
    dplyr::mutate(
      deming = purrr::map(
        deming,
        ~ deming::deming(value_y ~ value_x - 1, data = .) %$%
          c(
            list(slope = coefficients[2]),
            as.list(ci[2,]),
            list(r_squared = cor(model)[1,2])
          ) %>%
          tibble::as_tibble()
      )
    ) %>%
    tidyr::unnest(deming) %>%
    dplyr::mutate(
      slope_label = sprintf("m==%.2f(%.2f-%.2f)",
                            slope, `lower 0.95`, `upper 0.95`),
      r2_label = sprintf("R^2==%.2f", r_squared)
    )
}


make_alpha_plot <- function(rarefy_data, demingfits, cluster) {
  ggplot(rarefy_data, aes(x = value_x, y = value_y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray50") +
    geom_point(shape = 1, alpha = 0.8) +
    geom_abline(
      aes(slope = slope, intercept = 0),
      color = "blue",
      alpha = 0.8,
      size = 0.5,
      data = demingfits
    ) +
    geom_abline(
      aes(slope = `lower 0.95`, intercept = 0),
      color = "blue",
      alpha = 0.5,
      size = 0.2,
      data = demingfits
    ) +
    geom_abline(
      aes(slope = `upper 0.95`, intercept = 0),
      color = "blue",
      alpha = 0.5,
      size = 0.2,
      data = demingfits
    ) +
    geom_text(
      aes(x = 0, y = 60, label = slope_label),
      hjust = 0,
      # color = "blue",
      size = 2,
      parse = TRUE,
      data = demingfits
    ) +
    geom_text(
      aes(x = 0, y = 52, label = r2_label),
      hjust = 0,
      # color = "blue",
      size = 2,
      parse = TRUE,
      data = demingfits
    ) +
    facet_grid(strategy_y ~ strategy_x, switch = "both") +
    theme_bw() +
    coord_equal() +
    theme(strip.background = element_blank(),
          strip.placement = "outside") +
    xlab(paste(cluster, "richness")) +
    ylab(paste(cluster, "richness"))
}

make_accum_plot <- function(sample_iNEXT, cluster_type) {
  purrr::flatten(sample_iNEXT) %>%
    dplyr::bind_rows(.id = "sample") %>%
    tidyr::extract(
      sample,
      into = c("seq_run", "well"),
      regex = "([ipS][sbH][-_][0-9]{3,4})_([A-H]1?[0-9])",
      remove = FALSE
    ) %>%
    dplyr::left_join(datasets, by = "seq_run") %>%
    dplyr::mutate(strategy = paste(tech, amplicon)) %>%
    dplyr::filter(method != "extrapolated") %>%
    ggplot(aes(x = m, y = qD, color = strategy)) +
    # lines for each sample
    geom_line(aes(group = sample), alpha = 0.08) +
    # small points to represent the actual observed read depth and richness
    geom_point(data = ~dplyr::filter(., method == "observed"),
               alpha = 0.8, size = 1, shape = 1) +
    xlab("Number of rarefied reads") +
    ylab(paste(cluster_type, "richness")) +
    scale_color_strategy(name = NULL) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
}
