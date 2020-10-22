#### Output helper functions ####
# These functions format in-line results or help to generate figures and tables
# directly in the Rmarkdown for the paper, or in drake_plan_light.R
# author Brendan Furneaux

# Format a number with "k" or "M" suffic for thousands or millions, as appropriate
k_or_M <- function(x, ..., .sep = " ", .function = format) {
  case_when(
    abs(x) > 1e6 ~ paste(.function(x / 1e6, ...), "M", sep = .sep),
    abs(x) > 1e3 ~ paste(.function(x / 1e3, ...), "k", sep = .sep),
    is.na(x) ~ NA_character_,
    TRUE ~ paste(.function(x, ...), " ", sep = .sep)
  )
}

# Format a number as a percent
percent <- function(x, ...) {
  paste0(formatC(as.numeric(x) * 100, ...), "%")
}

# Format a vector of values to be used as in-line text
text_list <- function(x, ...) {
  x = format(x, ...)
  glue::glue_collapse(x, sep = ", ", last = if (length(x) > 2) ", and " else " and ")
}

# Get R^2 values for ASV or OTU read counts between two different sequencing runs
compR2 <- function(x, y, type) {
  xamp <- plyr::mapvalues(
    x,
    datasets$seq_run,
    datasets$amplicon,
    FALSE
  )
  yamp <- plyr::mapvalues(
    y,
    datasets$seq_run,
    datasets$amplicon,
    FALSE
  )
  x <- plyr::mapvalues(
    x,
    datasets$seq_run,
    paste(
      datasets$tech,
      datasets$machine),
    FALSE
  )
  y <- plyr::mapvalues(
    y,
    datasets$seq_run,
    paste(
      datasets$tech,
      datasets$machine),
    FALSE
  )
  formatC(
    filter(
      readd(comparisons),
      x == x_var,
      xamp == x_amplicon,
      y == y_var,
      yamp == y_amplicon,
      type == !!type
    )$r.squared,
    digits = 2,
    format = "f")
}

# Make plots showing relationship between ASV or OTU read counts in different
# sequencing runs
read_comparison <- function(multi_table, comparisons, type) {
  multi_table <- filter(multi_table, type == !!type)
  comparisons <- filter(comparisons, type == !!type)
  p <- multi_table %>%
    filter(x > 0, y > 0) %>%
    ggplot(aes(y = y, x = x)) +
    geom_point(alpha = 0.2, shape = 1) +
    scale_x_log10(
      breaks = c(1, 1000, 1000000),
      labels = c("1", "1k", "1M"),
      minor_breaks = c(10, 100, 10000, 100000),
      limits = c(1, 1e6)
    ) +
    scale_y_log10(
      breaks = c(1, 1000, 1000000),
      labels = c("1", "1k", "1M"),
      minor_breaks = c(10, 100, 10000, 100000),
      limits = c(1, 1e6)
    ) +
    coord_equal() +
    xlab(paste(type, "read count")) +
    ylab(paste(type, "read count")) +
    geom_rug(aes(y = y), sides = "l", alpha = 0.2,
             data = multi_table %>% filter(x == 0, y > 0)) +
    geom_rug(aes(x = x), sides = "b", alpha = 0.2,
             data = multi_table %>% filter(y == 0, x > 0)) +
    ggnomics::facet_nested(y_amplicon + y_var ~ x_amplicon + x_var, switch = "both", space = "free", nest_line = TRUE) +
    stat_smooth(method = "loess", formula = y ~ x) +
    geom_abline(aes(intercept = `(Intercept)`, slope = 1),
                linetype = 2L, alpha = 0.5, data = comparisons) +
    geom_text(
      aes(label = label),
      data = comparisons,
      x = 6,
      y = 0.7,
      color = "black",
      hjust = 1,
      parse = TRUE,
      inherit.aes = FALSE
    ) +
    theme(strip.placement = "outside", strip.background = element_blank())
  remove_empty_facets(p)
}

remove_empty_facets <- function(p) {
  g <- ggplotGrob(p)
  g$grobs[g$layout$name %in% c("panel-1-2", "panel-1-3", "panel-2-3")] <- NULL
  g$layout <- g$layout[!(g$layout$name %in% c("panel-1-2", "panel-1-3", "panel-2-3")),]
  g
}

# function to call inside group_map while making OTU/ASV read count comparison
# plots.
# Because of several levels of NSE, this didn't work well as an anonymous function.
choosevars <- function(d, g, .data) {
  .data %>%
    filter(type == g$type) %>%
    select(
      type,
      x = !!paste(g$x_var, "-", g$x_amplicon),
      y = !!paste(g$y_var, "-", g$y_amplicon)
    ) %>%
    filter(x > 0 | y > 0) %>%
    mutate(
      x_var = g$x_var,
      x_amplicon = g$x_amplicon,
      y_var = g$y_var,
      y_amplicon = g$y_amplicon
    )
}

scale_color_strategy <- function(...) {
  scale_color_manual(
    values = c(
      "PacBio Long" = "#e7298a",
      "PacBio Short" = "#1b9e77",
      "Illumina Short" = "#7570b3",
      "Ion Torrent Short" = "#d95f02"
    ),
    ...
  )
}

scale_color_tech <- function(...) {
  scale_color_manual(
    values = c(
      "PacBio" = "#e7298a",
      "Illumina" = "#7570b3",
      "Ion Torrent" = "#d95f02"
    ),
    ...
  )
}

scale_color_read <- function(...) {
  scale_color_manual(
    values = c(
      "PacBio Long" = "#e7298a",
      "PacBio Short" = "#1b9e77",
      "Illumina Short R1" = "#7570b3",
      "Illumina Short R2" = "#66a61e",
      "Ion Torrent Short" = "#d95f02"
    ),
    ...
  )
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
            scales::log_breaks(base = base),
            domain = c(1e-100, Inf))
}



# Present the results of venndata() as a table
venn_table <- function(venndata, var, caption, caption.short) {
  var <- enquo(var)

  linesep <- venndata %>%
    filter(!is.na(!!var)) %>%
    select(starts_with("reads")) %>%
    mutate_all(nchar) %>%
    as.matrix %>%
    is_greater_than(0) %>%
    rowSums %>%
    rle() %$%
    map(lengths, ~ c(rep("", . - 1), "\\addlinespace")) %>%
    unlist() %>%
    head(length(.) - 2) %>%
    c("\\midrule")

  empty_background <- if (options("knitr.table.format") == "html") "#777777" else "gray"

  venndata %>%
    rownames_to_column(" ") %>%
    filter(!is.na(!!var)) %>%
    mutate_all(
      ~ifelse(
        . == "",
        kableExtra::cell_spec(
          "  ",
          background = empty_background,
          background_as_tile = FALSE
        ),
        .
      )
    ) %>%
    mutate_at(" ", str_replace, "[01]+", "") %>%
    set_names(gsub("_[[:alpha:]]+[-._]\\d+$", "", names(.))) %>%
    set_names(kableExtra::linebreak(gsub("[_ ]", "\n", names(.)), align = "c")) %>%
    kable(
      booktabs = TRUE,
      caption = caption,
      caption.short = caption.short,
      linesep = linesep,
      align = ("rcrrcrrcrrcrr"),
      row.names = FALSE,
      escape = FALSE
    ) %>%
    kableExtra::add_header_above(
      header = c(
        " " = 2,
        names(venndata)[-(1)] %>%
          gsub(pattern = paste0("(", as_label(var), "_)?(reads|frac)_"), replacement = "") %>%
          rle() %$%
          set_names(lengths, plyr::mapvalues(values, datasets$seq_run, datasets$amplicon, FALSE))
      )
    ) %>%
    kableExtra::add_header_above(
      header = c(
        " " = 2,
        names(venndata)[-(1)] %>%
          gsub(pattern = paste0("(", as_label(var), "_)?(reads|frac)_"), replacement = "") %>%
          plyr::mapvalues(datasets$seq_run, paste(datasets$tech, datasets$machine), FALSE) %>%
          rle() %$%
          set_names(lengths, values)
      )
    )
}
