#### Output helper functions ####
# These functions format in-line results or help to generate figures and tables
# directly in the Rmarkdown for the paper, or in drake_plan_light.R

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

# compile data to make a table with the information for a Venn diagram
venndata <- function(data, vennvar, cols = names(data)[-1]) {
  vennvar <- enquo(vennvar)
  vennvarfrac <- paste0(as_label(vennvar), "_frac") %>% rlang::parse_quo(env = globalenv())
  colpattern <- paste0("(", paste(cols, collapse = "|"), ")_(reads|", as_label(vennvar), ")")
  data %>%
    select_at(c("seq", cols)) %>%
    mutate_at(cols, list(found = ~. > 0)) %>%
    group_by_at(vars(ends_with("_found"))) %>%
    summarize_at(cols, rlang::list2(!!vennvar := ~sum(. > 0), reads = sum)) %>%
    ungroup() %>%
    mutate_at(paste0(cols, "_found"), as.integer) %>%
    full_join(
      do.call(crossing, select_at(., paste0(cols, "_found"))),
      by =  paste0(cols, "_found")
    ) %>%
    mutate_if(is.integer, replace_na, 0L) %>%
    mutate_if(is.double, replace_na, 0L) %>%
    tidyr::unite(col = "set", seq_along(cols), sep = "") %>%
    filter(set != strrep(0, length(cols))) %>%
    pivot_longer(
      -1,
      names_to = c("seq_run", "what"),
      values_to = "value",
      names_pattern = colpattern
    ) %>%
    pivot_wider(names_from = "what", values_from = "value") %>%
    group_by(seq_run) %>%
    mutate_if(is.numeric, list(frac = ~./sum(.))) %>%
    bind_rows(
      group_by(., seq_run) %>%
        summarize_if(is.numeric, sum, na.rm = TRUE) %>%
        mutate(set = "Total")
    ) %>%
    mutate_at(
      vars(ends_with("_frac")),
      na_if,
      0
    ) %>%
    mutate_at(
      vars(ends_with("_frac")),
      formatC,
      digits = 2,
      format = "f",
      drop0trailing = FALSE
    ) %>%
    mutate_at(
      "reads",
      k_or_M,
      .sep = "",
      .function = formatC,
      format = "fg",
      digits = 2
    ) %>%
    mutate(
      !!vennvarfrac := if_else(grepl("^ *NA *$", !!vennvarfrac), "", !!vennvarfrac),
      !!vennvar := na_if(!!vennvar, 0),
      reads_frac = if_else(grepl("^ *[0.]+ *$", reads), "", reads_frac),
      reads = if_else(grepl("^ *[0.]+ *$", reads), "", reads)
    ) %>%
    select(set, seq_run, !!vennvar, !!vennvarfrac, reads, reads_frac) %>% {
      left_join(
        group_by(., set) %>%
          summarize(!!vennvar := max(!!vennvar, na.rm = TRUE) %>% ifelse(is.finite(.), ., NA)),
        select(., set, seq_run, !!vennvarfrac, reads, reads_frac) %>%
          pivot_wider(
            names_from = "seq_run",
            values_from = c(as_label(vennvarfrac), "reads", "reads_frac")
          ),
        by = "set"
      )
    } %>%
    mutate(n = str_count(set, "1") + 5 * (set == "Total")) %>%
    arrange(n, desc(set)) %>% {
      bind_cols(
        select(., set, !!vennvar),
        map(cols, ~select(.y, ends_with(.x)), .)
      )
    } %>%
    column_to_rownames("set") %>%
    inset("Total", 1, nrow(data))
}

# Put the results of venndata() in the format to actually make the diagram
vennplot_data <- function(venndata, var) {
  var <- enquo(var)
  venndata %>%
    rownames_to_column("set") %>%
    filter(set != "Total") %>%
    mutate(set = map_chr(
      set,
      ~ strsplit(., "")[[1]] %>%
        as.integer() %>%
        as.logical() %>%
        magrittr::extract(LETTERS[seq_along(.)], .) %>%
        paste(collapse = "&"))) %>%
    mutate(
      reads_frac = select(., starts_with("reads_frac")) %>%
        pmap_chr(paste, sep = "/") %>%
        gsub("/+", "/", .) %>%
        gsub("(^/|/$)", "", .),
      label = str_c(!!var, reads_frac, sep = "\n") %>%
        replace_na("-")
    )
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

# get the ENA accession number from the submission receipt
get_reads_accno <- function(sample, reads_dir = "output/ENA/reads") {
  file.path(reads_dir, sample, "submit", "receipt.xml") %>%
    xml2::read_xml() %>%
    xml2::xml_find_first("RUN") %>%
    xml2::xml_attr("accession")
}

lookup_ncbi_taxon <- function(taxon, Taxonomy, rank, label, ...) {
  uid <- taxize::get_uid_(
    taxon,
    key = readLines("ENTREZ_KEY"),
    messages = FALSE
  )[[1]]
  if (is.null(uid)) {
    return(tibble::tibble(label = label, taxon = taxon, targetTaxonomy = Taxonomy, targetRank = rank))
  }
  uid$label <- label
  uid$taxon <- taxon
  uid$targetTaxonomy <- Taxonomy
  uid$targetRank <- rank
  uid$Taxonomy <- taxize::classification(
    uid$uid,
    db = "ncbi"
  ) %>%
    purrr::map_chr(~paste(.$name, collapse = ";"))
  uid$dist <- stringdist::stringdist(uid$Taxonomy, Taxonomy, method = "jaccard", q = 5)
  if (nrow(uid) <= 1) {
    return(uid)
  }
  # if (rank %in% uid$rank) {
  #   uid <- uid[uid$rank == rank, , drop = FALSE]
  #   if (nrow(uid) <= 1) {
  #     return(uid)
  #   }
  # }
  uid
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

log_read_ratio <- function(abund1, abund2) {
  list(read_ratio = log10(mean(abund1) / mean(abund2)) %>% ifelse(is.nan(.), 0, .))
}

log_ASV_ratio <- function(abund1, abund2) {
  list(ASV_ratio = log10(mean(abund1) / mean(abund2)) %>% ifelse(is.nan(.), 0, .))
}
