# functions to make targets related to venn diagrams
# author Brendan Furneaux

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

# Put the results of venndata() in the format to actually make the diagram
vennplot_tech <- function(table, var, letter) {
  var <- enquo(var)
  venndata <- venndata(table, !!var, cols = c("pb_483", "SH-2257", "is_057"))
  plotdata <- vennplot_data(venndata, !!var)
  venn <- select(plotdata, set, !!var) %>%
    deframe() %>%
    replace_na(0) %>%
    eulerr::venn()
  plot(
    venn,
    legend = FALSE,
    labels = FALSE,
    quantities = plotdata$label,
    main = letter,
    fills = c("grey85", "lightblue", "lightcoral")
  )
}

vennplot_amplicon <- function(table, var, letter) {
  var <- enquo(var)
  venndata <- mutate_if(table, is.numeric, ~ . / sum(.)) %>%
    mutate(long = pb_500,
           short = rowMeans(.[,c("pb_483", "SH-2257", "is_057")])) %>%
    select(seq, long, short) %>%
    venndata(!!var, cols = c("long", "short"))
  plotdata <- vennplot_data(venndata, !!var)
  venn <- select(plotdata, set, !!var) %>%
    deframe() %>%
    replace_na(0) %>%
    eulerr::venn()
  plot(
    venn,
    legend = FALSE,
    labels = FALSE,
    quantities = plotdata$label,
    main = letter,
    fills = c("white", "slateblue1"))
}
