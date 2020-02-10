# compile data to make a table 
venndata <- function(data, vennvar, cols = names(data)[-1]) {
  vennvar <- enquo(vennvar)
  vennvarfrac <- paste0(as_label(vennvar), "_frac") %>% parse_quo(env = global_env())
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
      names_pattern = "([[:alpha:]]+[-._]\\d+)_(.+)"
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
      !!vennvarfrac := if_else(grepl("^ *[0.]+ *$", !!vennvar), "", !!vennvarfrac),
      !!vennvar := na_if(!!vennvar, 0),
      reads_frac = if_else(grepl("^ *[0.]+ *$", reads), "", reads_frac),
      reads = if_else(grepl("^ *[0.]+ *$", reads), "", reads)
    ) %>%
    select(set, seq_run, !!vennvar, !!vennvarfrac, reads, reads_frac) %>% {
      left_join(
        group_by(., set) %>%
          summarize(!!vennvar := max(!!vennvar, na.rm = TRUE) %>% ifelse(is.finite(.), ., NA)),
        select(., set, seq_run, !!vennvarfrac, reads, frac = reads_frac) %>%
          pivot_wider(
            names_from = "seq_run",
            values_from = c(as_label(vennvarfrac), "reads", "frac")
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
      reads_frac = select(., starts_with("frac")) %>%
        pmap_chr(paste, sep = "/") %>%
        gsub("/+", "/", .) %>%
        gsub("(^/|/$)", "", .),
      label = str_c(!!var, reads_frac, sep = "\n") %>%
        replace_na("-")
    )
}

venn_table <- function(venndata, var, caption) {
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
    set_names(gsub("_", " ", names(.))) %>%
    kable(
      booktabs = TRUE,
      caption = caption,
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
