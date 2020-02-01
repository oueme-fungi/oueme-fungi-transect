read_platemap <- function(filename, sheet, skip = 2) {
  # read cells individually
  cells <- tidyxl::xlsx_cells(filename, sheet)
  # find the header row
  header_line <- cells %>%
    dplyr::filter(character %in% c("Year", "Sample", "Plate", "row")) %$%
    row %>%
    table %>%
    {names(.)[which.max(.)]} %>%
    as.integer
  # read the data
  raw <- readxl::read_excel(filename, sheet = sheet, 
                            skip = header_line - 1,
                            col_names = TRUE, .name_repair = "minimal")
  # the data is not "tidy"; we have groups of columns for the different
  # primer pairs.  Find the cells which define the primer pairs.
  primer_pair <- cells %>%
    dplyr::filter(row < header_line,
           stringr::str_detect(character, "ITS")) %>%
    # parse out the primer names
    dplyr::select(col, character) %>%
    tidyr::extract(character, c("forward", "reverse"),
                   regex = "^([:alnum:]+) ?and ?([:alnum:]+).*") %>%
    tidyr::unite(primer_pair, forward, reverse, remove = FALSE) %>%
    # find the columns they refer to
    dplyr::mutate(col_end = dplyr::lead(col) - 1,
                  col_end = tidyr::replace_na(col_end, ncol(raw))) %>%
    dplyr::rowwise()
  # Get the data which applies to all primer pairs.
  out_all <- raw[1:(min(primer_pair$col) - 1)]
  # tidy the whole thing
  dplyr::do(primer_pair,
     dplyr::bind_cols(out_all,
               primer_pair = rep(.$primer_pair, nrow(raw)),
               raw[(.$col):(.$col_end)])) %>%
    # get rid of empty columns
    purrr::discard(~all(is.na(.))) %>%
    # make nice names
    dplyr::rename(comment = V1,
           well = row,
           dna_conc = "Conc. (ng/µl)",
           dna_tot = DNA,
           pcr_conc = "ng/µl",
           pcr_tot = "Total ng (39 µl)",
           pool_vol = "µl Product",
           year = "Year",
           plate = "Plate") %>%
    dplyr::select(-"gITS7 ITS4", -"No.") %>%
    # Sample name is actually three pieces of information
    tidyr::extract("Sample", c("site", "x", "qual"),
                   "(V1P[13]|T[12])[ -][ST](\\d+)(\\D*)$") %>%
    dplyr::mutate_at("x", as.integer) %>%
    # Sites and qualities were recorded differently in the two years
    dplyr::mutate(
      site = dplyr::case_when(site == "V1P1" ~ "Ang",
                              site == "V1P3" ~ "Gan",
                              site == "T1" ~ "Gan",
                              site == "T2" ~ "Ang",
                              TRUE ~ NA_character_)
    ) %>%
    # For some site/years, coordinates were measured from the edge of
    # the plot rather than numbered from 1--25 
    dplyr::group_by(year, site) %>%
    dplyr::mutate(x = x - if (!any(is.na(x)) && max(x) > 25) 12 else 0) %>%
    dplyr::ungroup() %>%
    # Well is also two pieces of information
    tidyr::extract("well", c("row", "column"),
                   "([A-H])(\\d+)", remove = FALSE) %>%
    dplyr::mutate_at("row", factor, levels = LETTERS[1:8]) %>%
    dplyr::mutate_at(c("column", "row"), as.integer) %>%
    # The comment column tells us about blanks and controls.
    dplyr::mutate(
      sample_type =
        factor(comment,
               levels = c("Blank", "Pos", "Pos. Kontroll")) %>%
        forcats::fct_collapse(Blank = "Blank",
                              Pos = c("Pos", "Pos. Kontroll")) %>%
        forcats::fct_explicit_na("Sample"),
      buffer = dplyr::case_when(
        year == "2015" & sample_type == "Sample" ~ "Xpedition",
        year == "2016" & qual == "X" ~ "Xpedition",
        year == "2016" & sample_type == "Sample" ~ "LifeGuard",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::mutate_at("plate", formatC, width = 3, flag = "0")
}
