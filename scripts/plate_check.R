read_platemap <- function(filename, sheet, skip = 2) {
  # read cells individually
  cells <- tidyxl::xlsx_cells(filename, sheet)
  # find the header row
  header_line <- cells %>%
    filter(character %in% c("Year", "Sample", "Plate", "row")) %$%
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
  primer.pair <- cells %>%
    filter(row < header_line,
           str_detect(character, "ITS")) %>%
    # parse out the primer names
    select(col, character) %>%
    tidyr::extract(character, c("Forward", "Reverse"),
                   regex = "^([:alnum:]+) ?and ?([:alnum:]+).*") %>%
    tidyr::unite(Primer.Pair, Forward, Reverse, remove = FALSE) %>%
    # find the columns they refer to
    mutate(col.end = dplyr::lead(col) - 1,
           col.end = replace_na(col.end, ncol(raw))) %>%
    rowwise()
  # Get the data which applies to all primer pairs.
  out_all <- raw[1:(min(primer.pair$col) - 1)]
  # tidy the whole thing
  do(primer.pair,
     bind_cols(out_all,
               Primer.Pair = rep(.$Primer.Pair, nrow(raw)),
               raw[(.$col):(.$col.end)])) %>%
    # get rid of empty columns
    discard(~all(is.na(.))) %>%
    # make nice names
    dplyr::rename(Comment = V1,
           Well = row,
           DNA.Conc = "Conc. (ng/µl)",
           DNA.Tot = DNA,
           PCR.Conc = "ng/µl",
           PCR.Tot = "Total ng (39 µl)",
           Pool.Vol = "µl Product") %>%
    select(-"gITS7 ITS4", -"No.") %>%
    # Sample name is actually three pieces of information
    tidyr::extract("Sample", c("Site", "X", "Qual"),
                   "(V1P[13]|T[12])[ -][ST](\\d+)(\\D*)$") %>%
    mutate_at("X", as.integer) %>%
    # Sites were recorded differently in the two years
    mutate(Site = case_when(Site == "V1P1" ~ "Ang",
                            Site == "V1P3" ~ "Gan",
                            Site == "T1" ~ "Gan",
                            Site == "T2" ~ "Ang",
                            TRUE ~ NA_character_)) %>%
    # For some site/years, coordinates were measured from the edge of
    # the plot rather than numbered from 1--25 
    group_by(Year, Site) %>%
    mutate(X = X - if (!any(is.na(X)) && max(X) > 25) 12 else 0) %>%
    ungroup %>%
    # Well is also two pieces of information
    tidyr::extract("Well", c("Row", "Column"),
                   "([A-H])(\\d+)", remove = FALSE) %>%
    mutate_at("Row", factor, levels = LETTERS[1:8]) %>%
    mutate_at(c("Column", "Row"), as.integer) %>%
    # The comment column tells us about blanks and controls.
    mutate(sample_type =
             factor(Comment,
                    levels = c("Blank", "Pos", "Pos. Kontroll")) %>%
             fct_collapse(Blank = "Blank", Pos = c("Pos", "Pos. Kontroll")) %>%
             fct_explicit_na("Sample")) %>%
    mutate_at("Plate", formatC, width = 3, flag = "0")
}
