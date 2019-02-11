library(magrittr)
library(tidyverse)
library(readxl)
library(tidyxl)

if (interactive()) {
  base.dir <- str_extract(getwd(), ".+oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- here("config")
  plate.file <- file.path(lab.dir, "Brendan_soil2.xlsx")
  out.file <- file.path(lab.dir, "short.platekey.csv")
}

xlsx_cells(plate.file, "PCR-Setup") %>%
  filter(data_type == "character")
  

plate <- read_xlsx(plate.file, "Concentration samples", skip = 2, col_names = TRUE) %>%
  rename(Comment = X__1, Well = row) %>%
  select (-X__2) %>%
  # Sample number consists of site, location along transect, and possible qualifier
  tidyr::extract("Sample", c("Site", "X", "Qual"), "(V1P[13]|T[12])[ -][ST](\\d+)(\\D*)$") %>%
  mutate_at("X", as.integer) %>%
  # Sites were recorded differently in the two years
  mutate(Site = case_when(Site == "V1P1" ~ "Ang",
                          Site == "V1P3" ~ "Gan",
                          Site == "T1" ~ "Gan",
                          Site == "T2" ~ "Ang",
                          TRUE ~ NA_character_)) %>%
  # For some site/years, coordinates were measured from the edge of the plot rather than numbered from 1--25 
  group_by(Year, Site) %>%
  mutate(X = X - if (!any(is.na(X)) && max(X) > 25) 12 else 0) %>%
  ungroup %>%
  tidyr::extract("Well", c("Row", "Column"), "([A-H])(\\d+)", remove = FALSE) %>%
  mutate_at("Row", factor, levels = LETTERS) %>%
  mutate_at(c("Column", "Row"), as.integer) %>%
  mutate(sample_type = factor(Comment, levels = c("Blank", "Pos", "Pos. Kontroll")) %>%
           fct_collapse(Blank = "Blank", Pos = c("Pos", "Pos. Kontroll")) %>%
           fct_explicit_na("Sample"))
