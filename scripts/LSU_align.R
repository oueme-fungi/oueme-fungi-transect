library(magrittr)
library(tidyverse)
library(Biostrings)
library(here)

if (interactive()) {
  base.dir <- getwd() %>%
    str_extract(".*oueme-fungi-transect")
  data.dir <- file.path(base.dir, "data")
  lab.dir <- here("config")
  ref.dir <- file.path(base.dir, "reference")
  lsuref <- file.path(ref.dir, "SILVA_132_LSURef_Fungi.fasta")
  lsuout <- file.path(ref.dir, "SILVA_132_LSURef_Fungi.align.fasta")
}

lsu <- readRNAStringSet(lsuref)
bases <- str_subset(alphabet(lsu), "\\w")
sublsu <- subseq(lsu, 149001, width = 1000, allow.nonnarrowing = TRUE) %>%
  as.matrix()
lsumatrix <-
  seq(1, min(width(lsu[1]), width(lsu[1]) %/% 1000 * 1000 + 1), by = 1000) %>%
  map(~ { cat(., ": ")
    subseq(lsu, ., width = min(1000, width(lsu[1]) - .)) %>%
        as.matrix %>%
        magrittr::extract(.,,plyr::alply(., .margins = 2, unique) %>%
            map(is_in, bases) %>%
            map_lgl(any) %T>%
              {cat(sum(.), "\n")},
          drop = FALSE)}) %>%
  keep(~ length(.) > 0) %>%
  do.call(cbind, .)

lsumatrix %>%
  apply(1, paste, collapse = "") %>%
  RNAStringSet %>%
  writeXStringSet(filepath = lsuout)

lsuprofile <- 
  apply(lsumatrix, 2, . %>% table %>% vegan::diversity() %>% exp)
lsupresence <-
  apply(lsumatrix, 2, . %>% is_in(bases) %>% sum) /
  nrow(lsumatrix)
  

pllsumatrix %<>% keep(~ length(.) > 0)
lsumatrix %<>% do.call(cbind, .)

length(lsumatrix[[1]])
subseq(lsu, 67001, width = 1000) %>%
  as.matrix %>%
  plyr::alply(2, unique) %>%
  map(is_in, bases) %>%
  map_lgl(any)
map(seq(1,  sublsu[,apply(sublsu, MARGIN = 2, . %>% unique %>% is_in(bases) %>% any),
  drop = FALSE] %>%
  dim
lsu %<>% as.matrix
width(lsu)
lsu[,2]
lsu %<>% maskGaps(min.fraction = 0.99, min.block.width = 1)
writeRNA