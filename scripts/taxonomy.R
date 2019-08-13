
patch_taxa <- function(c12n, patch_file) {
  assertthat::assert_that((assertthat::is.string(patch_file)
                           && file.exists(patch_file))
                          || is.null(patch_file))
  if (!is.null(patch_file)) {
    patch <- readr::read_csv(patch_file,
                             col_types = readr::cols(.default = readr::col_character()))
    assertthat::assert_that(assertthat::has_name(patch, "pattern"),
                            assertthat::has_name(patch, "replacement"))
    c12n <- stringi::stri_replace_all_regex(c12n,
                                            pattern = patch$pattern,
                                            replacement = patch$replacement,
                                            vectorize_all = FALSE)
  }
  return(c12n)
}

# take only up to n sequences from each group.
# if the group is a species (key$rank == "species") then take the single most
# representative sequence (minimum total manhattan 5-mer distance to the other
# sequences.)
# Would it be better to take the longest?  Filter based on which has early LSU
# primer sites?  We'll just have to trust that most of the sequences in the
# database start at the beginning of LSU.
trim_group <- function(table, key, indexcol, seq, n, ...) {
  cat("group with", nrow(table), "sequences\n")
  if (nrow(table) <= n) return(table)
  index <- table[[indexcol]]
  assertthat::assert_that(is.integer(index),
                          all(index > 0),
                          max(index) <= length(seq))
  seq <- ape::as.DNAbin(Biostrings::DNAStringSet(seq[index]))
  is_species <- !is.null(key[["rank"]]) && identical(key[["rank"]], "species")
  if (is_species) {
    kmers <- as.matrix(kmer::kdistance(seq))
    s <- which.min(rowSums(kmers))
  } else {
    s <- sample(1:nrow(table), n)
  }
  return(table[s,,drop = FALSE])
}

formatdb <- function(db = c("unite", "rdp", "silva"),
                     seq_file,
                     patch_file = NULL,
                     maxGroupSize = 10,
                     ...) {
  db = match.arg(db)
  switch(db,
         rdp = formatdb_rdp(seq_file = seq_file,
                                patch_file = patch_file, ...),
         silva = formatdb_silva(seq_file = seq_file,
                                patch_file = patch_file, ...),
         unite = formatdb_unite(seq_file = seq_file,
                                patch_file = patch_file, ...))
}

#### RDP ####
seq_file <- file.path("reference", "rdp.fasta.gz")
patch_file <- file.path("reference", "rdp_patch.csv")
formatdb_rdp <- function(seq_file, patch_file = NULL,
                         maxGroupSize = 10, ...) {
  # read the RDP dataset
  rdp_db <- Biostrings::readDNAStringSet(seq_file) %>% unique()
  
  # define the ranks.  The aim is only to search to genus level.
  rdp_ranks <-  c("rootrank", "kingdom", "phylum", "class", "order",
                  "family", "genus", "species")
  rdp_rankpattern <- paste0(";(", paste(rdp_ranks, collapse = "|"), ")")
  
  # classification = c(12 letters)n
  rdp_c12n <- tibble::tibble(name = names(rdp_db),
                             seqidx = seq_along(name)) %>%
    # parse the name string
    tidyr::extract(name, into = c("id", "species", "c12n"),
                   "([:alnum:]+) ([^;\\t]+)[^\\t]*\\tLineage=(.+)") %>%
    # RDP format is "taxon;rank;taxon;rank;"  Remove the ranks.
    dplyr::mutate_at("c12n", stringr::str_replace_all, rdp_rankpattern, "") %>%
    # Apply "patches" to the database
    dplyr::mutate_at("c12n", patch_taxa, patch_file) %>%
    # if the "species" is a species, then add it
    dplyr::mutate(
      genus = stringr::str_replace(c12n, ".+;", ""),
      c12n = ifelse(stringr::str_starts(species, paste0(genus, " ")),
                    paste0(c12n, ";", species),
                    c12n)
    )
  
  # find the unique taxa and give them an index
  rdp_taxa <-
    dplyr::select(rdp_c12n, c12n) %>%
    unique() %>%
    dplyr::mutate(taxid = seq_along(c12n))
  
  # map the index back into the classification list
  rdp_c12n %<>%
    dplyr::left_join(rdp_taxa, by = c("c12n")) %>%
    dplyr::select(-c12n)
  
  # break apart the lineage into ranks
  rdp_taxa %<>%
    dplyr::bind_cols(stringr::str_split_fixed(.$c12n, ";", length(rdp_ranks)) %>%
                       magrittr::set_colnames(rdp_ranks) %>%
                       tibble::as_tibble()) %>%
    dplyr::select(-c12n) %>%
    tidyr::gather(key = "Rank", value = "taxon", !!!rdp_ranks, factor_key = TRUE) %>%
    dplyr::arrange(taxid, Rank) %>%
    # remove "incertae sedis"
    dplyr::filter(stringr::str_detect(taxon, "incertae", negate = TRUE),
                  stringr::str_detect(taxon, "unclassified_", negate = TRUE),
                  stringr::str_detect(taxon, "uncultured_", negate = TRUE),
                  stringr::str_starts(taxon, "[A-Z]"),
                  nchar(taxon) > 0) %>%
    # assign parents
    dplyr::mutate(parent = ifelse(Rank == "rootrank",
                                  NA_character_,
                                  dplyr::lag(taxon))) %>%
    # The RDP database sometimes lists a higher taxon as the genus when the genus
    # is unknown. Remove these.
    dplyr::filter(!purrr::map_lgl(parent == taxon, isTRUE))
  
  # Put the corrected lineage back into the database
  rdp_c12n %<>%
    dplyr::left_join(rdp_taxa %>%
                       dplyr::group_by(taxid) %>%
                       dplyr::summarize(c12n = paste(taxon, collapse = ";")),
                     by = "taxid")
  
  # Check that each taxon has a unique rank and parent
  rdp_taxa %>% dplyr::group_by(taxon) %>%
    dplyr::arrange(taxon) %>%
    dplyr::filter(dplyr::n_distinct(Rank) > 1 | dplyr::n_distinct(parent) > 1) %>%
    {assertthat::assert_that(nrow(.) == 0)}
  
  # put together the lineage (i.e., higher ranks) for each taxon.
  # this is easy to do within each leaf taxon, because we already have them in
  # order
  rdp_taxa %<>%
    dplyr::group_by(taxid) %>%
    dplyr::mutate(c12n = purrr::accumulate(taxon, paste, sep = ";"),
                  Level = seq_along(taxon) - 1) %>%
    dplyr::ungroup() %>%
    # take only the unique taxa
    dplyr::select(-taxid) %>%
    unique() %>%
    # verify that we only have one of each taxon
    assertr::assert(assertr::is_uniq, taxon)
  
  # sometimes we can get more information from the species name than we were given
  # in the lineage; e.g. "uncultured Glomerales  Lineage=Root;Fungi;"
  # If the "species" entry gives us a more specific, but non-conflicting lineage,
  # then use it instead of the given lineage
  rdp_c12n %<>%
    dplyr::mutate_at("species", stringr::str_replace, "uncultured ", "") %>%
    dplyr::left_join(dplyr::select(rdp_taxa, species = taxon, c12n),
                     by = "species") %>%
    dplyr::mutate(c12n = dplyr::if_else(stringr::str_detect(pattern = c12n.x,
                                                            string = c12n.y),
                                        c12n.y,
                                        c12n.x,
                                        missing = c12n.x)) %>%
    dplyr::select(-c12n.x, -c12n.y) %>%
    #find the rank
    dplyr::left_join(dplyr::select(rdp_taxa, c12n, Rank), by = "c12n")
  
  # Format the taxon list as it is wanted for DECIPHER.
  rdp_taxa %<>%
    dplyr::mutate(Index = seq_along(taxon),
                  Rank = as.character(Rank)) %>%
    dplyr::left_join(dplyr::select(., parent = taxon, Parent = Index),
                     by = "parent") %>%
    dplyr::mutate_at("Parent", tidyr::replace_na, 0) %>%
    dplyr::rename(Name = taxon)
  
  if (is.finite(maxGroupSize)) {
    rdp_c12n %<>%
      dplyr::group_by(c12n, Rank) %>%
      dplyr::group_map(trim_group, indexcol = "seqidx",
                       seq = rdp_db, n = maxGroupSize) %>%
      dplyr::ungroup()
  }
  
  list(train = rdp_db[rdp_c12n$seqidx],
       taxonomy = rdp_c12n$c12n,
       rank = rdp_taxa)
}

#### RDP training set ####
seq_file <- file.path("reference", "rdp_train.LSU.fasta.gz")
tax_file <- file.path("reference", "rdp_train_tax.txt")
patch_file <- file.path("reference", "rdp_train_patch.csv")
formatdb_rdp <- function(seq_file, patch_file = NULL,
                         maxGroupSize = 10, ...) {
  # read the RDP dataset
  rdp_train_db <- Biostrings::readDNAStringSet(seq_file) %>% unique()
  rdp_train_data <- tibble(idx = seq_along(rdp_train_db),
                            c12n = names(rdp_train_db)) %>%
    assertr::assert(function(x) stringr::str_detect(x, "[^\\t ]+[\\t ]+([A-Z][A-Za-z0-9 \\-_]+;){3,6}[A-Za-z0-9 \\-_]+$"),
                    c12n) %>%
    dplyr::mutate_at("c12n", stringr::str_replace, "[^\\t ]+[\\t ]+", "")
  
  # define the ranks
  rdp_ranks <-  c("rootrank", "kingdom", "phylum", "class", "order",
                  "family", "genus") %>%
    lapply(taxa::taxon_rank)

  # parse the classification and build a taxonomy
  rdp_train_taxa <- rdp_train_data %>%
    taxa::parse_tax_data(class_cols = "c12n", class_sep = ";")
  
  rank_idx <- rdp_train_taxa$n_supertaxa() + 1
  
  for (i in seq_along(rank_idx)) {
    rdp_train_taxa$taxa[[i]]$rank <- rdp_ranks[[rank_idx[i]]]
  }
  
  rdp_train_taxa <- readr::read_delim(tax_file, delim = "*",
                                      col_names = c("ID", "Name", "Parent",
                                                    "Level", "Rank")) %>%
    dplyr::left_join(dplyr::select(., Parent = ID, ParentName = Name)) %>%
    dplyr::mutate(Ancestor = Parent,
                  c12n = paste0(Name, ";"))
  while (!all(is.na(rdp_train_taxa$Ancestor))) {
    rdp_train_taxa %<>%
      dplyr::left_join(dplyr::select(., Ancestor = ID,
                                     AncestorName = Name,
                                     NextAncestor = Parent), by = "Ancestor") %>%
      dplyr::mutate(c12n = paste(AncestorName, c12n, sep = ";"),
                    Ancestor = NextAncestor) %>%
      dplyr::select(-AncestorName, -NextAncestor)
  }
  
  rdp_train_taxa %<>%
    dplyr::mutate_at("c12n", stringr::str_replace_all,"NA;", "") %>%
    dplyr::filter(!stringr::str_detect(Name, "incertae")) %>%
    dplyr::mutate_at("c12n", stringr::str_replace_all, "[A-Z][a-z]+ incertae sedis;", "")
  
  # define the ranks.  The aim is only to search to genus level.
  rdp_ranks <-  c("rootrank", "domain", "phylum", "class", "order",
                  "family", "genus", "species")
  rdp_rankpattern <- paste0(";(", paste(rdp_ranks, collapse = "|"), ")")
  
  # classification = c(12 letters)n
  rdp_c12n <- tibble::tibble(name = names(rdp_db),
                             seqidx = seq_along(name)) %>%
    # parse the name string
    tidyr::extract(name, into = c("id", "species", "c12n"),
                   "([:alnum:]+) ([^;\\t]+)[^\\t]*\\tLineage=(.+)") %>%
    # RDP format is "taxon;rank;taxon;rank;"  Remove the ranks.
    dplyr::mutate_at("c12n", stringr::str_replace_all, rdp_rankpattern, "") %>%
    # Apply "patches" to the database
    dplyr::mutate_at("c12n", patch_taxa, patch_file) %>%
    # if the "species" is a species, then add it
    dplyr::mutate(
      genus = stringr::str_replace(c12n, ".+;", ""),
      c12n = ifelse(stringr::str_starts(species, paste0(genus, " ")),
                    paste0(c12n, ";", species),
                    c12n)
    )
  
  # find the unique taxa and give them an index
  rdp_taxa <-
    dplyr::select(rdp_c12n, c12n) %>%
    unique() %>%
    dplyr::mutate(taxid = seq_along(c12n))
  
  # map the index back into the classification list
  rdp_c12n %<>%
    dplyr::left_join(rdp_taxa, by = c("c12n")) %>%
    dplyr::select(-c12n)
  
  # break apart the lineage into ranks
  rdp_taxa %<>%
    dplyr::bind_cols(stringr::str_split_fixed(.$c12n, ";", length(rdp_ranks)) %>%
                       magrittr::set_colnames(rdp_ranks) %>%
                       tibble::as_tibble()) %>%
    dplyr::select(-c12n) %>%
    tidyr::gather(key = "Rank", value = "taxon", !!!rdp_ranks, factor_key = TRUE) %>%
    dplyr::arrange(taxid, Rank) %>%
    # remove "incertae sedis"
    dplyr::filter(stringr::str_detect(taxon, "incertae", negate = TRUE),
                  stringr::str_detect(taxon, "unclassified_", negate = TRUE),
                  stringr::str_detect(taxon, "uncultured_", negate = TRUE),
                  stringr::str_starts(taxon, "[A-Z]"),
                  nchar(taxon) > 0) %>%
    # assign parents
    dplyr::mutate(parent = ifelse(Rank == "rootrank",
                                  NA_character_,
                                  dplyr::lag(taxon))) %>%
    # The RDP database sometimes lists a higher taxon as the genus when the genus
    # is unknown. Remove these.
    dplyr::filter(!purrr::map_lgl(parent == taxon, isTRUE))
  
  # Put the corrected lineage back into the database
  rdp_c12n %<>%
    dplyr::left_join(rdp_taxa %>%
                       dplyr::group_by(taxid) %>%
                       dplyr::summarize(c12n = paste(taxon, collapse = ";")),
                     by = "taxid")
  
  # Check that each taxon has a unique rank and parent
  rdp_taxa %>% dplyr::group_by(taxon) %>%
    dplyr::arrange(taxon) %>%
    dplyr::filter(dplyr::n_distinct(Rank) > 1 | dplyr::n_distinct(parent) > 1) %>%
    {assertthat::assert_that(nrow(.) == 0)}
  
  # put together the lineage (i.e., higher ranks) for each taxon.
  # this is easy to do within each leaf taxon, because we already have them in
  # order
  rdp_taxa %<>%
    dplyr::group_by(taxid) %>%
    dplyr::mutate(c12n = purrr::accumulate(taxon, paste, sep = ";"),
                  Level = seq_along(taxon) - 1) %>%
    dplyr::ungroup() %>%
    # take only the unique taxa
    dplyr::select(-taxid) %>%
    unique() %>%
    # verify that we only have one of each taxon
    assertr::assert(assertr::is_uniq, taxon)
  
  # sometimes we can get more information from the species name than we were given
  # in the lineage; e.g. "uncultured Glomerales  Lineage=Root;Fungi;"
  # If the "species" entry gives us a more specific, but non-conflicting lineage,
  # then use it instead of the given lineage
  rdp_c12n %<>%
    dplyr::mutate_at("species", stringr::str_replace, "uncultured ", "") %>%
    dplyr::left_join(dplyr::select(rdp_taxa, species = taxon, c12n),
                     by = "species") %>%
    dplyr::mutate(c12n = dplyr::if_else(stringr::str_detect(pattern = c12n.x,
                                                            string = c12n.y),
                                        c12n.y,
                                        c12n.x,
                                        missing = c12n.x)) %>%
    dplyr::select(-c12n.x, -c12n.y) %>%
    #find the rank
    dplyr::left_join(dplyr::select(rdp_taxa, c12n, Rank), by = "c12n")
  
  # Format the taxon list as it is wanted for DECIPHER.
  rdp_taxa %<>%
    dplyr::mutate(Index = seq_along(taxon),
                  Rank = as.character(Rank)) %>%
    dplyr::left_join(dplyr::select(., parent = taxon, Parent = Index),
                     by = "parent") %>%
    dplyr::mutate_at("Parent", tidyr::replace_na, 0) %>%
    dplyr::rename(Name = taxon)
  
  if (is.finite(maxGroupSize)) {
    rdp_c12n %<>%
      dplyr::group_by(c12n, Rank) %>%
      dplyr::group_map(trim_group, indexcol = "seqidx",
                       seq = rdp_db, n = maxGroupSize) %>%
      dplyr::ungroup()
  }
  
  list(train = rdp_db[rdp_c12n$seqidx],
       taxonomy = rdp_c12n$c12n,
       rank = rdp_taxa)
}

#### Silva ####
tax_file <- file.path("reference", "silva_tax.txt")
seq_file <- file.path("reference", "silva.LSU.fasta.gz")
patch_file <- file.path("reference", "silva_patch.csv")
formatdb_silva <- function(tax_file, seq_file, patch_file = NULL,
                           maxGroupSize = 10, ...) {
  silva_ranks <- c("rootrank", "domain", "major_clade", "superkingdom", "kingdom", "subkingdom", "infrakingdom", "superphylum", "phylum", "subphylum", "infraphylum", "class", "subclass", "order", "family", "subfamily", "genus", "species")
  
  silva_outranks <- c("rootrank", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  silva_taxa <-
    readr::read_tsv(tax_file,
                    col_names = c("c12n", "Index", "Rank", "remark", "release"),
                    col_types = readr::cols(Index = readr::col_integer(),
                                            .default = readr::col_character())) %>%
    dplyr::filter(!stringr::str_detect(c12n, "Incertae Sedis;$"),
                  stringr::str_detect(c12n, "Eukaryota")) %>%
    dplyr::mutate_at("c12n", patch_taxa, patch_file) %>%
    dplyr::bind_rows(list(c12n = "", Index = 1, Rank = "rootrank")) %>%
    dplyr::mutate(c12n = paste0("Root;", c12n),
                  Name = stringr::str_match(c12n, "([^;]+);$")[,2],
                  Rank = factor(Rank, levels = silva_ranks),
                  c12n = stringr::str_replace_all(c12n, "[^;]*Incertae Sedis;", ""),
                  parentc12n = stringr::str_replace(c12n, "[^;]+;$", ""),
                  Level = lengths(strsplit(parentc12n, ";"))) %>%
    dplyr::filter(!duplicated(Name)) %>%
    dplyr::left_join(dplyr::select(., parentc12n = c12n, Parent = Index), by = "parentc12n") %>%
    assertr::verify(!is.na(Parent) | Index == 1) %>%
    assertr::assert(assertr::is_uniq, Index, c12n, Name) %>%
    dplyr::mutate_at("Parent", tidyr::replace_na, 0) %>%
    dplyr::arrange(Index)
  silva_taxa %>%
    dplyr::select(c12n, Index) %>%
    dplyr::mutate(Name = str_split(c12n, ";")) %>%
    tidyr::unnest(Name) %>%
    dplyr::filter(Name != "") %>%
    dplyr::group_by(Index) %>%
    dplyr::mutate(c12n = purrr::accumulate(Name, paste, sep = ";") %>%
                    paste0(";")) %>%
    dplyr::left_join(silva_taxa %>% select(c12n, Rank)) %>%
    dplyr::filter(Rank %in% silva_outranks) %>%
    dplyr::summarize(Rank = last(Rank),
                     c12n = paste(c(Name, ""), collapse = ";")) %>%
    View()
  
  silva_db <- Biostrings::readRNAStringSet(seq_file) %>% unique()
  
  silva_c12n <- 
    tibble::tibble(name = names(silva_db),
                   seqidx = seq_along(name)) %>%
    tidyr::extract(name, c("ID.prefix", "ID.num", "start", "end", "c12n", "species"),
                   "^([:alpha:]+)(\\d+)\\.(\\d+).(\\d+)\\s+(.+;)([^;]+)$") %>%
    dplyr::filter(stringr::str_detect(c12n, "Eukaryota")) %>%
    dplyr::mutate_at("c12n", patch_taxa, patch_file) %>%
    dplyr::mutate(c12n = paste0("Root;", c12n),
                  c12n = stringr::str_replace_all(c12n, "[^;]*Incertae Sedis;", "")) %>%
    assertr::verify(c12n %in% silva_taxa$c12n) %>%
    dplyr::left_join(dplyr::select(silva_taxa, c12n, Index, Rank, Name, Level),
                     by = "c12n") %>%
    dplyr::mutate(is_species = Rank == "genus" &
                    stringr::str_starts(species, paste0(Name, " ")),
                  c12n = ifelse(is_species,
                                paste0(c12n, species, ";"),
                                c12n),
                  Rank = dplyr::if_else(is_species,
                                factor("species", levels = silva_ranks),
                                Rank))
  
  silva_taxa %<>% dplyr::bind_rows(
    dplyr::filter(silva_c12n, is_species) %>%
      dplyr::select(Name, Index, Level, c12n, species, Rank) %>%
      unique() %>%
      dplyr::mutate(Name = species,
                    Parent = Index,
                    Index = seq_along(Parent) + max(silva_taxa$Index),
                    Level = Level + 1,
                    Rank = as.character(Rank)) %>%
      dplyr::select(c12n, Name, Parent, Index, Level, Rank))
  
  if (is.finite(maxGroupSize)) {
    silva_c12n %<>%
      dplyr::group_by(c12n, Rank) %>%
      dplyr::group_map(trim_group, indexcol = "seqidx",
                       seq = silva_db, n = maxGroupSize) %>%
      dplyr::ungroup()
  } 
  
  list(train = silva_db[silva_c12n$seqidx],
       taxonomy = silva_c12n$c12n,
       rank = silva_taxa)
}

#### UNITE ####
seq_file = file.path("reference", "unite.fasta.gz")
patch_file = file.path("reference", "unite_patch.csv")
formatdb_unite <- function(seq_file, patch_file = NULL,
                           maxGroupSize = 10, ...) {
  unite_db <- Biostrings::readDNAStringSet(seq_file) %>% unique()
  
  unite_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  unite_data <- tibble::tibble(idx = seq_along(unite_db),
                               name = names(unite_db)) %>%
    tidyr::separate(name, c("species", "accno", "sh", "type", "c12n"), sep = "\\|") %>%
    
    dplyr::mutate_at("c12n", patch_taxa, patch_file) %>%
    dplyr::mutate_at("c12n",
                     stringi::stri_replace_all_regex,
                     "[kpcofgs]__",
                     "",
                     vectorize_all = TRUE) %>%
    tidyr::separate(c12n, unite_ranks, sep = ";") %>%
    dplyr::mutate(rootrank = "Root")
  
  unite_taxa <- taxa::parse_tax_data(unite_data,
                                     class_cols = c("rootrank", unite_ranks),
                                     named_by_rank = TRUE)
  
  
  unite_c12n <- tibble::tibble(name = names(unite_db),
                               seqidx = seq_along(name)) %>%
    tidyr::separate(name, c("species", "accno", "sh", "type", "c12n"), sep = "\\|") %>%
    dplyr::mutate_at("c12n", ~ paste0("r__Root;", .))
  
  unite_taxa <- dplyr::select(unite_c12n, c12n) %>%
    unique() %>%
    dplyr::mutate(leafidx = seq_along(c12n))
  
  unite_c12n %<>% dplyr::left_join(unite_taxa, by = "c12n")
  
  unite_patch <- readr::read_csv(patch_file)
  unite_taxa %<>%
    dplyr::mutate_at("c12n", stringi::stri_replace_all_regex,
                     pattern = unite_patch$pattern,
                     replacement = unite_patch$replacement,
                     vectorize_all = FALSE) %>%
    dplyr::mutate_at("c12n", strsplit, ";") %>%
    tidyr::unnest(c12n) %>%
    tidyr::separate(c12n, c("Rank", "taxon"), "__") %>%
    dplyr::mutate_at("Rank", match.arg, unite_ranks, several.ok = TRUE) %>%
    dplyr::mutate_at("Rank", factor, levels = unite_ranks) %>%
    dplyr::filter(taxon != "unidentified",
                  stringr::str_detect(taxon, "_sp$", negate = TRUE),
                  stringr::str_detect(taxon, "Incertae_sedis", negate = TRUE)) %>%
    dplyr::mutate(parent = ifelse(Rank == "rootrank",
                                  NA_character_,
                                  dplyr::lag(taxon))) %>%
    dplyr::group_by(leafidx)
  
  unite_c12n %<>% dplyr::select(-c12n) %>%
    dplyr::left_join(unite_taxa %>%
                       dplyr::summarize(c12n = paste(taxon, collapse = ";"),
                                        Rank = dplyr::last(Rank)),
                     by = "leafidx") %>%
    dplyr::filter(c12n != "Root") %>%
    assertr::assert(function(x) startsWith(x, "Root;"), c12n)
  
  unite_taxa %<>%
    dplyr::mutate(Level = seq_along(Rank) - 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-leafidx) %>%
    unique() %>%
    assertr::assert(assertr::is_uniq, taxon)
  
  unite_taxa %<>% 
    dplyr::mutate(Index = seq_along(taxon)) %>%
    dplyr::rename(Name = taxon) %>%
    dplyr::left_join(dplyr::select(., parent = Name, Parent = Index)) %>%
    assertr::verify(!is.na(Parent) | Level == 0) %>%
    dplyr::mutate_at("Parent", tidyr::replace_na, 0)
  
  if (is.finite(maxGroupSize)) {
    unite_c12n %<>%
      dplyr::group_by(c12n, Rank) %>%
      dplyr::group_map(trim_group, indexcol = "seqidx",
                       seq = unite_db, n = maxGroupSize) %>%
      dplyr::ungroup()
  }
  
  list(train = unite_db[unite_c12n$seqidx],
       taxonomy = unite_c12n$c12n,
       rank = unite_taxa)
}

# Based on code from DECIPHER vignette "Classify Sequences"
refine_classifier <- function(trainset, taxonomy, rank, out_root,
                              maxIterations = 5,
                              allowGroupRemoval = FALSE,
                              ...) {
  if (!is.null(out_root) && !is.na(out_root)) {
    oldremoves <- list.files(path = dirname(out_root),
                             pattern = paste0(basename(out_root), "\\d+\\.csv"),
                             full.names = TRUE)
    unlink(oldremoves)
  }
  
  # count the sequences in each taxon
  groupCounts <- table(taxonomy)
  u_groups <- names(groupCounts)
  remove <- logical(length(taxonomy))
  
  probSeqsPrev <- integer() # suspected problem sequences from prior iteration
  
  for (i in seq_len(maxIterations)) {
    cat("Training iteration: ", i, "\n", sep = "")
    
    # Free up memory
    if (exists("trainingSet")) {
      remove(trainingSet)
      gc()
    }
    
    # train the classifier
    trainingSet <- DECIPHER::LearnTaxa(trainset[!remove],
                             taxonomy[!remove],
                             rank,
                             ...)
    # look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    # underannotated sequences are those which are placed with high confidence
    # in a subclade of the one they are annotated in.
    # These can be safely removed, even if they are the last sequence
    # with that annotation.
    underannot <- with(trainingSet$problemSequences,
                       startsWith(Predicted, Expected))
    
    if (length(probSeqs) == 0) {
      cat("No problem sequences remaining.\n")
      return(trainingSet)
    } else if (length(probSeqs) == length(probSeqsPrev) &&
               all(probSeqsPrev == probSeqs)) {
      cat("Iterations converged.\n")
      return(trainingSet)
    }
    if (i == maxIterations)
      return(trainingSet)
    
    # accept the new annotations for the underannotated sequences.
    #taxonomy[!remove][probSeqs][underannot] <-
     
    # trainingSet$problemSequences$Predicted[underannot]
    
    readr::write_csv(trainingSet$problemSequences[!underannot,],
                     paste0(out_root, i, ".csv"))
    
    probSeqsPrev <- probSeqs#[!underannot]
    
    # remove any remaining problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE # remove all problem sequences
    index <- index[!underannot] # only re-add the ones which are not simply underannotated
    if (!allowGroupRemoval) {
      # replace any removed groups
      missing <- !(u_groups %in% taxonomy[!remove])
      missing <- u_groups[missing]
      if (length(missing) > 0) {
        index <- index[taxonomy[index] %in% missing]
        remove[index] <- FALSE # don't remove
      }
    }
  }
  # We shouldn't ever get here, but just in case.
  return(trainingSet)
}

taxonomy <- function(seq.table, reference, multithread = FALSE) {
  assertthat::assert_that(file.exists(reference),
                          assertthat::is.readable(reference))
  # we can take a community matrix (in which case the sequences are the column
  # names) or the sequences
  if (is.matrix(seq.table)) {
    nreads <- Matrix::colSums(seq.table)
    seq <- colnames(seq.table)
  } else {
    nreads <- 1
    seq <- as.character(seq.table)
  }
  is.RNA <- stringr::str_detect(seq[1], "[Uu]")
  # seq <- if (is.RNA) Biostrings::RNAStringSet(seq) else Biostrings::DNAStringSet(seq)
  
  tax <- dada2::assignTaxonomy(seqs = seq, 
                               refFasta = reference,
                               multithread = multithread) %>%
    tibble::as_tibble(rownames = "seq") %>%
    # remove taxon rank prefixed from Unite reference
    dplyr::mutate_at(dplyr::vars(-seq), stringr::str_replace, "^[kpcofgs]__", "") %>%
    dplyr::mutate(
      # add Species to RDP reference
      Species = if ("Species" %in% names(.)) Species else NA_character_,
      Species = ifelse(is.na(Genus) | is.na(Species),
                       NA_character_,
                       paste(Genus, Species)),
      # Put the whole classification in one field for FUNGuild
      Taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species,
                       sep = ";") %>%
        stringr::str_replace_all(stringr::fixed(";NA"), ""),
      # Use the finest classification available as a short "name" for trees
      Name = dplyr::coalesce(Species, Genus, Family, Order, Class,
                             Phylum, Kingdom, "unknown") %>%
        stringr::str_replace_all("\\s", "_")) %>%
    # If there are duplicate names, number them.
    dplyr::group_by(Name) %>%
    dplyr::mutate(newname = if (dplyr::n() > 1) {
      paste(Name, seq_along(Name), sep = "_")
    } else {
      Name
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Name) %>%
    dplyr::rename(Name = newname) %>%
    
    # add in the reads if we had them.
    dplyr::left_join(tibble::tibble(seq = seq,
                                    nreads = nreads),
                     by = "seq")
}

sintax <- function(seq, db, sintax_cutoff = NULL, multithread = FALSE) {
  args <- character()
  if (assertthat::is.string(seq) && file.exists(seq)) {
    seqfile <- seq
    seq <- seqinr::read.fasta(
      file = seqfile,
      as.string = TRUE) %>%
      {tibble(seq = unlist(seqinr::getSequence(., as.string = TRUE)),
              label = seqinr::getName(.))}
  } else {
    seqfile <- tempfile("seq", fileext = "fasta")
    on.exit(file.remove(seqfile))
    if (methods::is(seq, "XStringSet")) {
      Biostrings::writeXStringSet(seq, seqfile)
      seq <- tibble::tibble(seq = as.character(seq, use.names = FALSE),
                            label = names(seq))
    } else if (methods::is(seq, "ShortRead")) {
      ShortRead::writeFasta(seq, seqfile)
      seq <- tibble::tibble(seq = as.character(seq@sread, use.names = FALSE),
                            label = as.character(seq@id, use.names = FALSE))
    } else if (is.character(seq)) {
      if (is.null(names(seq))) names(seq) <- tzara::seqhash(seq)
      is.RNA <- stringr::str_detect(seq[1], "[Uu]")
      if (is.RNA) {
        seq <- Biostrings::RNAStringSet(seq)
      } else {
        seq <- Biostrings::DNAStringSet(seq)
      }
      Biostrings::writeXStringSet(seq, seqfile)
      seq <- tibble::tibble(label = names(seq),
                            seq = as.character(seq, use.names = FALSE))
    }
  }
  args <- c("--sintax", seqfile)
  
  if (assertthat::is.string(db) && file.exists(db)) {
    dbfile <- db
  } else {
    dbfile <- tempfile("db", fileext = "fasta")
    on.exit(file.remove(dbfile))
    
    if (methods::is(db, "XStringSet")) {
      Biostrings::writeXStringSet(db, dbfile)
    } else if (methods::is(db, "ShortRead")) {
      ShortRead::writeFasta(db, dbfile)
    } else if (is.character(db)) {
      if (is.null(names(db))) stop("Taxonomy database given as character vector must be named.")
      is.RNA <- stringr::str_detect(db[1], "[Uu]")
      if (is.RNA) {
        db <- Biostrings::RNAStringSet(db)
      } else {
        db <- Biostrings::DNAStringSet(db)
      }
      Biostrings::writeXStringSet(db, dbfile)
    }
  }
  args <- c(args, "--db", dbfile)
  tablefile <- tempfile("table", fileext = "tsv")
  args <- c(args, "--tabbedout", tablefile)
  on.exit(file.remove(tablefile))
  
  if (!missing(sintax_cutoff)) args <- c(args, "--sintax_cutoff", sintax_cutoff)
  if (!missing(multithread)) {
    if (isFALSE(multithread)) multithread <- 1
    args <- c(args, "--threads", multithread)
  }
  system2("vsearch", args = args)
  # vsearch outputs an extra tab when it cannot place the sequence
  system2("sed", args = c("--in-place", "s/\\t\\t\\t\\t/\\t\\t\\t/", tablefile))
  readr::read_tsv(tablefile,
                  col_names = c("label", "hit", "strand", "c12n"),
                  col_types = "cccc") %>%
    dplyr::left_join(seq, ., by = "label")
}

sintax_format <- function(tax) {
  tax <- dplyr::mutate_at(tax, "c12n", sub, pattern = "(,?[dkpcofgs]:unidentified)+$", replacement = "") %>%
    dplyr::mutate(name = sub(c12n, pattern = ".*,?[dkpcofgs]:", replacement = "")) %>%
    tidyr::separate_rows(c12n, sep = ",") %>%
    dplyr::mutate_at("c12n", dplyr::na_if, "") %>%
    tidyr::separate(c12n, into = c("rank", "taxon"), sep = ":") %>%
    dplyr::mutate_at("taxon", dplyr::na_if, "unidentified") %>%
    dplyr::mutate_at("rank", factor,
                     levels = c("d", "k", "p", "c", "o", "f", "g", "s"),
                     labels = c("domain", "kingdom",
                                "phylum", "class", "order", "family",
                                "genus", "species")) %>%
    tidyr::spread(key = "rank", value = "taxon") %>%
    dplyr::select(-`<NA>`) %>%
    dplyr::mutate_at("name", rlang::`%|%`, "unknown") %>%
    tidyr::unite(col = "Taxonomy",
                 dplyr::one_of("domain", "kingdom", "phylum", "class", "order",
                               "family", "genus", "species"),
                 sep = ";",
                 remove = FALSE) %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "(;?NA)+", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "[|].*", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "_sp^", replacement = "") %>%
    dplyr::mutate_at("Taxonomy", sub, pattern = "([^;]+);\\1", replacement = "\\1")
  
  if ("species" %in% names(tax)) {
    tax <- dplyr::mutate(tax, name = paste0(name, ifelse(is.na(species),
                                                                     "_sp", "")))
  } else {
    tax <- dplyr::mutate(tax, name = paste0(name, "_sp"))
  }
  tax
}
