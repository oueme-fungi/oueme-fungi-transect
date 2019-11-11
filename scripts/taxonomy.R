read_header <- function(header_file, patch_file, format) {
  switch(format,
    rdp = read_header_rdp(header_file, patch_file),
    unite = read_header_unite(header_file, patch_file),
    stop("unknown reference database format: ", format)
  )
}

read_header_rdp <- function(header_file, patch_file) {
  header_file <- if (is.character(header_file)) {
    Biostrings::readDNAStringSet(header_file) 
  } else {
    header_file
  }
  header_file %>%
    names() %>%
    patch_taxonomy(patch_file) %>%
    stringr::str_match("(gi\\|\\d+\\|e?[gm]b\\|)?([A-Z]+_?[0-9]+)[.\\d|-]*[:space:]+Root;(.+)") %>%
    {tibble::tibble(index = seq_len(nrow(.)),
                    accno = .[,3],
                    classifications = .[,4])}
}

read_header_unite <- function(header_file, patch_file) {
  header_file <- if (is.character(header_file)) {
    Biostrings::readDNAStringSet(header_file) 
  } else {
    header_file
  }
  header_file %>%
    names() %>%
    patch_taxonomy(patch_file) %>%
    stringr::str_match("[^|]\\|([A-Z]+_?[0-9]+)\\|[^|]+\\|re[pf]s(_singleton)?\\|(([kpcofgs]__[-\\w.]+;?){7})") %>%
    {tibble::tibble(index = seq_len(nrow(.)),
                    accno = .[,2],
                    classifications = .[,4])} %>%
    dplyr::mutate_at("classifications",
                     stringr::str_replace_all,
                     "[kpcofgs]__", "") %>%
    dplyr::mutate_at("classifications", stringr::str_replace_all,
                     ";unidentified", "")
}

replace_header <- function(in_fasta, out_fasta, new_header,
                           in_format, out_format, 
                           patch_file) {
  fasta <- Biostrings::readDNAStringSet(in_fasta)
  old_header <- read_header(fasta, patch_file, in_format)
  assertthat::assert_that(all(old_header$accno %in% new_header$accno))
  new_header <- dplyr::select(old_header, accno) %>%
    dplyr::left_join(new_header, by = "accno")
  
  write_taxonomy(taxonomy = new_header, fasta = fasta, outfile = out_fasta,
                 format = out_format)
}


accno_c12n_table <- function(tax_data) {
  dplyr::left_join(tibble::tibble(accno = tax_data$data$query_data,
                                  taxon_ids = names(tax_data$data$query_data)),
                   tax_data$get_data_frame(c("taxon_ids", "classifications")),
                   by = "taxon_ids")
}

reduce_ncbi_taxonomy <- function(taxonomy, taxa, ranks, keytaxa) {
  dplyr::filter(taxa, !ncbi_rank %in% ranks,
                !ncbi_name %in% keytaxa) %$%
  dplyr::mutate_at(taxonomy, "classifications",
                   stringi::stri_replace_all_fixed,
                   pattern = paste0(ncbi_name, ";"),
                   replacement = "",
                   vectorize_all = FALSE) %>%
  dplyr::mutate_at("classifications", stringr::str_replace_all, ";$", "")
}

reduce_taxonomy <- function(taxonomy) {
    # take only the first capitalized word at each taxonomic level
    stringr::str_replace_all(taxonomy,
                         "(^|;)[^A-Z;]*([A-Z]+[a-z0-9]+)[^;]*",
                         "\\1\\2") %>%
    # remove repeated taxa (due to removed incertae sedis or species epithet)
    stringr::str_replace_all("([A-Z]+[a-z0-9]+)(;\\1)+(;|$)", "\\1\\3")
}

translate_taxonomy <- function(taxonomy, c12n, reference) {
  c12n <- c12n %>%
    dplyr::mutate(n_supertaxa = stringr::str_count(classifications, ";")) %>%
    dplyr::group_split(n_supertaxa)
  
  for (i in seq_along(c12n)) {
    flog.info("Translating %i taxa at level %i.", nrow(c12n[[i]]), i)
    replacements <- c12n[[i]] %>%
      dplyr::select(taxon_names, classifications) %>%
      dplyr::inner_join(reference, by = "taxon_names",
                        suffix = c("_src", "_ref")) %>%
      dplyr::group_by(classifications_src) %>%
      dplyr::mutate(mismatch = all(classifications_src != classifications_ref)) %>%
      dplyr::filter(mismatch) %>%
      dplyr::mutate(dist = stringdist::stringdist(classifications_src,
                                                  classifications_ref)) %>%
      dplyr::arrange(dist, .by_group = TRUE) %>%
      dplyr::summarize(classifications_ref = dplyr::first(classifications_ref)) %>%
      dplyr::mutate_at("classifications_src", paste0, "(;|$)") %>%
      dplyr::mutate_at("classifications_ref", paste0, "$1")
    if (nrow(replacements)) {
      taxonomy$classifications <-
        stringi::stri_replace_all_regex(taxonomy$classifications,
                                        replacements$classifications_src,
                                        replacements$classifications_ref,
                                        vectorize_all = FALSE)
      for (j in seq_along(c12n)) {
        if (j < i) next
        c12n[[j]]$classifications <-
          stringi::stri_replace_all_regex(c12n[[j]]$classifications,
                                          replacements$classifications_src,
                                          replacements$classifications_ref,
                                          vectorize_all = FALSE)
      }
    }
  }
  flog.info("Uniquifying taxonomy.")
  uniquify_taxonomy(taxonomy, dplyr::bind_rows(c12n))
}

uniquify_taxonomy <- function(taxonomy, c12n) {
  # finding duplicates is much faster than checking whether the
  # classification is present in the taxonomy, so do that first.
  taxdupes <- unique(c12n$taxon_names[duplicated(c12n$taxon_names)])
  flog.info("Found %i initial duplicated taxon names. Searching reference...",
            length(taxdupes))
  c12n <-
    dplyr::filter(c12n, taxon_names %in% taxdupes) %>%
    dplyr::filter(purrr::map_lgl(classifications,
                                 ~any(stringi::stri_detect_fixed(
                                   taxonomy$classifications, .))))
  # now we only need the ones where duplicates are actually present.
  taxdupes <- c12n$taxon_names[duplicated(c12n$taxon_names)]
  flog.info("Found %i duplicated taxon names in reference.", length(taxdupes))
  if (length(taxdupes) > 0) {
    dplyr::filter(c12n, taxon_names %in% taxdupes) %>%
      dplyr::mutate(kingdom = stringr::str_extract(classifications, "^[^;]+")) %>%
      dplyr::filter(kingdom == "Metazoa") %>%
      dplyr::mutate(replacement = paste0(classifications, "(Metazoa)$1"),
                    pattern = paste0(classifications, "(;|$)")) %$%
      dplyr::mutate_at(taxonomy,
                       "classifications",
                       stringi::stri_replace_all_regex,
                       pattern = pattern,
                       replacement = replacement,
                       vectorize_all = FALSE)
  } else {
    taxonomy
  }
  
}

patch_taxonomy <- function(taxonomy, patch_file) {
  assertthat::assert_that((assertthat::is.string(patch_file)
                           && file.exists(patch_file))
                          || is.null(patch_file))
  if (is.null(patch_file)) return(taxonomy)
  
  replace <- readLines(patch_file) %>%
    stringr::str_subset("^s/") %>%
    stringr::str_match("s/(.+)/(.+)/g?")
  if (nrow(replace) == 0) return(taxonomy)
  for (n in seq_len(nrow(replace))) {
    taxonomy <- gsub(replace[n, 2], replace[n,3], taxonomy)
  }
  taxonomy
}

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

# Removes unnecessary ranks from taxonomy, ensures intermediate missing ranks are "incertae sedis" and trailing missing ranks are "unidentified"

regularize_taxonomy <- function(taxonomy, rank_in,
                                rank_out = c("kingdom", "phylum", "class",
                                             "order", "family", "genus"),
                                sep = ";") {
  taxonomy <-
    stringr::str_split_fixed(taxonomy,
                             n = length(rank_in),
                             pattern = sep) %>%
    set_colnames(rank_in) %>%
    tibble::as_tibble() %>%
    dplyr::select(!!!rank_out) %>%
    dplyr::mutate_all(gsub, pattern = ".+[Ii]ncertae[_ ]sedis", replacement = "") %>%
    dplyr::mutate_all(dplyr::na_if, "")
  for (i in 1:(length(rank_out) - 2)) {
    taxonomy[[i + 1]] <- dplyr::coalesce(taxonomy[[i + 1]],
                                         paste0(taxonomy[[i]],
                                                "_Incertae_sedis"))
  }
  taxonomy <-
    dplyr::mutate_all(taxonomy,
                      stringr::str_replace,
                      "(_Incertae_sedis)+",
                      "_Incertae_sedis")
  for (i in (length(rank_out) - 1):2) {
    taxonomy[[i]] <- ifelse(is.na(taxonomy[[i + 1]]) &
                              endsWith(taxonomy[[1]], "_Incertae_sedis"),
                            NA_character_,
                            taxonomy[[i]])
  }
  for (i in 1:(length(rank_out) - 1)) {
    taxonomy[[i + 1]] <- dplyr::coalesce(taxonomy[[i + 1]],
                                         paste("unidentified",
                                               taxonomy[[i]],
                                               sep = "_"))
  }
  taxonomy %>%
    dplyr::mutate_all(stringr::str_replace,
                      "(unidentified_)+",
                      "unidentified_") %>%
    purrr::pmap_chr(paste, sep = ";")
}

read_taxonomy_unite <- function(file) {
  readLines(file) %>%
  stringr::str_subset("^>") %>%
  stringr::str_extract("k__.+") %>%
  stringr::str_replace("s__.+", "") %>%
  stringr::str_replace_all("[kpcofgs]__", "") %>%
  stringr::str_replace("(;unidentified)+$", "") %>%
  unique()
}

read_taxonomy_rdp <- function(file) {
  
}

write_taxonomy <- function(taxonomy, fasta, outfile, format) {
  switch(format,
    dada2 = write_taxonomy_dada2(taxonomy, fasta, outfile),
    sintax = write_taxonomy_sintax(taxonomy, fasta, outfile),
    stop("unknown taxonomy format: ", format)
  )
}

write_taxonomy_dada2 <- function(taxonomy, fasta, outfile) {
  if (is.character(fasta)) {
    fasta <- Biostrings::readDNAStringSet(fasta)
  }
  fasta <- fasta[!is.na(taxonomy$classifications)]
  
  taxonomy %>%
    dplyr::filter(!is.na(classifications)) %$%
    set_names(fasta, classifications) %>%
    Biostrings::writeXStringSet(outfile, compress = endsWith(outfile, ".gz"))
}

write_taxonomy_sintax <- function(taxonomy, fasta, outfile) {
  if (is.character(fasta)) {
    fasta <- Biostrings::readDNAStringSet(fasta)
  }
  fasta <- fasta[!is.na(taxonomy$classifications)]
  
  taxonomy %>%
    dplyr::filter(!is.na(classifications)) %>%
    dplyr::mutate_at("classifications", stringr::str_replace,
                     "([^;]+);([^;]+);([^;]+);([^;]+);([^;]+);([^;]+)",
                     "tax=k:\\1,p:\\2,c:\\3,o:\\4,f:\\5,g:\\6") %$%
    set_names(fasta, classifications) %>%
    Biostrings::writeXStringSet(outfile, compress = endsWith(outfile, ".gz"))
}

read_classification_tedersoo <- function(file) {
    readxl::read_xlsx(file) %>%
    dplyr::select(-subdomain) %>%
    # Remove duplicate taxa within one kingdom
    # Unless noted otherwise the choice of which to remove is based on
    # Index Fungorum (for Fungi) or GBIF (for other organisms)
    dplyr::filter(
      !(genus == "Rhodotorula" & family == "unspecified"),
      !(genus == "Verticillium" & family == "unspecified")#,
      # !(genus == "Automolus" & class == "Insecta"),
      # !(genus == "Clania" & order = "Lepidoptera"),
      # !(genus == "Euxinia" & class == "Malacostraca"),
      # !(genus == "Keijia" & class == "Arachnida"),
      # !(genus == "Napo" & order == "Hymenoptera"),
      # !(genus == "Oxyacanthus" & order == "Amphipoda"),
      # !(genus == "")
    ) %>%
    # dplyr::mutate(
    #   genus = ifelse((genus == "Eutrapela" & order = "Coleoptera"),
    #                  "Chromomoea",
    #                  genus) %>%
    #     ifelse((genus == "Ichthyophaga" & phylum = "Platyhelminthes"),
    #            "Piscinquilinus",
    #            genus)
    # )
    dplyr::mutate_all(dplyr::na_if, "unspecified") %>%
    dplyr::mutate(
      kingdom = dplyr::coalesce(kingdom, 
                                "Eukaryota_reg_Incertae_sedis"),
      subkingdom = dplyr::coalesce(subkingdom, 
                                   paste0(kingdom, "_subreg_Incertae_sedis")),
      phylum = dplyr::coalesce(phylum, 
                               paste0(subkingdom, "_phy_Incertae_sedis")),
      subphylum = dplyr::coalesce(subphylum, 
                                  paste0(phylum, "_subphy_Incertae_sedis")),
      class = dplyr::coalesce(class, 
                              paste0(subphylum, "_cl_Incertae_sedis")),
      order = dplyr::coalesce(order, 
                              paste0(class, "_ord_Incertae_sedis")),
      family = dplyr::coalesce(family, 
                               paste0(order, "_fam_Incertae_sedis"))
    ) %>%
    dplyr::mutate_all(sub,
                      pattern = "(_(sub)?(reg|phy|cl|ord)_Incertae_sedis)+(_(sub)?(reg|phy|cl|ord|fam)_Incertae_sedis$)",
                      replacement = "\\4") %>%
    # The subphylum of Variosea is sometimes missing.
    # The main text says that it should be in Mycetozoa
    dplyr::mutate(subphylum = ifelse(class == "Variosea",
                                     "Mycetozoa",
                                     subphylum),
                  # Craniata is a phylum;
                  # Craniatea is a class in Brachiopoda
                  class = ifelse(class == "Craniata",
                                 "Craniatea",
                                 class)) %>%
    unique() %>%
    taxa::parse_tax_data(class_cols = 1:8,
                         named_by_rank = TRUE) %>%
    taxa::get_data_frame(c("taxon_names", "classifications")) %>%
    dplyr::filter(!stringr::str_detect(taxon_names, "[Ii]ncertae[_ ]sedis"))
  
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
    s <- sample(seq_len(nrow(table)), n)
  }
  return(table[s,,drop = FALSE])
}

formatdb <- function(db = c("unite", "rdp", "silva"),
                     seq_file,
                     patch_file = NULL,
                     maxGroupSize = 10,
                     ...) {
  db = match.arg(db)
  f <- switch(db,
         rdp = formatdb_rdp,
         silva = formatdb_silva,
         unite = formatdb_unite,
         stop("unknown database format: ", db))
  f(seq_file = seq_file, patch_file = patch_file, ...)
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

taxonomy <- function(seq.table, reference, method, multithread = FALSE, ...) {
  # we can take a community matrix (in which case the sequences are the column
  # names) or the sequences
  if (is.matrix(seq.table)) {
    nreads <- Matrix::colSums(seq.table)
    seq <- colnames(seq.table)
  } else {
    nreads <- 1
    seq <- as.character(seq.table)
    names(seq) <- names(seq.table)
  }
  is.RNA <- stringr::str_detect(seq[1], "[Uu]")
  # seq <- if (is.RNA) Biostrings::RNAStringSet(seq) else Biostrings::DNAStringSet(seq)
  
  f <- switch(method,
         dada2 = taxonomy_dada2,
         sintax = taxonomy_sintax,
         idtaxa = taxonomy_idtaxa,
         stop("unknown taxonomy assignment method: ", method)
  )
  f(seq = seq, reference = reference, multithread = multithread, ...)
}

taxonomy_dada2 <- function(seq, reference, multithread = FALSE, threshold = 50, 
                           tryRC = FALSE,
                           outputBootstraps = TRUE,
                           verbose = TRUE, ...) {
  dada2::assignTaxonomy(seqs = chartr("Uu", "Tt", seq),
                        refFasta = reference,
                        tryRC = tryRC,
                        outputBootstraps = outputBootstraps,
                        multithread = multithread,
                        verbose = verbose,
                        minBoot = threshold,
                        ...)
}

taxonomy_sintax <- function(seq, reference, min_confidence = NULL, multithread = FALSE, ...) {
  args <- character()
  if (assertthat::is.string(seq) && file.exists(seq)) {
    seqfile <- seq
    seq <- seqinr::read.fasta(
      file = seqfile,
      as.string = TRUE) %>%
      {tibble(seq = unlist(seqinr::getSequence(., as.string = TRUE)),
              label = seqinr::getName(.))}
  } else {
    seqfile <- tempfile("seq", fileext = ".fasta")
    on.exit(file.remove(seqfile))
    if (methods::is(seq, "XStringSet")) {
      Biostrings::writeXStringSet(seq, seqfile)
      seq <- tibble::tibble(label = names(seq),
                            seq = as.character(seq, use.names = FALSE))
    } else if (methods::is(seq, "ShortRead")) {
      ShortRead::writeFasta(seq, seqfile)
      seq <- tibble::tibble(label = as.character(seq@id, use.names = FALSE),
                            seq = as.character(seq@sread, use.names = FALSE))
    } else if (is.character(seq)) {
      if (is.null(names(seq))) names(seq) <- tzara::seqhash(seq)
      is.RNA <- any(stringr::str_detect(seq, "[Uu]"))
      if (is.RNA) {
        seq <- Biostrings::RNAStringSet(chartr("Tt", "Uu", seq))
      } else {
        seq <- Biostrings::DNAStringSet(seq)
      }
      Biostrings::writeXStringSet(seq, seqfile)
      seq <- tibble::tibble(label = names(seq),
                            seq = as.character(seq, use.names = FALSE))
    }
  }
  args <- c("--sintax", seqfile)
  
  if (assertthat::is.string(reference) && file.exists(reference)) {
    dbfile <- reference
  } else {
    dbfile <- tempfile("db", fileext = ".fasta")
    on.exit(file.remove(dbfile))
    
    if (methods::is(reference, "XStringSet")) {
      Biostrings::writeXStringSet(reference, dbfile)
    } else if (methods::is(reference, "ShortRead")) {
      ShortRead::writeFasta(reference, dbfile)
    } else if (is.character(reference)) {
      if (is.null(names(reference))) stop("Taxonomy database given as character vector must be named.")
      is.RNA <- stringr::str_detect(reference[1], "[Uu]")
      if (is.RNA) {
        reference <- Biostrings::RNAStringSet(reference)
      } else {
        reference <- Biostrings::DNAStringSet(reference)
      }
      Biostrings::writeXStringSet(reference, dbfile)
    }
  }
  args <- c(args, "--db", dbfile)
  tablefile <- tempfile("table", fileext = ".tsv")
  args <- c(args, "--tabbedout", tablefile)
  on.exit(file.remove(tablefile))
  
  if (!missing(min_confidence)) args <- c(args, "--sintax_cutoff", min_confidence)
  if (!missing(multithread)) {
    if (isFALSE(multithread)) multithread <- 1
    args <- c(args, "--threads", multithread)
  }
	  system2("vsearch", args = args)
  if (missing(min_confidence)) {
    # vsearch outputs an extra tab when it cannot place the sequence
    system2("sed", args = c("--in-place", "'s/\\t\\t\\t/\\t\\t/'", tablefile))
    readr::read_tsv(tablefile,
                    col_names = c("label", "hit", "strand"),
                    col_types = "ccc") %>%
      dplyr::left_join(seq, ., by = "label")
  } else {
    # vsearch outputs an extra tab when it cannot place the sequence
    system2("sed", args = c("--in-place", "'s/\\t\\t\\t\\t/\\t\\t\\t/'", tablefile))
    readr::read_tsv(tablefile,
                    col_names = c("label", "hit", "strand", "c12n"),
                    col_types = "cccc") %>%
      dplyr::left_join(seq, ., by = "label")
  }
}

taxonomy_idtaxa <- function(seq, reference, multithread = FALSE, strand = "top", min_confidence = 40, ...) {
  if (isTRUE(multithread)) multithread <- NULL
  if (isFALSE(multithread)) multithread <- 1
  if (!methods::is(seq, "XStringSet")) {
    seq <- Biostrings::DNAStringSet(chartr("Uu", "Tt", seq))
  }
  DECIPHER::IdTaxa(test = seq,
                   trainingSet = reference,
                   strand = strand,
                   processors = multithread,
                   threshold = min_confidence,
                   ...)
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

taxtable <- function(tax, ...) {
  if (methods::is(tax, "Taxa") && methods::is(tax, "Test")) {
    taxtable_idtaxa(tax, ...)
  } else if (is.list(tax) && "tax" %in% names(tax) && "boot" %in% names(tax)) {
    taxtable_dada2(tax, ...)
  } else if (is.data.frame(tax) && all(rlang::has_name(tax, c("label", "hit")))) {
    taxtable_sintax(tax, ...)
  } else {
    stop("Unknown taxonomy table format.")
  }
}

taxtable_sintax <- function(tax, min_confidence = 0, ...) {
  tidyr::separate_rows(tax, hit, sep = ",") %>%
    dplyr::mutate_at("hit", dplyr::na_if, "") %>%
    tidyr::extract(hit, into = c("rank", "taxon", "confidence"),
                   regex = "([dkpcofgs]):([^(]+)\\(([01]\\.\\d+)\\)") %>%
    dplyr::mutate_at("taxon", dplyr::na_if, "unidentified") %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae",
                     replacement = "incertae") %>%
    dplyr::mutate_at("confidence", as.numeric) %>%
    dplyr::mutate_at("rank", factor,
                     levels = c("d", "k", "p", "c", "o", "f", "g", "s"),
                     labels = c("domain", "kingdom",
                                "phylum", "class", "order", "family",
                                "genus", "species")) %>%
    dplyr::select(label, rank, taxon, confidence) %>%
    dplyr::filter(confidence >= min_confidence, !is.na(taxon))
}

taxtable_idtaxa <- function(tax, min_confidence = 0, names = NULL, ...) {
  if (!missing(names) && !is.null(names)) names(tax) <- names
  purrr::imap_dfr(tax, ~tibble::tibble(label = .y,
                                       rank = factor(.x$rank,
                                                     levels = c("rootrank",
                                                                "domain",
                                                                "kingdom",
                                                                "phylum",
                                                                "class",
                                                                "order",
                                                                "family",
                                                                "genus",
                                                                "species")),
                                       taxon = gsub(" ", "_", .x$taxon),
                                       confidence = .x$confidence / 100)) %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae", replacement = "incertae") %>%
    dplyr::filter(rank != "rootrank", !startsWith(taxon, "unclassified_"),
                  confidence >= min_confidence) %>%
    dplyr::mutate_at("rank", factor, levels = c("domain",
                                                "phylum",
                                                "class",
                                                "order",
                                                "family",
                                                "genus",
                                                "species"),
                     labels = c("kingdom",
                                "phylum",
                                "class",
                                "order",
                                "family",
                                "genus",
                                "species"))
}

taxtable_dada2 <- function(tax, names = rownames(tax$tax),
                           min_confidence = 0, ...) {
  taxa <- tax$tax %>% magrittr::set_rownames(names) %>%
    tibble::as_tibble(rownames = "label") %>%
    tidyr::gather(key = "rank", value = "taxon", -1)
  conf <- (tax$boot / 100) %>% magrittr::set_rownames(names) %>%
    tibble::as_tibble(rownames = "label") %>%
    tidyr::gather(key = "rank", value = "confidence", -1)
  dplyr::full_join(taxa, conf, by = c("label", "rank")) %>%
    dplyr::arrange(label) %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae", replacement = "incertae") %>%
    dplyr::filter(!is.na(taxon), confidence >= min_confidence) %>%
    dplyr::mutate_at("rank", tolower)
}

# fits an IDTAXA model to a fasta file with [UV]SEARCH/SINTAX-style taxonomy
# info
fit_idtaxa <- function(fasta) {
  seqdata <- Biostrings::readDNAStringSet(fasta)
  taxonomy <- names(seqdata)
  
  #IDTAXA does not need intermediate/trailing placeholder nodes
  taxonomy <- gsub("[dkpcofgs]:unidentified[^,]+,?", "", taxonomy)
  taxonomy <- gsub("[dkpcofgs]:[^,]+incertae_sedis[^,;]*,?", "", taxonomy)
  
  #it does need a root node
  taxonomy <- gsub(".*tax=", "r:Root,", taxonomy)
  taxonomy <- gsub(";.*", "", taxonomy)
  
  taxdata <- taxa::parse_tax_data(taxonomy,
                                  class_sep = ",",
                                  class_regex = "([rdkpcofgs]):(.+)",
                                  class_key = c("taxon_rank", "taxon_name"))
  edgelist <- taxdata$edge_list
  ranklist <- c("rootrank", "domain", "kingdom", "phylum", "class", "order",
                "family", "genus", "species")
  ranks <- data.frame(Index = taxdata$taxon_indexes()[edgelist$to],
                      Name = taxdata$taxon_names()[edgelist$to],
                      Parent = taxdata$taxon_indexes()[edgelist$from],
                      Level = taxdata$n_supertaxa()[edgelist$to],
                      Rank = match.arg(taxdata$taxon_ranks()[edgelist$to],
                                       ranklist,
                                       several.ok = TRUE),
                      stringsAsFactors = FALSE)
  ranks$Parent[is.na(ranks$Parent)] <- 0
  
  taxonomy <- gsub("[rdkpcofgs]:", "", taxonomy)
  taxonomy <- gsub(",", ";", taxonomy)
  
  trainSet <- DECIPHER::LearnTaxa(train = seqdata,
                                  taxonomy = taxonomy,
                                  rank = ranks,
                                  verbose = TRUE)
  
}
