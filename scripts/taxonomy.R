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

rank_factor <- function(r,
                        ranks = c("rootrank", "domain", "kingdom", "phylum",
                                  "class", "order", "family", "genus",
                                  "species"),
                        abbrev = FALSE) {
  if (abbrev) {
    factor(r, levels = substr(ranks, 1, 1), labels = ranks, ordered = TRUE)
  } else {
    factor(r, levels = ranks, ordered = TRUE)
  }
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
    dplyr::mutate_at("rank", rank_factor, abbrev = TRUE) %>%
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
    dplyr::mutate_at("rank", rank_factor, abbrev = TRUE) %>%
    dplyr::select(label, rank, taxon, confidence) %>%
    dplyr::filter(confidence >= min_confidence, !is.na(taxon))
}

taxtable_idtaxa <- function(tax, min_confidence = 0, names = NULL, ...) {
  if (!missing(names) && !is.null(names)) names(tax) <- names
  purrr::imap_dfr(tax, ~tibble::tibble(label = .y,
                                       rank = rank_factor(.x$rank),
                                       taxon = gsub(" ", "_", .x$taxon),
                                       confidence = .x$confidence / 100)) %>%
    dplyr::mutate_at("taxon", gsub, pattern = "Incertae", replacement = "incertae") %>%
    dplyr::filter(rank != "rootrank", !startsWith(taxon, "unclassified_"),
                  confidence >= min_confidence)
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
    dplyr::mutate_at("rank", tolower) %>%
    dplyr::mutate_at("rank", rank_factor)
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

# combine a named list of taxonomy tables
combine_taxon_tables <- function(tables, allseqs) {
  tibble::enframe(tables) %>%
    tidyr::extract(
      col = name,
      into = c("region", "reference", "ref_region", "method"),
      regex = "taxon_(32S|5_8S|ITS[12]?|LSU|long)_(unite|warcup|rdp_train)_(ITS[12]?|LSU)_(idtaxa|dada2|sintax)",
      remove = FALSE
    ) %>%
    tidyr::unnest("value") %>%
    dplyr::filter(
      !is.na(rank),
      !is.na(label),
      confidence >= 0.5,
      !startsWith(taxon, "unidentified"),
      !endsWith(taxon, "incertae_sedis"),
      reference != "warcup" | rank != "kingdom"
    ) %>%
    dplyr::mutate_at("rank", rank_factor) %>%
    dplyr::group_by(label, rank) %>%
    dplyr::mutate(
      n_tot = dplyr::n(),
      n_diff = dplyr::n_distinct(taxon, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      dplyr::select(allseqs, label = "hash", n_reads = "nreads"),
      by = "label"
    )
}

#### taxon_labels ####
# make labels summarizing the taxonomy of each sequence
make_taxon_labels <- function(t) {
    dplyr::group_by(t, label, rank, n_reads) %>%
    dplyr::summarize(
      taxon = 
        table(taxon) %>%
        paste0(names(.), collapse = "/") %>%
        gsub(pattern = "(.+/.+)", replacement = "<\\1>") %>%
        gsub(pattern = "(mycota|mycetes|ales|aceae)", replacement = "") %>%
        gsub(pattern = "incertae_sedis", replacement = "i_s") %>%
        gsub(pattern = "Fungi\\b", replacement = "F") %>%
        gsub(pattern = "Basidio\\b", replacement = "B") %>%
        gsub(pattern = "Asco\\b", replacement = "A") %>%
        gsub(pattern = "Chytridio\\b", replacement = "Chy") %>%
        gsub(pattern = "Zygo\\b", replacement = "Z")
    ) %>%
    dplyr::group_by(label, n_reads) %>%
    dplyr::arrange(rank) %>%
    dplyr::summarize(tip_label = paste(label[1],
                                       format(n_reads[1], width = 5),
                                       paste0(taxon, collapse = "-")))
}

#### relabel_tree ####
# replaces tree tip labels from old with labels from new
relabel_tree <- function(tree, old, new) {
  tree$tip.label <-
    plyr::mapvalues(tree$tip.label, old, new, warn_missing = FALSE)
  tree
}

# If all members of the clade are assigned uniquely to one taxon, then assign
# that taxon
# Otherwise, if there exists one clade that is at least one of the possibilities
# for all of the members, then assign that one.
clade_taxon <- function(tree, tax, node, rank) {
  tips <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
  taxa <- dplyr::filter(tax, label %in% tips, rank == !!rank) %>%
    dplyr::group_by(label)
  
  # If there are no relevant taxon assignments, then we can't do anything.
  if (nrow(taxa) == 0) return(NA_character_)
  # If only one thing is assigned, then assign that.
  if (dplyr::n_distinct(taxa$taxon) == 1) return(unique(taxa$taxon))
  best_taxon <- unique(taxa$taxon[taxa$n_diff == 1])
  if (length(best_taxon) > 1) return(NA_character_)
  consensus_taxon <- 
    dplyr::group_map(taxa, ~ unique(.$taxon)) %>%
    purrr::reduce(intersect)
  if (length(consensus_taxon) != 1) return(NA_character_)
  consensus_taxon
}


# create an environment to hold the current state of the taxonomic assignment
# between branching recursive calls to phylotax_
new_phylotax_env <- function(tree, taxa, parent = parent.frame()) {
  rlang::child_env(
    .parent = parent,
    node_taxa = tibble::tibble(
      node = integer(),
      rank = taxa$rank[FALSE],
      taxon = character()
    ),
    tip_taxa = dplyr::filter(taxa, label %in% tree$tip.label)
  )
}


phylotax_ <- function(tree, taxa, node, ranks, e) {
  if (length(ranks) == 0) return()
  
  parents <- phangorn::Ancestors(tree, node, type = "all")
  for (rank in ranks) {
    r <- rank_factor(rank)
    if (any(e$node_taxa$node %in% parents & e$node_taxa$rank == rank)) next
    taxon <- clade_taxon(tree, e$tip_taxa, node, r)
    if (is.na(taxon)) {
      futile.logger::flog.debug("Could not assign a %s to node %d.", rank, node)
      for (n in phangorn::Children(tree, node)) {
        if (n > length(tree$tip.label))
          phylotax_(tree, e$tip_taxa, n, ranks, e)
      }
      break
    } else {
      futile.logger::flog.info("Assigned node %d to %s %s.", node, rank, taxon)
      ranks <- ranks[-1]
      e$node_taxa <- dplyr::bind_rows(
        e$node_taxa,
        tibble::tibble(node = node, rank = r, taxon = taxon)
      )
      tips <- tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
      wrongTaxa <- e$tip_taxa %>%
        dplyr::filter(
          label %in% tips,
          rank == r,
          taxon != !!taxon
        ) %>%
        dplyr::select(-rank, -taxon, -n_tot, -n_diff, -n_reads, -confidence)
      # remove assignments which are not consistent with the one we just chose
      e$tip_taxa <- dplyr::bind_rows(
        dplyr::filter(e$tip_taxa, rank < r),
        dplyr::filter(e$tip_taxa, rank >= r) %>%
          dplyr::anti_join(wrongTaxa, by = names(wrongTaxa))
      )
    }
  }
}

# Assign taxon labels to nodes in a tree when there is a consensus of IDs on
# descendent tips.
phylotax <- function(tree, taxa) {
  e <- new_phylotax_env(tree, taxa)
  ranks <- sort(unique(taxa$rank))
  phylotax_(tree, taxa, phangorn::getRoot(tree), ranks, e)
  as.list(e)
}
