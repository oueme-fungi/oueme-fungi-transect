library(drake)
library(magrittr)
loadd(conseq)
conseq %<>% dplyr::filter(!is.na(long)) %>%
  dplyr::mutate(label = tzara::seqhash(long))
lsu_consensus <- dplyr::filter(conseq, !is.na(LSU)) %$%
  set_names(LSU, label)
its_consensus <- dplyr::filter(conseq, !is.na(ITS)) %$%
  set_names(chartr("U", "T", ITS), label)
long_consensus <- conseq %$%
  set_names(long, label)
idtaxaenv <- new.env()
load("reference/rdp_train.idtaxa.Rdata", envir = idtaxaenv)
rdp_train_idtaxa <- idtaxaenv$trainingSet
remove("idtaxaenv")
tax_idtaxa <- DECIPHER::IdTaxa(Biostrings::RNAStringSet(lsu_consensus),
                               trainingSet = rdp_train_idtaxa,
                               strand = "top",
                               threshold = 40,
                               processors = 3)
tictoc::tic()
tax_sintax <- sintax(seq = chartr("U", "T", lsu_consensus),
                     db = "reference/rdp_train.LSU.sintax.fasta.gz",
                     multithread = 3)
tictoc::toc()
tictoc::tic()
tax_dada2 <- dada2::assignTaxonomy(chartr("U", "T", lsu_consensus),
                                   refFasta = "reference/rdp_train.LSU.dada2.fasta.gz",
                                   tryRC = FALSE,
                                   outputBootstraps = TRUE,
                                   verbose = TRUE,
                                   multithread = 3)
tictoc::toc()

idtaxaenv <- new.env()
load("reference/warcup.idtaxa.Rdata", envir = idtaxaenv)
warcup_idtaxa <- idtaxaenv$trainingSet
tax_its_idtaxa <- DECIPHER::IdTaxa(Biostrings::DNAStringSet(its_consensus),
                                   trainingSet = warcup_idtaxa,
                                   strand = "top",
                                   threshold = 40,
                                   processors = 3)
tictoc::tic()
tax_its_sintax <- sintax(seq = its_consensus,
                         db = "reference/warcup.ITS.sintax.fasta.gz",
                         multithread = 3)
tictoc::toc()
tictoc::tic()
tax_its_dada2 <- dada2::assignTaxonomy(its_consensus,
                                       refFasta = "reference/warcup.ITS.dada2.fasta.gz",
                                       tryRC = FALSE,
                                       outputBootstraps = TRUE,
                                       verbose = TRUE,
                                       multithread = 3)
tictoc::toc()

# idtaxaenv <- new.env()
# load("reference/unite.idtaxa.Rdata", envir = idtaxaenv)
# unite_idtaxa <- idtaxaenv$trainingSet
# tax_unite_idtaxa <- DECIPHER::IdTaxa(Biostrings::DNAStringSet(its_consensus),
#                                      trainingSet = unite_idtaxa,
#                                      strand = "top",
#                                      threshold = 40,
#                                      processors = 3)
tictoc::tic()
tax_unite_sintax <- sintax(seq = its_consensus,
                           db = "reference/unite.ITS.sintax.fasta.gz",
                           multithread = 3)
tictoc::toc()

taxmap <- tibble::tibble(name = c("idtaxa_rdp", "sintax_rdp", "dada2_rdp",
                                  "idtaxa_warcup", "sintax_warcup", "dada2_warcup",
                                  "sintax_unite"),
                         tax = list(tax_idtaxa, tax_sintax, tax_dada2,
                                    tax_its_idtaxa, tax_its_sintax, tax_its_dada2,
                                    tax_unite_sintax),
                         seq = list(lsu_consensus, lsu_consensus, lsu_consensus,
                                    its_consensus, its_consensus, its_consensus,
                                    its_consensus),
                         min_confidence = c(0.4, 0.5, 0.5, 0.4, 0.5, 0.5, 0.5)) %>%
  dplyr::mutate(names = lapply(seq, names)) %>%
  dplyr::mutate(tax = purrr::pmap(., taxtable)) %>%
  purrr::pmap(function(tax, name, ...) dplyr::rename_at(tax, c("taxon", "confidence"),
                                                        paste, name, sep = "_")) %>%
  purrr::reduce(dplyr::full_join, by = c("label", "rank")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(n = dplyr::n_distinct(c(taxon_idtaxa_rdp, taxon_sintax_rdp,
                                        taxon_dada2_rdp, taxon_idtaxa_warcup,
                                        taxon_sintax_warcup, taxon_dada2_warcup,
                                        taxon_sintax_unite),
                                      na.rm = TRUE)) %>%
  dplyr::left_join(dplyr::select(conseq, label, nreads), by = "label")

taxmap_labels <- taxmap %>%
  dplyr::filter(rank != "species") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(taxon_consensus = table(c(taxon_idtaxa_rdp, taxon_sintax_rdp,
                                          taxon_dada2_rdp, taxon_idtaxa_warcup,
                                          taxon_sintax_warcup, taxon_dada2_warcup,
                                          taxon_sintax_unite)) %>%
                  paste0(names(.), collapse = "/") %>%
                  gsub(pattern = "(.+/.+)", replacement = "<\\1>") %>%
                  gsub(pattern = "(mycota|mycetes|ales|aceae)", replacement = "") %>%
                  gsub(pattern = "incertae_sedis", replacement = "i_s") %>%
                  gsub(pattern = "Fungi\\b", replacement = "F") %>%
                  gsub(pattern = "Basidio\\b", replacement = "B") %>%
                  gsub(pattern = "Asco\\b", replacement = "A") %>%
                  gsub(pattern = "Chytridio\\b", replacement = "Chy") %>%
                  gsub(pattern = "Zygo\\b", replacement = "Z")) %>%
  dplyr::group_by(label, nreads) %>%
  dplyr::summarize(tip_label = paste(label[1],
                                     format(nreads[1], width = 5),
                                     paste0(taxon_consensus, collapse = "-")))

longtree <- castor::read_tree(file = "RAxML_bipartitions.long",
                              include_node_labels = TRUE) %>%
  inset2("tip.label",
         plyr::mapvalues(.$tip.label, 
                         taxmap_labels$label,
                         taxmap_labels$tip_label,
                         warn_missing = FALSE))

castor::write_tree(longtree, "output/labeled_long.tree")

lsutree <- castor::read_tree(file = "RAxML_bipartitions.lsu",
                              include_node_labels = TRUE) %>%
  inset2("tip.label",
         plyr::mapvalues(.$tip.label, 
                         taxmap_labels$label,
                         taxmap_labels$tip_label,
                         warn_missing = FALSE))

castor::write_tree(lsutree, "output/labeled_lsu.tree")

taxon_consensus

dplyr::filter(taxmap, n > 1) %>% View


warcup <- readLines("reference/warcup.fasta.gz") %>%
  stringr::str_subset("^>") %>%
  stringr::str_extract("Fungi;.+") %>%
  stringr::str_replace(";[^;]+;?$", "") %>%
  unique()

warcup_replace <- readLines("reference/warcup.sed") %>%
  stringr::str_subset("^s/") %>%
  stringr::str_match("s/(.+)/(.+)/g?")
for (n in 1:nrow(warcup_replace)) {
  warcup <- gsub(warcup_replace[n, 2], warcup_replace[n,3], warcup)
}

warcup_taxa <- taxa::parse_tax_data(warcup)
warcup_class <- warcup_taxa$get_data_frame(c("taxon_names", "classifications")) %>%
  dplyr::filter(!stringr::str_detect(taxon_names, "[Ii]ncertae[_ ]sedis"))

unite <- readLines("reference/unite.fasta.gz") %>%
  stringr::str_subset("^>") %>%
  stringr::str_extract("k__.+") %>%
  stringr::str_replace("s__.+", "") %>%
  stringr::str_replace_all("[kpcofgs]__", "") %>%
  stringr::str_replace("(;unidentified)+$", "") %>%
  unique()

unite_replace <- readLines("reference/unite.sed") %>%
  stringr::str_subset("^s/") %>%
  stringr::str_match("s/(.+)/(.+)/g?")
for (n in 1:nrow(unite_replace)) {
  unite <- gsub(unite_replace[n, 2], unite_replace[n,3], unite)
}

unite_taxa <- taxa::parse_tax_data(unite)
unite_class <- unite_taxa$get_data_frame(c("taxon_names", "classifications")) %>%
  dplyr::filter(!stringr::str_detect(taxon_names, "[Ii]ncertae[_ ]sedis"),
                !stringr::str_detect(taxon_names, "unidentified"))

warcup_mismatches <- dplyr::left_join(warcup_class, unite_class,
                                      by = "taxon_names", suffix = c("_warcup", "_unite")) %>%
  dplyr::filter(classifications_warcup != classifications_unite) %>%
  dplyr::mutate(root1_warcup = stringr::str_extract(classifications_warcup, "^[^;]+;"),
                root2_warcup = stringr::str_extract(classifications_warcup, "[^;]+;[^;]+;"),
                leaf1_warcup = stringr::str_extract(classifications_warcup, "[^;]+$"),
                leaf2_warcup = stringr::str_extract(classifications_warcup, "[^;]+;[^;]+$"),
                root1_unite = stringr::str_extract(classifications_unite, "^[^;]+;"),
                root2_unite = stringr::str_extract(classifications_unite, "[^;]+;[^;]+;"),
                leaf1_unite = stringr::str_extract(classifications_unite, "[^;]+$"),
                leaf2_unite = stringr::str_extract(classifications_unite, "[^;]+;[^;]+$"),
                rootmatch = (root2_unite == root2_warcup) %>% tidyr::replace_na(FALSE),
                leafmatch = (leaf2_unite == leaf2_warcup)  %>% tidyr::replace_na(FALSE)) %>%
  dplyr::select(-taxon_names)

start <- TRUE
while (start || any(warcup_mismatches$rootmatch) || any(warcup_mismatches$leafmatch)) {
  start <- FALSE
  warcup_mismatches %<>%
    dplyr::mutate(
      classifications_warcup = dplyr::if_else(
        rootmatch,
        stringr::str_replace(classifications_warcup, "^[^;]+;", ""),
        classifications_warcup),
      classifications_warcup = dplyr::if_else(
        leafmatch,
        stringr::str_replace(classifications_warcup, ";[^;]+$", ""),
        classifications_warcup),
      classifications_unite = dplyr::if_else(
        rootmatch,
        stringr::str_replace(classifications_unite, "^[^;]+;", ""),
        classifications_unite),
      classifications_unite = dplyr::if_else(
        leafmatch,
        stringr::str_replace(classifications_unite, ";[^;]+$", ""),
        classifications_unite),
      root1_warcup = stringr::str_extract(classifications_warcup, "^[^;]+;"),
      root2_warcup = stringr::str_extract(classifications_warcup, "[^;]+;[^;]+;"),
      leaf1_warcup = stringr::str_extract(classifications_warcup, "[^;]+$"),
      leaf2_warcup = stringr::str_extract(classifications_warcup, "[^;]+;[^;]+$"),
      root1_unite = stringr::str_extract(classifications_unite, "^[^;]+;"),
      root2_unite = stringr::str_extract(classifications_unite, "[^;]+;[^;]+;"),
      leaf1_unite = stringr::str_extract(classifications_unite, "[^;]+$"),
      leaf2_unite = stringr::str_extract(classifications_unite, "[^;]+;[^;]+$"),
      rootmatch = (root2_unite == root2_warcup) %>% tidyr::replace_na(FALSE),
      leafmatch = (leaf2_unite == leaf2_warcup) %>% tidyr::replace_na(FALSE)
    ) %>%
    unique() 
}

incertae_entries <- gsub("(;[^;]*[Ii]ncertae[ _]sedis)+", "\\1", warcup) %>%
  stringr::str_extract("Fungi(;[^;]+)*;[^;]+[Ii]ncertae[ _]sedis;[^;]+") %>%
  unique() %>%
  stringr::str_subset("(.+mycetes);\\1[_ ][Ii]ncertae[ _]sedis;.+ales", negate = TRUE)

incertae_terminal <- stringr::str_extract(incertae_entries, "[^;]+$")
incertae_parent <- stringr::str_match(incertae_entries, "^(.+?);[^;]+_Incertae sedis")[,2]
library(rgbif)
rgbif_terminal <- purrr::map_dfr(incertae_terminal, rgbif::name_backbone)
incertae_replacements <-
  rgbif_terminal  %>%
  dplyr::mutate(genus = dplyr::coalesce(genus,
                                        paste0(dplyr::coalesce(family, order,
                                                               class, phylum,
                                                               kingdom),
                                               "_Incertae sedis")),
                family = dplyr::coalesce(family,
                                         paste0(dplyr::coalesce(order, class,
                                                                phylum, kingdom),
                                                "_Incertae sedis")),
                order = dplyr::coalesce(order,
                                        paste0(dplyr::coalesce(class,
                                                               phylum, kingdom),
                                               "_Incertae sedis")),
                class = dplyr::coalesce(class,
                                        paste0(dplyr::coalesce(phylum, kingdom),
                                               "_Incertae sedis")),
                phylum = dplyr::coalesce(phylum,
                                         paste0(kingdom,"_Incertae sedis"))) %>%
  with(paste(kingdom, phylum, class, order, family, genus, sep = ";"))
incertae_patterns <- paste0(incertae_parent, ";([^;]+_Incertae sedis;)+", incertae_terminal)

incertae_table <- tibble::tibble(entry = incertae_entries,
                                 parent = incertae_parent,
                                 terminal = incertae_terminal,
                                 pattern = incertae_patterns,
                                 replacement = incertae_replacements) %>%
  dplyr::filter(!is.na(entry), !stringr::str_detect(replacement, pattern)) %>%
  dplyr::mutate(pattern2 = paste0(parent, ";.+;", terminal),
                replacement2 = stringr::str_extract(replacement, pattern2))
stringr::str_detect(incertae_replacements, incertae_patterns)

rdp <- readLines("reference/rdp_train.fasta.gz") %>%
  stringr::str_subset("^>")

rdp_nonfungi <- stringr::str_subset(rdp, "Fungi", negate = TRUE) %>%
  stringr::str_match(">([^.]+)(\\.\\d+\\.\\d+)?[:space:]+Root;.+") %>%
  magrittr::extract(,2)
nonfungi_summ <- split(rdp_nonfungi, seq_along(rdp_nonfungi) %/% 50) %>%
  purrr::map(rentrez::entrez_summary, db = "nucleotide")

nonfungi_taxids <- purrr::map(nonfungi_summ, rentrez::extract_from_esummary, "taxid") 
taxize::use_entrez()
nonfungi_class <- purrr::map(nonfungi_taxids, ~ { Sys.sleep(1); taxize::classification(., db = "ncbi")})
rdp_nonfungi_taxdata <- list()
nonfungi_split <- split(rdp_nonfungi, seq_along(rdp_nonfungi) %/% 50)
for (i in seq_along(nonfungi_split)) {
  if (i > length(rdp_nonfungi_taxdata)) {
    rdp_nonfungi_taxdata[[i]] <- taxa::lookup_tax_data(nonfungi_split[[i]], type = "seq_id")
  }
}

ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")

# rdp_nonfungi_newheader <- purrr::map_dfr(rdp_nonfungi_taxdata, function(d) {
#   dplyr::left_join(tibble::tibble(accno = d$data$query_data,
#                                   taxon_ids = names(d$data$query_data)),
#                    d$get_data_frame(c("taxon_ids", "classifications")),
#                    by = "taxon_ids")
# }) %>%
#   dplyr::mutate_at("classifications", paste0, ";")

rdp_nonfungi_newheader <- dplyr::select(ncbi_taxmap, accno, classifications) %>%
  dplyr::mutate_at("classifications", stringr::str_replace, ".*Eukaryota;", "")

rdp_nonfungi_taxa <- purrr::map_dfr(rdp_nonfungi_taxdata, ~ .$data$tax_data) %>%
  unique()

rdp_nonfungi_reducedheader <-
  dplyr::filter(rdp_nonfungi_taxa, !ncbi_rank %in% ranks,
                !ncbi_name %in% unite_class$taxon_names) %$%
  dplyr::mutate_at(rdp_nonfungi_newheader, "classifications",
                   stringi::stri_replace_all_fixed,
                   pattern = paste0(ncbi_name, ";"),
                   replacement = "",
                   vectorize_all = FALSE) %>%
  dplyr::mutate_at("classifications", stringr::str_replace_all, ";$", "")

rdp_nonfungi_class <- rdp_nonfungi_newheader$classifications %>%
  unique() %>%
  taxa::parse_tax_data() %>%
  taxa::get_data_frame(c("taxon_names", "classifications")) %>% 
  dplyr::mutate(n_supertaxa = stringr::str_count(classifications, ";")) %>%
  dplyr::group_split(n_supertaxa)

# rdp_nonfungi_class <- purrr::map_dfr(rdp_nonfungi_taxdata, function(d) {
#   taxa::filter_taxa(d, ncbi_rank %in% ranks |
#                       taxon_names %in% unite_class$taxon_names,
#                     reassign_obs = FALSE)$get_data_frame(c("taxon_names", "classifications"))
# }) %>%
#   unique() %>% 
#   dplyr::mutate(n_supertaxa = stringr::str_count(classifications, ";")) %>%
#   dplyr::group_split(n_supertaxa)

purrr::map_dfr(rdp_nonfungi_taxdata, taxa::get_data_frame, c("taxon_names", "ncbi_rank", "classifications")) %>%
  dplyr::filter(stringr::str_detect(rdp_nonfungi_class[[21]]$classifications, taxon_names)) %>%
  unique() %>%
  View

for (i in seq_along(rdp_nonfungi_class)) {
  # if (i <= 3) {
    replacements <- rdp_nonfungi_class[[i]] %>%
      dplyr::select(taxon_names, classifications) %>%
      dplyr::inner_join(tedersoo_class, by = "taxon_names",
                        suffix = c("_rdp", "_tedersoo")) %>%
      dplyr::group_by(classifications_rdp) %>%
      dplyr::mutate(mismatch = all(classifications_rdp != classifications_tedersoo)) %>%
      dplyr::filter(mismatch) %>%
      dplyr::mutate(dist = stringdist::stringdist(classifications_rdp, classifications_tedersoo)) %>%
      dplyr::arrange(dist, .by_group = TRUE) %>%
      dplyr::summarize(classifications_tedersoo = dplyr::first(classifications_tedersoo))
    if (nrow(replacements)) {
      rdp_nonfungi_reducedheader$classifications <-
        stringi::stri_replace_all_fixed(rdp_nonfungi_reducedheader$classifications,
                                        replacements$classifications_rdp,
                                        replacements$classifications_tedersoo,
                                        vectorize_all = FALSE)
      for (j in seq_along(rdp_nonfungi_class)) {
        if (j < i) next
        rdp_nonfungi_class[[j]]$classifications <-
          stringi::stri_replace_all_fixed(rdp_nonfungi_class[[j]]$classifications,
                                          replacements$classifications_rdp,
                                          replacements$classifications_tedersoo,
                                          vectorize_all = FALSE)
      }
    }
  # } else {
  #   kingdoms <- rdp_nonfungi_reducedheader$classifications %>%
  #     stringr::str_extract("^[^;]+") %>%
  #     unique()
  #   for (k in kingdoms) {
  #     replacements <- rdp_nonfungi_class[[i]] %>%
  #       dplyr::select(taxon_names, classifications) %>%
  #       dplyr::filter(startsWith(classifications, k)) %>%
  #       dplyr::inner_join(tedersoo_class %>%
  #                           dplyr::filter(startsWith(classifications, k)),
  #                         by = "taxon_names",
  #                         suffix = c("_rdp", "_tedersoo")) %>%
  #       dplyr::group_by(classifications_rdp) %>%
  #       dplyr::mutate(mismatch = all(classifications_rdp != classifications_tedersoo)) %>%
  #       dplyr::filter(mismatch) %>%
  #       dplyr::mutate(dist = stringdist::stringdist(classifications_rdp, classifications_tedersoo)) %>%
  #       dplyr::arrange(dist, .by_group = TRUE) %>%
  #       dplyr::summarize(classifications_tedersoo = dplyr::first(classifications_tedersoo))
  #     if (nrow(replacements)) {
  #       ki <- which(startsWith(rdp_nonfungi_reducedheader$classifications, k))
  #       rdp_nonfungi_reducedheader$classifications[ki] %<>%
  #         stringi::stri_replace_all_fixed(replacements$classifications_rdp,
  #                                         replacements$classifications_tedersoo,
  #                                         vectorize_all = FALSE)
  #       for (j in seq_along(rdp_nonfungi_class)) {
  #         if (j < i) next
  #         ki <- which(startsWith(rdp_nonfungi_class[[j]]$classifications, k))
  #         rdp_nonfungi_class[[j]]$classifications[ki] %<>%
  #           stringi::stri_replace_all_fixed(replacements$classifications_rdp,
  #                                           replacements$classifications_tedersoo,
  #                                           vectorize_all = FALSE)
  #       }
  #     }
  #   }
  # }
}

rdp_nonfungi_reducedheader %<>%
  stringr::str_replace("Apicomplexa")

dplyr::filter(rdp_nonfungi_reducedheader, stringr::str_count(classifications, ";") < 7) %>% View
stringr::str_count(rdp_nonfungi_reducedheader$classifications, ";") %>% table

all(rdp_nonfungi_reducedheader$classifications %in% tedersoo_class$classifications)

dplyr::filter(rdp_nonfungi_reducedheader,
              !classifications %in% tedersoo_class$classifications) %>%
  dplyr::left_join(rdp_nonfungi_newheader, by = "accno") %>%
  dplyr::select(classifications.x, classifications.y) %>%
  unique() %>%
  View()

rdp_mismatches[[1]]

nonfungi_test_classifications <- nonfungi_test$get_data_frame(c("taxon_names", "classifications"))


tedersoo_raw <-
  readxl::read_xlsx("reference/Tedersoo_Eukarya_classification.xlsx")

tedersoo_taxmap <- tedersoo_raw %>%
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
                     named_by_rank = TRUE)

tedersoo_class <- tedersoo_taxmap %>%
  taxa::get_data_frame(c("taxon_names", "classifications")) %>%
  dplyr::filter(!stringr::str_detect(taxon_names, "[Ii]ncertae[_ ]sedis"))

dplyr::group_by(tedersoo_raw, genus) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  View

rdp_nonfungi %>%
  taxa::lookup_tax_data(type = "seq_id", database = "ncbi")
rdp_nonfungi_summ <- rentrez::entrez_summary("nucleotide", rdp_nonfungi)
testtax <- rentrez::entrez_summary("taxonomy", "99159")
testtax2 <- rentrez::entrez_
#testing:
fasta <- "reference/rdp_train.LSU.sintax.fasta.gz"
