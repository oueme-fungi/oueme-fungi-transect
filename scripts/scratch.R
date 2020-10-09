library(drake)
library(dplyr)
library(tidyverse)
library(magrittr)
library(Biostrings)


seq_file <- file.path("reference", "rdp_train.LSU.fasta.gz")

rdp_train_names <- Biostrings::readDNAStringSet(seq_file) %>%
  names()
tax <-
  rdp_train_names %>%
  taxa::extract_tax_data(key = c(c12n = "class"),
                         regex = "[^\\t ]+[\\t ]+(.+)",
                         class_sep = ";",
                         )

library(taxa)

sessionInfo()

taxdata <-
  c("Root;Fungi;Basidiomycota;Agaricomycetes;Hymenochaetales;Hymenochaetaceae;Fuscoporia",
    "Root;Fungi;Basidiomycota;Microbotryomycetes;Microbotryales;Microbotryaceae;Microbotryum",
    "Root;Fungi;Ascomycota;Dothideomycetes;Botryosphaeriales;Botryosphaeriaceae;Microdiplodia")

tax <- parse_tax_data(taxdata, class_sep = ";")

ranks <- lapply(c("rootrank", "kingdom", "phylum", "class", "order", "family", "genus"),
                taxon_rank)

rank_idx <- tax$n_supertaxa() + 1

for (i in seq_along(rank_idx)) {
  tax$taxa[[i]]$rank <- ranks[[rank_idx[i]]]
}

tax$taxa



tax_file <- file.path("reference", "rdp_train_tax.txt")

rdp_train_taxonomy2 <-
  readr::read_delim(tax_file, delim = "*",
                  col_names = c("ID", "Name", "Parent",
                                "Level", "Rank")) %>%
  taxa::parse_edge_list(taxon_id = "ID",
                        supertaxon_id = "Parent",
                        taxon_name = "Name",
                        taxon_rank = "Rank")

rdp_train_taxonomy$taxon_ranks() <-
  str_count(rdp_train_taxonomy$classifications(), ";") %>%
  factor(labels = c("rootrank", "kingdom", "phylum", "class", "order", "family", "genus"))

rdp_train_taxonomy %>%
  filter_taxa(str_detect(taxon_names, "incertae"), invert = TRUE)


patch_file <- file.path("reference", "rdp_train_patch.csv")
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

# define the ranks.  The aim is only to search to genus level.
rdp_ranks <-  c("rootrank", "domain", "phylum", "class", "order",
                "family", "genus", "species")

rdp_train_taxa %<>%
  dplyr::mutate_at("c12n", stringr::str_replace_all,"NA;", "") %>%
  dplyr::filter(!stringr::str_detect(Name, "incertae")) %>%
  dplyr::mutate_at("c12n", stringr::str_replace_all, "[A-Z][a-z]+ incertae sedis;", "") %>%
  mutate_at("Rank", factor, levels = rdp_ranks) %>%
  mutate_at("Rank", recode, domain = "kingdom")
DECIPHER::AlignSeqs()

seq_file = file.path("reference", "unite.fasta.gz")
patch_file = file.path("reference", "unite_patch.csv")
unite_taxonomy <- Biostrings::readDNAStringSet(seq_file) %>% unique() %>%
  names() %>%
  taxa::extract_tax_data(regex = ".+\\|([^|]+)",
                         key = "class",
                         class_regex = "([kpcofgs])__(.+)",
                         class_key = c("taxon_rank", "taxon_name"),
                         class_sep = ";")

unite_ranks <- c("rootrank", "kingdom", "phylum", "class", "order", "family", "genus", "species")
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

unite_taxa %<>%

filecombine <- left_join(select(rdp_train_taxa, Name, Parent = ParentName,
                                Rank = Rank) %>%
                           mutate_at("Rank", factor, levels = rdp_ranks) %>%
                           mutate_at("Rank", recode, domain = "kingdom") %>%
                           filter(Rank != "species") %>%
                           unique() %>%
                           left_join(select(., Parent = Name, ParentRank = Rank)),
                         select(unite_taxa, Name, Parent = parent,
                                Rank = Rank) %>%
                           filter(Rank != "species") %>%
                           unique() %>%
                           left_join(select(., Parent = Name, ParentRank = Rank)),
                         by = "Name",
                         suffix = c(".rdp", ".unite")) %>%
  filter(is.na(Parent.unite) | !(Parent.rdp == Parent.unite & Rank.rdp == Rank.unite)) %>%
  anti_join(select(unite_taxa, Parent.rdp = Name, ParentRank.rdp = Rank))

alignedrdp <- readRNAMultipleAlignment("reference/rdp.LSU.st.gz", format = "stockholm")

allcached <- drake::cached()
loadd("cons_tax")
loadd("combined_pb_500")

loadd("classifier_rdp_LSU")
loadd("classifier_silva_LSU")
loadd("classifier_unite_ITS2")

silvatax <- as_tibble(classifier_silva_LSU[c("taxonomy", "taxa", "ranks")])
rdptax <- as_tibble(classifier_rdp_LSU[c("taxonomy", "taxa", "ranks")])
unitetax <- as_tibble(classifier_unite_ITS2[c("taxonomy", "taxa", "ranks")])
loadd("conseq")

loadd("dbprep_unite_ITS2")
loadd("dbprep_unite_ITS")
loadd("dbprep_silva_LSU")

posttax <- function(tax, threshhold = 0.6) {
  tax$confidence <- tax$confidence/tax$confidence[1]
  tax <- tibble::as_tibble(tax)
  tax <- dplyr::filter(tax, confidence >= threshhold)
  tax <- dplyr::slice(tax, -1)
  tax <- dplyr::select(tax, rank, taxon)
  tax <- tidyr::spread(tax, key = rank, value = taxon)
}
testLSU <- Biostrings::RNAStringSet(dplyr::filter(conseq, complete.cases(conseq))$LSU[1:10 + 100])
testITS <- Biostrings::RNAStringSet(dplyr::filter(conseq, complete.cases(conseq))$ITS[1:10 + 100])
testITS2 <- Biostrings::DNAStringSet(dplyr::filter(conseq, complete.cases(conseq))$ITS2[1:10 + 100])

nrow(conseq)


tax_silva <- DECIPHER::IdTaxa(testLSU, classifier_silva_LSU, strand = "top", processors = 3, threshold = 30)
tax_silva2 <- DECIPHER::IdTaxa(testLSU, classifier_silva_LSU, strand = "top", processors = 3, threshold = 30, samples = 2 * L^0.47)
tax_silva3 <- DECIPHER::IdTaxa(testLSU, classifier_silva_LSU, strand = "top", processors = 3, threshold = 30, samples = 0.5 * L^0.47)
dbprep_silva_LSU$train %>%
  set_names(dbprep_silva_LSU$taxonomy)
dadatax_silva <- dada2::assignTaxonomy(testLSU, )
tax_rdp <- DECIPHER::IdTaxa(testLSU, classifier_rdp_LSU, strand = "top", processors = 3, threshold = 50)
tax_unite <- DECIPHER::IdTaxa(testITS2, classifier_unite_ITS2, strand = "top", processors = 3, threshold = 30)

tax_rdp2 <- dada2::assignTaxonomy(as.character(Biostrings::DNAStringSet(testLSU)), "reference/rdp.LSU.dada2.fasta.gz")
tax_unite2 <- dada2::assignTaxonomy(as.character(Biostrings::DNAStringSet(testITS)), "reference/unite.ITS2.fasta.gz")

tax_dada_unite <- dada2::assignTaxonomy(testITS, unite_file)

inocybe_rdp <- dplyr::filter(rdp_c12n, grepl("Inocybe", c12n))
unique(inocybe_rdp$c12n)

penicillium_rdp <- dplyr::filter(rdp_c12n, grepl("Penicillium", c12n))
unique(penicillium_rdp$species)

silva_db[c(250, 389,9594, 10092,10368, 14830, 14837,15414, 15501, 21447, 22943)]
testseq <- ape::as.DNAbin(Biostrings::DNAStringSet(silva_db[c(1291,
                                                              107650,
                                                              108150,
                                                              110387,
                                                              113291,
                                                              113331,
                                                              113356,
                                                              114795,
                                                              115896,
                                                              115972,
                                                              115990

)]))
testkmers <- kmer::kcount(testseq)
testkmeans <- kmeans(testkmers, 10)
dplyr::n_distinct(testkmeans$cluster)

aga_long <- combined_pb_500 %>%
  dplyr::group_by(ITS2) %>%
  dplyr::mutate(nreads = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nreads == max(nreads), !is.na(long)) %$%
  rlang::set_names(long, seqhash(long)) %>%
  unique() %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet()

aga_aln <- DECIPHER::AlignSeqs(aga_long, iterations = 10, refinements = 10)

aga_phydat <- phangorn::phyDat(ape::as.DNAbin(Biostrings::DNAStringSet(aga_aln)))
aga_dist <- phangorn::dist.hamming(aga_phydat)
aga_upgma <- phangorn::NJ(aga_dist)
aga_mp <- phangorn::optim.parsimony(aga_upgma, aga_phydat)
aga_boot <- phangorn::bootstrap.phyDat(aga_phydat, phangorn::pratchet, bs = 100)
optim

phyloseq::otu_table()

groups <- its_join(combined_pb_500)
cons_chim <- cons_tax %>%
  dplyr::select(sequence = long, abundance = nreads) %>%
  dada2::isBimeraDenovo()

qlong <- Biostrings::RNAStringSet(cons_tax$long) %>%
  Biostrings::letterFrequency(., "MRWSYKVHDBN")
table(qlong)
cons_tax$qlong <- qlong

maxabund_ambig <- cons_tax %>% dplyr::filter(qlong > 1) %>% dplyr::arrange(dplyr::desc(nreads)) %$% ITS2[1]
worst_ambig <- cons_tax %>% dplyr::arrange(dplyr::desc(as.integer(qlong))) %$% ITS2[1]

dplyr::filter(combined_pb_500, ITS2 == maxabund_ambig, !is.na(long)) %$%
  rlang::set_names(long, seqhash(long)) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet() %>%
  Biostrings::writeXStringSet("data/maxabund.fasta")

dplyr::filter(combined_pb_500, ITS2 == worst_ambig, !is.na(long)) %$%
  rlang::set_names(long, seqhash(long)) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet() %>%
  Biostrings::writeXStringSet("data/worst.fasta")


bigmaps <- drake::cached() %>%
  str_subset("dada_map_pb_500") %>%
  map(readd, character_only = TRUE)



pbraw <- drake::cached() %>%
  str_subset("^regions_pb_500") %>%
  map_dfr(function(name) {
    x = readd(name, character_only = TRUE)
    if (!is(x, "ShortReadQ")) return(tibble())
    tibble(file = name,
           data = list(tibble(
             seq.id = as.character(x@id),
             seq = as.character(x@sread),
             ee = rowSums(10^-(as(x@quality, "matrix")/10), na.rm = TRUE))))
  }) %>%
  tidyr::extract(
    file,
    c("Seq.Run", "Plate", "Well", "Direction", "Region"),
    "regions_([:alpha:]{2}_\\d{3})(\\d{3})([A-H]1?\\d)([rf]?)([:alnum:]+)") %>%
  unnest(data)

hiQ_LSU <- pbraw %>%
  dplyr::filter(Region == "LSU",
                ee <= 3) %$%
  rlang::set_names(seq, seq.id) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet()

hiQ_LSU_aln <- DECIPHER::AlignSeqs(hiQ_LSU, processors = 3)

combine_bigmaps <- function(bigmaps, rawdata) {
  purrr::map_dfr(bigmaps, ~tibble::tibble(file = names(.), data = .)) %>%
    tidyr::extract(
      col = "file",
      into = c("Seq.Run", "Plate", "Well", "Direction", "Region"),
      regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.qfilt\\.fastq\\.gz") %>%
    tidyr::unnest(data) %>%
    dplyr::full_join(rawdata) %>%
    dplyr::group_by(seq.id) %>%
    dplyr::filter(any(!is.na(asv.idx))) %>%
    dplyr::mutate(seq = dplyr::coalesce(asv.seq, derep.seq, seq)) %>%
    dplyr::select(-derep.seq, -derep.idx, -asv.seq, -asv.idx) %>%
    tidyr::spread(key = Region, value = seq)
}

combined <- combine_bigmaps(bigmaps, pbraw)

chimeras <- list()
chimeras[["ITS"]] <- table(combined$ITS) %>%
  tibble::tibble(sequence = names(.), abundance = .) %>%
  dada2::isBimeraDenovo(multithread = 3)

chimeras[["LSU"]] <- table(combined$LSU) %>%
  tibble::tibble(sequence = names(.), abundance = .) %>%
  dada2::isBimeraDenovo(multithread = 3)

pen_testlist <-
  combined %>% group_by(ITS2) %>%
  mutate(nreads = n()) %>%
  ungroup() %>%
  filter(nreads == max(nreads)) %$%
  rlang::set_names(LSU, seq.id) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet()
aga_testlist <- combined %>% group_by(ITS2) %>%
  mutate(nreads = n()) %>%
  dplyr::ungroup() %>%
  filter(nreads == 683) %$%
  rlang::set_names(LSU, seq.id) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet()
russ_testlist <- combined %>% group_by(ITS2) %>%
  mutate(nreads = n()) %>%
  dplyr::ungroup() %>%
  filter(nreads == 342) %$%
  rlang::set_names(LSU, seq.id) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet()

pen_dmat <- suppressWarnings(DECIPHER::DistanceMatrix(pen_testlist, verbose = FALSE))
outliers1 <- odseq::odseq_unaligned(dmat, type = "distance")

pen_aln <- DECIPHER::AlignSeqs(pen_testlist, processors = 3)
pen_dmat <- DECIPHER::DistanceMatrix(pen_aln, correction = "Jukes-Cantor")
pen_tree <- DECIPHER::IdClusters(pen_dmat, method = "ML", myXStringSet = pen_aln, processors = 3, type = "dendrogram")
DECIPHER::WriteDendrogram(pen_tree, file = "Penicillium.tree")

russ_aln <- DECIPHER::AlignSeqs(russ_testlist, processors = 3)
russ_dmat <- DECIPHER::DistanceMatrix(russ_aln, correction = "Jukes-Cantor")
russ_tree <- DECIPHER::IdClusters(russ_dmat, method = "ML", myXStringSet = russ_aln, processors = 3, type = "dendrogram")
DECIPHER::WriteDendrogram(russ_tree, file = "Russula.tree")

aga_aln <- DECIPHER::AlignSeqs(aga_testlist, processors = 3)
aga_dmat <- DECIPHER::DistanceMatrix(aga_aln, correction = "Jukes-Cantor")
aga_tree <- DECIPHER::IdClusters(aga_dmat, method = "ML", myXStringSet = aga_aln, processors = 3, type = "dendrogram")
DECIPHER::WriteDendrogram(aga_tree, file = "Agaricus.tree")



aln <- aln %>%
  Biostrings::RNAMultipleAlignment() %>%
  Biostrings::maskGaps(min.fraction = 0.5, min.block.width = 1) %>%
  as("RNAStringSet")
dmat <- DECIPHER::DistanceMatrix(aln, correction = "Jukes-Cantor")
testtree <- DECIPHER::IdClusters(dmat, method = "NJ", showPlot = TRUE, myXStringSet = aln, processors = 3, type = "dendrogram")
DECIPHER::WriteDendrogram(testtree, file = "Penicillium.tree")

outliers <- odseq::odseq(RNAMultipleAlignment(aln))
c(aln[!outliers], aln[outliers]) %>% writeXStringSet("Russula-outliers.fasta")
summary(c(dmat))


calculate_consensus2 <- function(seq, names) {
  seq <- rlang::set_names(seq, names)
  seq <- stats::na.omit(seq)
  if (length(seq) < 3) return(NA_character_)
  cat("Calculating consensus of", length(seq), "sequences...\n")
  tictoc::tic("total")
  on.exit(tictoc::toc())
  seq <- Biostrings::DNAStringSet(seq) %>%
    Biostrings::RNAStringSet()

  cat(" Aligning...\n")
  tictoc::tic("  alignment")
  aln <- DECIPHER::AlignSeqs(seq, processors = 3, verbose = FALSE)
  tictoc::toc()

  cat(" Removing outliers...\n")
  tictoc::tic("  outliers")
  outliers <- odseq::odseq(Biostrings::RNAMultipleAlignment(aln))
  cat("  -removed", sum(outliers), "/", length(outliers),
      "sequences as outliers.\n")
  aln <- aln[!outliers]
  tictoc::toc()

  cat(" Masking gaps...\n")
  tictoc::tic("  masking")
  aln <- aln %>%
    Biostrings::RNAMultipleAlignment() %>%
    Biostrings::maskGaps(min.fraction = 0.5, min.block.width = 1) %>%
    as("RNAStringSet")
  tictoc::toc()

  cat(" Calculating consensus...\n")
  tictoc::tic("  consensus")
  on.exit(tictoc::toc(), add = TRUE)
  DECIPHER::ConsensusSequence(aln,
                              threshold = 0.5,
                              ambiguity = TRUE,
                              ignoreNonBases = TRUE,
                              includeTerminalGaps = FALSE)
}

conseq <-
  combined %>%
  group_by(ITS2) %>%
  filter(!is.na(ITS2), n() >= 3) %>%
  summarize(nreads = n(),
            ITS = as.character(calculate_consensus2(ITS, seq.id)),
            LSU = as.character(calculate_consensus2(LSU, seq.id)))

conseq_filt <- conseq %>%
  dplyr::mutate(qITS = Biostrings::RNAStringSet(ITS) %>%
                  Biostrings::letterFrequency(., "MRWSYKVHDBN"),
                qLSU = Biostrings::RNAStringSet(LSU) %>%
                  Biostrings::letterFrequency(., "MRWSYKVHDBN")) %>%
  dplyr::filter(qITS < 3, qLSU < 3)


taxITS2 <- readd("taxon_ITS2_unite")

taxLSU <- conseq_filt$LSU %>%
  unique() %>%
  stringr::str_replace_all("U", "T") %>%
  dada2::assignTaxonomy(refFasta <- "reference/rdp.fasta.gz",
                        taxLevels = c("Kingdom", "Phylum", "Class", "Order",
                                      "Family", "Genus", "Species"),
                        multithread = FALSE) %>%
  tibble::as_tibble(rownames = "LSU")

cons_tax  <-
  dplyr::left_join(conseq_filt, select(taxITS2, -nreads),
                   by = c("ITS2" = "seq")) %>%
  dplyr::mutate(name.ITS = dplyr::coalesce(Genus, Family, Order, Class, Phylum,
                                           Kingdom, "Unknown"),
                name.ITS = str_replace_all(name.ITS, " ", "_")) %>%
  dplyr::group_by(name.ITS) %>%
  dplyr::mutate(newname = if (n() > 1) {
    paste(name.ITS, seq_along(name.ITS), sep = "_")
  } else {
    name.ITS
  }) %>%
  dplyr::ungroup() %>%
  dplyr::select(-name.ITS) %>%
  dplyr::rename(name.ITS = newname) %>%
  dplyr::left_join(taxLSU %>% mutate_at("LSU", stringr::str_replace_all, "T", "U"),
                   suffix = c("_ITS2", "_LSU"),
                   by = "LSU") %>%
  dplyr::mutate(name.LSU = dplyr::coalesce(Genus_LSU, Family_LSU, Order_LSU,
                                           Class_LSU, Phylum_LSU, Kingdom_LSU,
                                           "Unknown"),
                name.LSU = str_replace_all(name.LSU, " ", "_")) %>%
  dplyr::group_by(name.LSU) %>%
  dplyr::mutate(newname = if (n() > 1) {
    paste(name.LSU, seq_along(name.LSU), sep = "_")
  } else {
    name.LSU
  }) %>%
  dplyr::ungroup() %>%
  dplyr::select(-name.LSU) %>%
  dplyr::rename(name.LSU = newname) %>%
  tidyr::unite("name", name.ITS, name.LSU, sep = "/")

cons_tax %>%
  filter(nreads >= 3) %$%
  rlang::set_names(ITS, name) %>%
  RNAStringSet() %>%
  letterFrequency(., "MRWSYKVHDBN") %>%
  tibble::tibble(n = .) %>%
  ggplot(aes(n)) +
  stat_ecdf()

lsualn <-  cons_tax %$%
  rlang::set_names(LSU, name) %>%
  Biostrings::RNAStringSet() %>%
  DECIPHER::AlignSeqs(iterations = 10, refinements = 10)

lsudist <- DECIPHER::DistanceMatrix(lsualn, includeTerminalGaps = FALSE,
                                    correction = "Jukes-Cantor")
lsutree <- DECIPHER::IdClusters(lsudist, myXStringSet = lsualn,
                             method = "ML", showPlot = TRUE,
                             type = "dendrogram", processors = 3)

DECIPHER::WriteDendrogram(lsutree, "lsutree.tree")

cons_tax %$%
  rlang::set_names(paste0(ITS, LSU), name) %>%
  Biostrings::RNAStringSet() %>%
  Biostrings::writeXStringSet("longASV.fasta")

testset <- combined %>% group_by(ITS2) %>%
  group_split() %>%
  extract2(4)


ij <- its_join(bigmaps, verbose = TRUE)

calculate_consensus <- function(data) {
  assertthat::assert_that(is.data.frame(data),
                          assertthat::has_name(data, "seq"),
                          assertthat::has_name(data, "reads"))
  # hash the sequences to make unique names
  data <- data %>%
    dplyr::mutate(hash = map_chr(seq, digest::digest) %>%
                    stringr::str_sub(end = 8))
  # align using ClustalW
  # MUSCLE gives errors with over 100 sequences?
  aln <- set_names(data$seq, data$hash) %>%
    Biostrings::DNAStringSet() %>%
    msa::msaClustalW() %>%
    as.matrix()
  # Put the alignment in the data
  data <- data %>%
    dplyr::mutate(seq = purrr::map_chr(hash, ~seqinr::c2s(aln[.,])))
  # repeat each sequence equal to the number of times it was found, and
  # find a preliminary consensus; gaps in this consensus are removed.
  gaps <- rep(data$seq, data$reads) %>%
    Biostrings::DNAMultipleAlignment() %>%
    Biostrings::consensusString() %>%
    stringr::str_locate_all("-")
  gaps <- gaps[[1]]
  if (nrow(gaps) >= 1) {
    nongaps <- -gaps[,"start"]
  } else {
    nongaps <- TRUE
  }
  # Remove the gap columns and create a new consensus sequence, ignoring gaps.
  # Put the threshhold at 0,0: we always want the most common base.
  data %>%
    dplyr::mutate(seq = purrr::map_chr(hash, ~seqinr::c2s(aln[.,nongaps]))) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(
      seq = rep(seq, reads) %>%
        Biostrings::DNAMultipleAlignment() %>%
        msa::msaConsensusSequence(type = "upperlower",
                                  ignoreGaps = TRUE,
                                  thresh = c(0, 0)),
      reads = sum(reads))
}

group_cons <- ij %>%
  tidyr::gather(key = "region", value = "seq", ITS, ITS1, ITS2, LSU, long) %>%
  filter(!is.na(seq), chimera == FALSE) %>%
  group_by(group, region, seq) %>%
  summarize(reads = sum(reads)) %>%
  group_by(group, region) %>%
  group_map(
    function(data, key) {
      if (nrow(data) == 1) return(data)
      cat("Calculating consensus of", nrow(data), "different", key$region,
          "sequences for group", key$group, "\n")
      calculate_consensus(data)
    }
  )


taxLSU <- readd("taxon_LSU_rdp")
ijtax <- left_join(ij, taxITS2, by = c("ITS2" = "seq")) %>%
  group_by(group) %>%
  dplyr::summarize(name =
                     if (n_distinct(Species, na.rm = TRUE) == 1) {
                       unique(na.omit(Species))
                     } else if ( n_distinct(Genus, na.rm = TRUE) == 1) {
                       unique(na.omit(Genus))
                     } else if (n_distinct(Family, na.rm = TRUE) == 1) {
                       unique(na.omit(Family))
                     } else if (n_distinct(Order, na.rm = TRUE) == 1) {
                       unique(na.omit(Order))
                     } else if (n_distinct(Class, na.rm = TRUE) == 1) {
                       unique(na.omit(Class))
                     } else if (n_distinct(Phylum, na.rm = TRUE) == 1) {
                       unique(na.omit(Phylum))
                     } else {
                       "Unknown"
                     })

aln <- group_cons %>%
  dplyr::filter(region %in% c("ITS2", "LSU")) %>%
  dplyr::select(-reads) %>%
  tidyr::spread(key = region, value = seq) %>%
  dplyr::ungroup() %>%
  dplyr::filter(complete.cases(.)) %>%
  tidyr::unite("seq", ITS2, LSU, sep = "") %$%
  rlang::set_names(seq, group) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet() %>%
  msa::msaClustalW()

aln <- group_cons %>%
  dplyr::filter(region %in% c("LSU"), reads > 1) %>%
  dplyr::select(-reads) %>%
  tidyr::spread(key = region, value = seq) %>%
  dplyr::ungroup() %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::left_join(ijtax) %$%
  # tidyr::unite("seq", ITS2, LSU, sep = "") %$%
  rlang::set_names(LSU, name) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet() %>%
  DECIPHER::AlignSeqs()

dist <- DECIPHER::DistanceMatrix(aln, includeTerminalGaps = FALSE, correction = "Jukes-Cantor")
tree <- DECIPHER::IdClusters(dist, myXStringSet = aln,
                     method = "ML", showPlot = TRUE)
library(DECIPHER)
ij %>%
  tidyr::gather(key = "region", value = "seq", ITS, ITS1, ITS2, LSU, long) %>%
  filter(!is.na(seq), chimera == FALSE) %>%
  group_by(group, region, seq) %>%
  summarize(reads = sum(reads)) %>%
  group_by(group, region) %>%
  filter(group == 14, region == "LSU") %>%
  select(seq, reads) %>%
  calculate_consensus()


ij %>% group_by(group) %>% summarize_at(c("ITS", "ITS1", "ITS2", "LSU", "long"), n_distinct, na.rm = TRUE) %>% View()

tax <- readd(taxon_LSU_rdp)
tax2 <- readd(taxon_ITS2_unite)
sum(tax2$nreads)
sum(tax$nreads)
tax3 <- readd(taxon_long_rdp)
sum(tax3$nreads)
tax4 <- readd(taxon_ITS_unite)
tax1 <- readd(taxon_ITS1_unite)

table(tax$Genus) %>%
  sort(decreasing = TRUE) %>%
  magrittr::extract(. > 10) %>%
  names() %>%
  tibble::tibble(Genus = .,
                 data = map(., . %>% {filter(tax, Genus == .)} %$% seq %>%
                              set_names(map_chr(., digest::digest) %>%
                                          substr(start = 1, stop = 8)) %>%
                              DNAStringSet)) %>%
  pwalk(function(Genus, data) writeFasta(data, glue("output/{Genus}.fasta")))

#### distances ####

loadd(big_seq_table_ITS2)
loadd(platemap)
loadd(datasets)
loadd(conseq_filt)

# rename the rows and columns of a sequence abundance table
# convert rownames from a fileame to {Seq.Run}{Plate}{Well}
# convert colnames from the full sequence to an 8-character hash.
# also pools forward and reverse reads from same well.

relabel_seqtable <- function(seqtable, hashlen = 8) {
  seqtable %>%
    # set the new column names
    magrittr::set_colnames(seqhash(colnames(.), len = hashlen)) %>%
    # convert to a tibble for easier name column operations
    tibble::as_tibble(rownames = "file") %>%
    # parse the filename
    tidyr::extract(
      col = "file",
      into = c("Seq.Run", "Plate", "Well", "Direction", "Region"),
      regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.qfilt\\.fastq\\.gz") %>%
    # create the ID
    tidyr::unite("ID", Seq.Run, Plate, Well, sep = "") %>%
    # remove unnecessary columns
    dplyr::select(-Direction, -Region) %>%
    # pool all reads from the same well
    dplyr::group_by(ID) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::ungroup() %>%
    # convert back to a matrix
    tibble::column_to_rownames("ID") %>%
    as.matrix
}

assemble_physeq <- function(platemap, datasets, seqtable, tree) {
  samp <- platemap %>%
    dplyr::mutate_at("Primer.Pair", tolower) %>%
    dplyr::left_join(
      datasets %>%
        dplyr::select(Dataset, Seq.Run, Tech, Forward, Reverse) %>%
        dplyr::mutate(
          Primer.Pair = paste(stringr::str_replace(Forward, "_tag.*$", ""),
                              stringr::str_replace(Reverse, "_tag.*$", ""),
                              sep = "_")),
      by = "Primer.Pair") %>%
    dplyr::mutate_at(c("Year", "Site", "Qual", "Plate", "Well", "Primer.Pair"),
                     factor) %>%
    tidyr::unite("ID", Seq.Run, Plate, Well, sep = "", remove = FALSE) %>%
    tibble::column_to_rownames("ID") %>%
    phyloseq::sample_data()

  asvtab <- seqtable %>%
    phyloseq::otu_table(taxa_are_rows = FALSE)

  phyloseq::phyloseq(samp, asvtab, tree)
}

%>%
  phyloseq::prune_samples(samples = phyloseq::sample_data(.)$Dataset == "long-pacbio"
                          & phyloseq::sample_data(.)$Qual == ""
                          & phyloseq::sample_data(.)$sample_type == "Sample")

st <- relabel_seqtable(big_seq_table_ITS2)

loadd(longtree)

tax_table(physeq) <- taxon_ITS2_unite %>%
  mutate_at("seq", seqhash) %>%
  select(-Taxonomy, -nreads) %>%
  column_to_rownames("seq") %>%
  as.matrix() %>%
  phyloseq::tax_table()

treetax <- taxon_ITS2_unite %>%
  mutate_at("seq", seqhash) %>%
  select(-Taxonomy, -nreads)

library(ggtree)
p <- ggtree(tr = longtree) %<+% treetax



outgroups <- treetax %>% filter(!is.na(Kingdom), Kingdom %in% c("Alveolata", "Chromista", "Plantae"), seq %in% longtree$tip.label)
MRCA(p, outgroups$seq)
reroot(longtree, MRCA(p, outgroups$seq))
ggtree(groupOTU(longtree, split(treetax$seq, treetax$Kingdom))) + aes(color = group) + scale_color_brewer(type = "qual") +
theme(legend.position = NULL)
p + aes(color = Kingdom)

reconstruct_ITS <-
  readd(dada_map_pb_500_001_ITS1) %>%
  select(seq.id, ITS1 = dada.seq) %>%
  inner_join(
    readd(dada_map_pb_500_001_5_8S) %>%
               select(seq.id, `5_8S` = dada.seq),
             by = "seq.id") %>%
  inner_join(
    readd(dada_map_pb_500_001_ITS2) %>%
               select(seq.id, ITS2 = dada.seq),
             by = "seq.id") %>%
  full_join(
    readd(dada_map_pb_500_001_ITS) %>%
      select(seq.id, ITS = dada.seq),
    by = "seq.id") %>%
  mutate(ITS = coalesce(ITS, paste0(ITS1, `5_8S`, ITS2))) %>%
  select(seq.id, ITS)

is_bimera_ITS <- reconstruct_ITS %>%
  group_by(ITS) %>%
  summarize(abundance = n()) %>%
  ungroup() %>%
  dplyr::rename(sequence = ITS) %>%
  dada2::isBimeraDenovo()

reconstruct_LSU <-
  readd(dada_map_pb_500_001_LSU1) %>%
  select(seq.id, LSU1 = dada.seq) %>%
  inner_join(
    readd(dada_map_pb_500_001_D1) %>%
      select(seq.id, D1 = dada.seq),
    by = "seq.id") %>%
  inner_join(
    readd(dada_map_pb_500_001_LSU2) %>%
      select(seq.id, LSU2 = dada.seq),
    by = "seq.id") %>%
  inner_join(
    readd(dada_map_pb_500_001_D2) %>%
      select(seq.id, D2 = dada.seq),
    by = "seq.id") %>%
  inner_join(
    readd(dada_map_pb_500_001_LSU3) %>%
      select(seq.id, LSU3 = dada.seq),
    by = "seq.id") %>%
  inner_join(
    readd(dada_map_pb_500_001_D3) %>%
      select(seq.id, D3 = dada.seq),
    by = "seq.id") %>%
  inner_join(
    readd(dada_map_pb_500_001_LSU4) %>%
      select(seq.id, LSU4 = dada.seq),
    by = "seq.id") %>%
  full_join(
    readd(dada_map_pb_500_001_LSU) %>%
      select(seq.id, LSU = dada.seq),
    by = "seq.id") %>%
  mutate(LSU = coalesce(LSU, paste0(LSU1, D1, LSU2, D2, LSU3, D3, LSU4))) %>%
  select(seq.id, LSU)

is_bimera_LSU <- reconstruct_LSU %>%
  group_by(LSU) %>%
  summarize(abundance = n()) %>%
  ungroup() %>%
  dplyr::rename(sequence = LSU) %>%
  dada2::isBimeraDenovo()

reconstruct_long <-
  inner_join(
    reconstruct_ITS %>%
      filter(!ITS %in% names(is_bimera_ITS)[is_bimera_ITS]),
    reconstruct_LSU %>%
      filter(!LSU %in% names(is_bimera_LSU)[is_bimera_LSU]),
    by = "seq.id"
  ) %>%
  full_join(
    readd(dada_map_pb_500_001_long) %>%
      select(seq.id, long = dada.seq),
    by = "seq.id"
  ) %>%
  mutate(long = coalesce(long, paste0(ITS, LSU))) %>%
  select(seq.id, long)

is_bimera_long <- reconstruct_long %>%
  group_by(long) %>%
  summarize(abundance = n()) %>%
  ungroup() %>%
  dplyr::rename(sequence = long) %>%
  dada2::isBimeraDenovo()

dada2::getUniques(reconstruct_ITS$ITS)
reconstruct_ITS %>%
  group_by(ITS) %>%
  summarize(abundance = n()) %>%
  ungroup() %>%
  dplyr::rename(sequence = ITS) %>%
  dada2::isBimeraDenovo()

reconstruct <- function(
  seqtabs,
  regions,
  output = "concat",
  order = setdiff(regions, output),
  read_column = "seq.id",
  asv_column = "dada.seq",
  raw_column = NULL,
  sample_column = NULL,
  sample_regex = NULL,
  sample_replace = NULL,
  chimera_offset = 0,
  ...)
{
  assertthat::assert_that(assertthat::is.string(read_column),
                          assertthat::is.string(asv_column),
                          is.null(raw_column) || is.na(raw_column) ||
                            assertthat::is.string(raw_column),
                          is.null(sample_column) || is.na(sample_column) ||
                            assertthat::is.string(sample_column),
                          is.null(sample_regex) || is.na(sample_regex) ||
                            assertthat::is.string(sample_regex),
                          is.null(sample_replace) || is.na(sample_replace) ||
                            assertthat::is.string(sample_replace))

  if (is.null(raw_column) || is.na(raw_column)) raw_column <- NULL
  if (is.null(sample_column) || is.na(sample_column)) sample_column <- NULL
  if (is.null(sample_regex) || is.na(sample_regex)) sample_regex <- NULL
  if (is.null(sample_replace) || is.na(sample_replace)) sample_replace <- NULL
  seqtabs <- purrr::map2(
    seqtabs,
    regions,
    function(st, reg) magrittr::set_names(
      st[,c(sample_column, read_column, asv_column, raw_column)],
      c(sample_column, read_column, reg, if (is.null(raw_column)) NULL else paste0(reg, "_raw"))))
  seqtabs <- purrr::map(
    unique(regions) %>% magrittr::set_names(., .),
    function(r) {
      dplyr::bind_rows(seqtabs[regions == r])
    }
  )
  if (!is.null(sample_regex)) {
    if (is.null(sample_column)) stop("If sample_regex is give, sample_column must also be given.")
    if (is.null(sample_replace)) {
      seqtabs <- purrr::map(
        seqtabs,
        dplyr::mutate_at,
        sample_column,
        stringr::str_extract,
        sample_regex
      )
    } else {
      seqtabs <- purrr::map(
        seqtabs,
        dplyr::mutate_at,
        sample_column,
        stringr::str_replace,
        sample_regex,
        sample_replace,
      )

    }
  }
  out <- purrr::reduce(
    seqtabs[order],
    dplyr::full_join,
    by = c(sample_column, read_column)
  )
  if (output %in% regions) {
    out <-
      dplyr::full_join(out, seqtabs[[output]], by = c(sample_column, read_column))
  }
  combos <- dplyr::group_by_at(out, regions) %>%
    dplyr::summarize(nreads = dplyr::n()) %>%
    dplyr::ungroup()

  for (i in seq(1, max(1, length(order) - 2), 2)) {
    chimset <- order[i:min(length(order), i + 2)]
    seqs <- do.call(stringr::str_c, out[,chimset])
    chims <-
      if (!is.null(sample_column)) {
        tibble::tibble(sample = out[[sample_column]], seq = seqs) %>%
          dplyr::filter(!is.na(seq)) %>% dplyr::group_by(sample, seq) %>%
          dplyr::summarize(nread = dplyr::n()) %>%
          dplyr::ungroup() %>%
          tidyr::spread(seq, nread, fill = 0L) %>%
          tibble::column_to_rownames("sample") %>%
          as.matrix() %>%
          dada2::isBimeraDenovoTable(...)
      } else {
        table(seqs) %>%
        {tibble::tibble(abundance = ., sequence = names(.))} %>%
          dada2::isBimeraDenovo(...)
      }
    chims <- names(chims)[chims]
    futile.logger::flog.info(
      "Removing %d/%d reads of %d/%d chimeric sequences from domains: %s.",
      sum(seqs %in% chims),
      length(seqs),
      length(chims),
      dplyr::n_distinct(seqs, na.rm = TRUE),
      paste(chimset, collapse = ", ")
    )
    out <- out[!seqs %in% chims,]
  }
  if (output %in% regions) {
    out[["_concat_"]] <- do.call(stringr::str_c, out[,order])
    out[[output]] <- dplyr::coalesce(out[[output]], out[["_concat_"]])
    out[["_concat_"]] <- NULL
  } else {
    out[[output]] <- do.call(stringr::str_c, out[,order])
  }
  out
}


debugonce(reconstruct)
reconstruct_long <- reconstruct(
  list(
    readd(dada_map_pb_500_001_ITS1),
    readd(dada_map_pb_500_001_5_8S),
    readd(dada_map_pb_500_001_ITS2),
    readd(dada_map_pb_500_001_LSU1),
    readd(dada_map_pb_500_001_D1),
    readd(dada_map_pb_500_001_LSU2),
    readd(dada_map_pb_500_001_D2),
    readd(dada_map_pb_500_001_LSU3),
    readd(dada_map_pb_500_001_D3),
    readd(dada_map_pb_500_001_LSU4),
    readd(dada_map_pb_500_001_long)),
  regions = c("ITS1", "5_8S", "ITS2", "LSU1", "D1", "LSU2", "D2", "LSU3", "D3",
              "LSU4", "long"),
  output = "long")

reconstruct_ITS <- reconstruct(
  list(
    readd(dada_map_pb_500_001_ITS1),
    readd(dada_map_pb_500_001_5_8S),
    readd(dada_map_pb_500_001_ITS2),
    readd(dada_map_pb_500_001_ITS)),
  regions = c("ITS1", "5_8S", "ITS2", "ITS"),
  output = "ITS")

reconstruct_ITS <- reconstruct(
  list(
    readd(dada_map_pb_500_001_ITS1),
    readd(dada_map_pb_500_001_5_8S),
    readd(dada_map_pb_500_001_ITS2),
    readd(dada_map_pb_500_001_ITS),
    readd(dada_map_pb_500_002_ITS1),
    readd(dada_map_pb_500_002_5_8S),
    readd(dada_map_pb_500_002_ITS2),
    readd(dada_map_pb_500_002_ITS)),
  regions = c("ITS1", "5_8S", "ITS2", "ITS", "ITS1", "5_8S", "ITS2", "ITS"),
  output = "ITS",
  sample_column = "name",
  sample_regex = "pb_500_00[12]-[A-H]1?[0-9]"
)

targets <- drake_cache()$list()
dada_maps <- targets[startsWith(targets, "dada_map_pb_500")]
dada_map_regions <- stringr::str_replace(dada_maps, "dada_map_pb_500_00[12]_", "")

reconstruct_long <-
  lapply(dada_maps, readd, character_only = TRUE) %>%
  reconstruct(
    regions = dada_map_regions,
    output = "long",
    order = c("ITS1", "5_8S", "ITS2", "LSU1", "D1", "LSU2", "D2", "LSU3", "D3",
              "LSU4"),
    sample_column = "name",
    sample_regex = "pb_500_00[12]-[A-H]1?[0-9]"
  )

reconstruct_32S <-
  dplyr::mutate(
    reconstruct_long,
    `32S` = stringr::str_c(`5_8S`, ITS2, LSU1, D1, LSU2, D2, LSU3, D3, LSU4),
    long = dplyr::coalesce(stringr::str_c(ITS1, `32S`), long)
  ) %>%
  dplyr::select(long, `32S`) %>%
  dplyr::filter(complete.cases(.)) %>%
  unique()

reconstruct_32s_aln <-
  reconstruct_32S %$%
  magrittr::set_names(`32S`, tzara::seqhash(long)) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet() %T>%
  Biostrings::writeXStringSet("temp/32S_recon.fasta") %>%
  cmalign(cmfile = 'reference/fungi_32S_LR5.cm', seq = ., glocal = TRUE, cpu = 4)

reconstruct_long_aln <-
  reconstruct_long %>%
  dplyr::select(long, ITS1) %>%
  dplyr::filter(complete.cases(.)) %>%
  unique() %>%
  dplyr::mutate(hash = tzara::seqhash(long)) %>%
  dplyr::inner_join(
    tibble::tibble(
      hash = reconstruct_32s_aln$names,
      aln = as.character(reconstruct_32s_aln$alignment@unmasked)
    ),
    by = "hash"
  ) %>%
  dplyr::mutate(
    ITS1 = chartr("T", "U", ITS1),
    ITS1 = stringr::str_pad(ITS1, max(nchar(ITS1)), "right", "-"),
    aln = stringr::str_c(ITS1, aln)
  )

long_recon_aln <- reconstruct_long_aln %$%
  magrittr::set_names(aln, hash) %>%
  Biostrings::RNAStringSet()

write_clustalw_ss(
  aln = long_recon_aln,
  sec_str = paste0(
    stringr::str_pad("xxxx", max(nchar(reconstruct_long_aln$ITS1)), "right", "-"),
    reconstruct_32s_aln$SS_cons
  ),
  ref = stringr::str_pad(chartr("v", ".", reconstruct_32s_aln$RF),
                         max(nchar(reconstruct_long_aln$aln)), "left", "-"),
  seq_names = reconstruct_long_aln$hash,
  file = "temp/long_recon.aln")


long_recon <-
  unique(reconstruct_long$long) %>%
  set_names(tzara::seqhash(.)) %>%
  purrr::discard(is.na) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::RNAStringSet()

long_recon_aln <- DECIPHER::AlignSeqs(long_recon, iterations = 10, refinements = 10, processors = 4)

Biostrings::writeXStringSet(long_recon_aln, "temp/long_recon_aln.fasta")

for (region in c("ITS1", "ITS2", "ITS", "LSU"))
  for (database in c("unite", "rdp_train", "warcup"))
    for (method in c("dada2", "sintax", "idtaxa")) {
      tname <- paste("taxon", region, database, region, method, sep = "_")
      if (tname %in% targets) {
        taxa <- readd(tname, character_only = TRUE)
        pcount <- dplyr::n_distinct(
          taxa$label[!is.na(taxa$taxon) &
                        taxa$rank == "phylum" &
                        taxa$confidence >= 0.5],
          na.rm = TRUE
        )
        lall <- dplyr::n_distinct(taxa$label, na.rm = TRUE)
        gcount <- dplyr::n_distinct(
          taxa$label[!is.na(taxa$taxon) &
                        taxa$rank == "genus" &
                        taxa$confidence >= 0.5],
          na.rm = TRUE
        )
        cat(region, database, method, "phylum: ", pcount, "/", lall, "genus: ", gcount, "/", lall,"\n")
      }
    }

seq <- Biostrings::DNAStringSet(
  Biostrings::readQualityScaledDNAStringSet(
    "sequences/trim/pb_500_001/pb_500_001-H6f.trim.fastq.gz"
  )
)
LSUx(seq, cm_5.8S = "reference/RF00002.cm", cm_32S = "reference/fungi_32S_LR5.cm", cpu = 3)
loadd(preconseq)

# combine_diskframe test ----
library(drake)
library(magrittr)
clean(dataset, observations, all_data, types, means)

combine_diskframe <- function(..., dir = drake_tempfile()) {
  dir.create(dir)
  shards <- rlang::ensyms(...)
  cache <- drake_cache()
  for (i in seq_along(shards)) {
    shardfile <- cache$file_return_key(rlang::expr_text(shards[[i]]))
    file.link(shardfile, file.path(dir, paste0(i, ".fst")))
  }
  disk.frame::disk.frame(dir)
}

create_tabular_data <- function(x, n) {
  data.frame(
    type = sample(letters[1:3], n, replace = TRUE),
    size = runif(n),
    stringsAsFactors = FALSE
  )
}

dataset = 1:3
plan <- drake_plan(
  observations = target(
    create_tabular_data(dataset, 5),
    transform = map(dataset = !!dataset),
    format = "fst"
  ),

  all_data = target(
    combine_diskframe(observations),
    transform = combine(observations),
    format = "diskframe"
  ),

  result = target(
    all_data %>%
      as.data.frame() %>%
      dplyr::group_by(type) %>%
      dplyr::summarize(mean = mean(size))
  )
)

make(plan)
readd(observations_1L)
readd(observations_2L)
readd(observations_3L)
readd(all_data)
readd(result)

# combine_dynamic_diskframe test ----
library(drake)
library(magrittr)
library(disk.frame)

combine_dynamic_diskframe <- function(dd, dir = drake_tempfile()) {
  dir.create(dir)
  cache <- drake_cache()
  for (i in seq_along(dd)) {
    shardfile <- file.path(cache$path_return, dd[i])
    stopifnot(file.exists(shardfile))
    file.link(shardfile, file.path(dir, paste0(i, ".fst")))
  }
  disk.frame::disk.frame(dir)
}

create_tabular_data <- function(x, n) {
  data.frame(
    type = sample(letters[1:3], n, replace = TRUE),
    size = runif(n),
    stringsAsFactors = FALSE
  )
}

plan <- drake_plan(
  dataset = 1:3,
  observations = target(
    create_tabular_data(dataset, 5),
    dynamic = map(dataset),
    format = "fst"
  ),

  all_data = target(
    combine_dynamic_diskframe(observations),
    format = "diskframe"
  ),

  result = target(
    all_data %>%
      as.data.frame() %>%
      dplyr::group_by(type) %>%
      dplyr::summarize(mean = mean(size))
  )
)

make(plan)
readd(all_data)
readd(result)

# just a diskframe test ----
library(drake)
library(magrittr)
n <- 200
observations = data.frame(
  type = sample(letters[1:3], n, replace = TRUE),
  size = runif(n),
  stringsAsFactors = FALSE
)

plan <- drake_plan(
  all_data = target(
    observations,
    format = "diskframe"
  ),
  result = target(
    all_data %>%
      disk.frame::chunk_group_by(type) %>%
      disk.frame::chunk_summarize(mean = mean(size)) %>%
      as.data.frame()
  )
)

make(plan)
readd(result)

long_aln <- Biostrings::readDNAMultipleAlignment("data/raxml/decipher/decipher_long.phy", format = "phylip") %>%
  Biostrings::DNAStringSet()

short_reads <- Biostrings::readDNAStringSet("data/raxml/decipher/short.fasta")

full_aln4 <- mafft_add(
  long_aln,
  short_reads[1:1000],
  add = "add",
  method = "10merpair",
  maxiterate = "1000",
  thread = 3,
  quiet = FALSE,
  exec = Sys.which("mafft")
)

full_aln2 <- ips::mafft(
  long_aln,
  short_reads[1:100],
  add = "addfragments",
  gt = ape::read.tree("data/raxml/decipher/RAxML_bestTree.decipher_long"),
  method = "localpair",
  maxiterate = 1000,
  thread = 3,
  quiet = FALSE,
  exec = Sys.which("mafft")
)


tree <- ape::read.tree("data/raxml/decipher/RAxML_bestTree.decipher_long")
aln <- Biostrings::readDNAStringSet("data/raxml/decipher/long_aln.fasta") %>%
  magrittr::extract(names(.) %in% tree$tip.label)
tree <- ape::keep.tip(tree, names(aln)) %T>%
  ape::write.tree("data/raxml/decipher/test_long.tree")
Biostrings::writeXStringSet(aln, "data/raxml/decipher/long_aln.fasta")
short_aln <- Biostrings::readDNAStringSet("data/raxml/decipher/short_aln.fasta")

newtree <- gappa_graft(epa_ng_result)

polytree <- grafts_to_polytomies(newtree, tree)

rax3 <- ips::raxml(as.matrix(ape::as.DNAbin(c(aln, short_aln))),
                   m = "GTRGAMMA", f = "a", N = "autoMRE_IGN",
                   p = 12345, x = 12345, backbone = polytree, file = "full",
                   exec = "raxmlHPC-PTHREADS-SSE3")

readd(allseqs) %>%
  dplyr::select(hash, full) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::arrange(dplyr::desc(nchar(full))) %>%
  dplyr::filter(!duplicated(full)) %>%
  dplyr::filter(!duplicated(hash)) %>%
  dplyr::filter(!hash %in% names(readd(aln_decipher_long_trim))) %$%
  set_names(full, hash) %>%
  chartr("U", "T", .) %>%
  Biostrings::DNAStringSet()

mafft_add(
  x = readd(aln_decipher_long_trim),
  y = shortseqs[1:100],
  add = "add",
  method = "10merpair",
  maxiterate = "1000",
  thread = 1,
  quiet = FALSE,
  exec = Sys.which("mafft")
)


library(drake)

cache1 <- new_cache("cache1")
cache2 <- new_cache("cache2")

plan1 <- drake_plan(
  x = target(
    rcauchy(1e6),
    transform = map(l = !!letters[1:3], .id = l)
  ),
  y = target(
    mean(x),
    transform = map(x, .id = l)
  )
)
make(plan1, cache = cache1)

cache2$import(from = cache1, list = plan1$target[startsWith(plan1$target, "y")])
# now copy cache2 to another machine, which does not have access to cache1

plan2 = drake_plan(
  y = target(trigger = trigger(mode = "blacklist"), transform = map(l = !!letters[1:3], .id = l)),
  z = target(c(y), transform = combine(y))
)

make(plan2, cache = cache2)
readd(z, cache = cache2)

library(drake)
options(clustermq.scheduler = "multisession")

t = c(1, 1, 1, 10)

plan = drake_plan(
  x = target({
    Sys.sleep(t)
  },
  transform = map(t = !!t, .id = FALSE)
  )
)
log <- tempfile()
make(plan, parallelism = "clustermq", jobs = 4, console_log_file = log, template = list(log_file = log))
readLines(log)

c(
  "taxon_table",
  dplyr::filter(plan, step == "taxon")$target,
  "raxml_decipher_LSU",
  "raxml_decipher_long",
  "raxml_epa_full",
  dplyr::filter(plan, step == "big_seq_table")$target,
  dplyr::filter(plan, step == "err")$target,
  dplyr::filter(plan, step == "chimeras")$target
)

cache <- drake_cache(".light")

epa_backbone <- ape::read.tree("data/raxml/epa_long-backbone.tre")
epa_short <-
  Biostrings::readDNAStringSet("data/raxml/epa_short_unaln.fasta") %>%
  Biostrings::RNAStringSet()

epa_short_aln <- DECIPHER::AlignSeqs(
  epa_short,
  iterations = 10,
  refinements = 10,
  processors = 4
)
Biostrings::writeXStringSet(epa_short_aln, "data/raxml/epa_short_realn.fasta")
epa_short_aln <- Biostrings::readRNAStringSet("data/raxml/epa_short_realn.fasta")
long_aln <- Biostrings::readDNAMultipleAlignment(
  "data/raxml/decipher/decipher_long.phy",
  format = "phylip"
)@unmasked %>%
  Biostrings::RNAStringSet()

full_aln <- DECIPHER::AlignProfiles(epa_short_aln, long_aln) %>%
  DECIPHER::AdjustAlignment()

mafft_full_aln <- Biostrings::readDNAMultipleAlignment(
  "data/raxml/epa_long.phy", format = "phylip")

loadd(raxml_decipher_long, cache = cache)

epa_decipher_full <- epa_ng(full_aln[names(long_aln)], raxml_decipher_long$bestTree, full_aln[setdiff(names(full_aln), names(long_aln))], model = raxml_decipher_long$info, threads = 4)

epa_decipher_full_graft <- gappa_graft(epa_decipher_full, threads = 4)

loadd(taxon_labels, cache = cache)
relabel_tree(epa_decipher_full_graft, taxon_labels$label, paste0('"', taxon_labels$tip_label, '"'), chimeras = readd(allchimeras_ITS2, cache = cache)) %>%
  castor::write_tree("data/trees/labeled_epa_decipher.tree")


loadd(physeq_all, cache = cache)
coph <- max_cophenetic(tree)[-seq_along(tree$tip.label),] %>%
  dplyr::arrange(max_coph) %>%
  dplyr::group_by(max_coph) %>%
  dplyr::summarize(n = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ASVs = length(tree$tip.label) - cumsum(n))
ggplot(coph, aes(max_coph, ASVs)) +
  geom_vline(xintercept = 0.0042, linetype = 2, color = "gray") +
  geom_step() +
  xlim(c(0, 0.05)) +
  xlab("cophenetic distance threshold")

readd(taxon_table, cache = cache)
readd(taxon_phy_long, cache = cache)$tip_taxa %>%
  dplyr::group_by_at(dplyr::vars(-n_reads)) %>%
  dplyr::summarize(n_reads = sum(n_reads)) %>%
  dplyr::ungroup()

node_taxa <- readd(taxon_phy_long, cache = cache)$node_taxa
tip_taxa <- readd(taxon_phy_long, cache = cache)$tip_taxa %>%
  dplyr::group_by_at(dplyr::vars(-n_reads)) %>%
  dplyr::summarize_at("n_reads", sum) %>%
  dplyr::group_by(label, rank) %>%
  dplyr::mutate(
    n_diff = dplyr::n_distinct(taxon),
    n_tot = dplyr::n(),
    n_method = dplyr::n_distinct(method, na.rm = TRUE),
    n_reference = dplyr::n_distinct(reference, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()
tree <- ape::keep.tip(
  readd(fungi_tree_epa_full, cache = cache),
  intersect(
    readd(fungi_tree_epa_full, cache = cache)$tip.label,
    readd(raxml_decipher_long, cache = cache)$bestTree$tip.label
  )
)

taxon_table <- readd(taxon_table, cache = cache) %>%
  dplyr::group_by_at(dplyr::vars(-n_reads)) %>%
  dplyr::summarize_at("n_reads", sum) %>%
  dplyr::ungroup() %>%
  dplyr::select(-(region:method)) %>%
  tidyr::extract(
    col = name,
    into = c("region", "reference", "ref_region", "method"),
    regex = "taxon_(32S|5_8S|ITS[12]?|LSU|long|short)_(unite|warcup|rdp_train)_(ITS[12]?|LSU)_(idtaxa|dada2|sintax)",
    remove = FALSE
  ) %>%
  dplyr::group_by(label, rank) %>%
  dplyr::filter((!"ITS" %in% region) | region != "short") %>%
  dplyr::mutate(
    n_tot = dplyr::n(),
    n_diff = dplyr::n_distinct(taxon, na.rm = TRUE),
    n_method = dplyr::n_distinct(method, na.rm = TRUE),
    n_reference = dplyr::n_distinct(reference, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

confident_IDs <- taxon_table %>%
  dplyr::filter(
    n_method >= 2,
    n_reference >= 2,
    n_diff == 1
  ) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(
    is_fungi = "Fungi" %in% taxon | !"kingdom" %in% rank
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    ifelse(
      is_fungi,
      rank == "genus",
      rank == "kingdom"
    )
  ) %$%
  label %>%
  unique()

taxon_table %>%
  dplyr::filter(label %in% confident_IDs) %>%
  dplyr::select(label, rank, taxon) %>%
  unique() %>%
  dplyr::arrange(rank) %>%
  dplyr::group_by(label) %>%
  dplyr::summarize(Taxonomy = paste(taxon, collapse = ";")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Taxonomy = ifelse(
      stringr::str_detect(Taxonomy, "^[A-Z][a-z]+mycota"),
      paste0("Fungi;", Taxonomy),
      Taxonomy
    )
  )

# When a group is poly-/paraphyletic in the tree, there are several possibilities:
# 1. The tree topology is wrong (solution: constrain tree to match taxonomy?)
# 2. Something is misannotated in the database (unlikely if the ID comes from multiple databases)
# 3. Taxonomy assignment is incorrect (unlikely if multiple methods agree)
# 4. Taxon is really para/polyphyletic (most likely for smaller/less studied taxa)

paraphy_groups <- node_taxa %>%
  dplyr::group_by(taxon, rank) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::group_modify(~tidyr::crossing(node1 = .x$node, node2 = .x$node) %>%
                        dplyr::filter(node2 > node1)) %>%
  dplyr::mutate(
    mrca = phangorn::mrca.phylo(
      x = tree,
      node = c(node1, node2),
      full = TRUE
    ),
    n1 = phangorn::Descendants(tree, node1) %>%
      purrr::map_int(length),
    n2 = phangorn::Descendants(tree, node2) %>%
      purrr::map_int(length),
    tips = phangorn::Descendants(tree, mrca),
    n3 = tips %>%
      purrr::map_int(length)

  ) %>%
  tidyr::unnest(tips) %>%
  dplyr::mutate(
    label = tree$tip.label[tips]
  ) %>%
  dplyr::left_join(tip_taxa, by = c("rank", "label"), suffix = c("", "_tip")) %>%
  dplyr::filter(taxon != taxon_tip)
  dplyr::mutate(
    tips = tip_taxa %>%
      dplyr::filter(
        label %in% tree$tip.label[phangorn::Descendants(
          tree,
          node = mrca
        )[[1]]]
      ) %>%
      dplyr::rename(tip_rank = rank, tip_taxon = taxon) %>%
      list()
  ) %>%
  tidyr::unnest(tips) %>%
  dplyr::filter(tip_rank == rank, tip_taxon != taxon)

readd(taxon_phy_long, cache = cache)$tip_taxa

guilds_long <-
  tip_taxa %>%
  dplyr::select(label, rank, n_diff, taxon) %>%
  dplyr::mutate_at("rank", rank_factor) %>%
  dplyr::arrange(rank) %>%
  dplyr::group_by(label) %>%
  dplyr::filter(dplyr::cumall(n_diff == 1)) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(., key = rank, value = taxon) %>%
  dplyr::mutate(
    kingdom = dplyr::if_else(
      !is.na(phylum) & endsWith(phylum, "mycota"),
      "Fungi",
      kingdom
    )
  ) %>%
  dplyr::left_join(
    .,
    tidyr::gather(., key = "rank", value = "taxon", kingdom:genus) %>%
      dplyr::group_by(label) %>%
      dplyr::summarize(Taxonomy = paste(taxon, collapse = ";")),
    by = "label"
  ) %>%
  FUNGuildR::funguild_assign(readd(funguild_db, cache = cache))

ecm_tips <- guilds_long %>%
  dplyr::filter(stringr::str_detect(guild, "Ectomycorrhizal"))

guilds_all <-
  taxon_table %>%
  dplyr::select(label, rank, taxon, n_diff) %>%
  unique() %>%
  dplyr::arrange(rank) %>%
  dplyr::group_by(label) %>%
  dplyr::filter(dplyr::cumall(n_diff == 1)) %>%
  dplyr::ungroup() %>%
  tidyr::spread(., key = rank, value = taxon) %>%
  dplyr::mutate(
    kingdom = dplyr::if_else(
      !is.na(phylum) & endsWith(phylum, "mycota"),
      "Fungi",
      kingdom
    )
  ) %>%
  dplyr::left_join(
    .,
    tidyr::gather(., key = "rank", value = "taxon", kingdom:genus) %>%
      dplyr::group_by(label) %>%
      dplyr::summarize(Taxonomy = paste(taxon, collapse = ";")),
    by = "label"
  ) %>%
  FUNGuildR::funguild_assign(readd(funguild_db, cache = cache))

ecm_tips_all <- guilds_all %>%
  dplyr::filter(stringr::str_detect(guild, "Ectomycorrhizal"))

longtree <- ape::read.tree("data/trees/labeled_fungi_decipher_long.tree")
longtree$tip.label <- substr(longtree$tip.label, start = 1, stop = 8)

####
library(magrittr)
full_aln <- Biostrings::readDNAMultipleAlignment("data/raxml/epa_long.phy", format = "phylip")@unmasked %>% Biostrings::RNAStringSet()
subject_tree <- ape::read.tree("data/raxml/decipher/RAxML_bestTree.decipher_long")
subject <- full_aln[!startsWith(as.character(full_aln), "-") & names(full_aln) %in% subject_tree$tip.label]
subject_tree <- ape::keep.tip(subject_tree, names(subject))
query <- full_aln[!names(full_aln) %in% subject_tree$tip.label]
subject_model <- readLines("data/raxml/decipher/RAxML_info.decipher_long")
threads = 3

####
n <- 10000
labels <- paste0("t", seq.int(n))

# generate a random alignment
aln <- sample(Biostrings::DNA_BASES, 100 * n, replace = TRUE)
aln <- matrix(aln, nrow = n)
aln <- apply(aln, 1, paste, collapse = "")
names(aln) <- labels
aln <- Biostrings::DNAStringSet(aln)

# create a very deep tree (each node has one singleton child)
x <- rnorm(n, sd = 0.001)
x <- rev(x)
names(x) <- labels
dist <- dist(x)
tree <- hclust(dist)
tree <- as.dendrogram(tree)

DECIPHER::AlignSeqs(
  myXStringSet = aln,
  guideTree = tree,
  iterations = 0,
  refinements = 0
)

library(magrittr)
library(drake)
dplyr::inner_join(
  dplyr::bind_rows(
    readd(soil_samples),
    readd(control_samples)
  ) %>% select(sample_alias, "project name", plate, well, amplicon, sample_type),
  dplyr::bind_rows(
  readd(file_meta_gits7its4) %>%
    mutate(plate = ifelse(tech == "Ion Torrent", list(c("001", "002")), as.list(plate))) %>%
    unnest(plate),
  readd(file_meta_its1lr5)
  ),
  by = c("plate", "well", "amplicon")
) %>%
  dplyr::transmute(
    STUDY = `project name`,
    SAMPLE = sample_alias,
    NAME = paste(sample_alias, tech, direction, sep = "-") %>% str_replace("-$", "") %>% str_replace(" ", ""),
    INSTRUMENT = paste(tech, machine) %>% str_replace("Ion S5", "S5"),
    LIBRARY_SOURCE = ifelse(sample_type == "Pos", "GENOMIC", "METAGENOMIC"),
    LIBRARY_SELECTION = "PCR",
    LIBRARY_STRATEGY = "AMPLICON",
    FASTQ = trim_file,
    DESCRIPTION = case_when(
      direction == "f" ~ "Reads originally in forward orientation",
      direction == "r" ~ "Reads originally in reverse orientation",
      TRUE ~ "")
  ) %>%
  chop(cols = c(SAMPLE, NAME, LIBRARY_SOURCE)) %>%
  mutate_at("LIBRARY_SOURCE", map, unique) %>%
  rowwise() %>%
  group_split() %>%
  walk(
    function(d) {
      d <- map(as.list(d), unlist)
      d <- map(as.list(d), unique)
      if (length(d$SAMPLE) > 1) {
        d$DESCRIPTION <- paste("also contains samples from", d$SAMPLE[2], "due to incomplete multiplexing")
        d$SAMPLE <- d$SAMPLE[1]
        d$NAME <- d$NAME[1]
        d$LIBRARY_SOURCE <- d$LIBRARY_SOURCE[1]
      }
      if (length(d$NAME) > 1) print(d)
      sink(file.path("output", "ENA", paste0(d$NAME, ".manifest")), append = FALSE)
      iwalk(d, ~ for (x in .x) {cat(.y, x); cat("\n")})
      sink()
    }
  )

dplyr::inner_join(
  dplyr::bind_rows(
    readd(soil_samples),
    readd(control_samples)
  ) %>% select(sample_alias, "project name", plate, well, amplicon, sample_type),
  readd(illumina_group_SH.2257),
  by = c("plate", "well", "amplicon")
) %>%
  dplyr::transmute(
    STUDY = `project name`,
    SAMPLE = sample_alias,
    NAME = paste(sample_alias, tech, direction, sep = "-") %>% str_replace("-$", "") %>% str_replace(" ", ""),
    INSTRUMENT = paste(tech, machine),
    LIBRARY_SOURCE = ifelse(sample_type == "Pos", "GENOMIC", "METAGENOMIC"),
    LIBRARY_SELECTION = "PCR",
    LIBRARY_STRATEGY = "AMPLICON",
    INSERT_SIZE = "350",
    FASTQ = trim_file_R1,
    FASTQ2 = trim_file_R2,
    DESCRIPTION = case_when(
      direction == "f" ~ "Reads in forward orientation",
      direction == "r" ~ "Reads in reverse orientation",
      TRUE ~ "")
  ) %>%
  slice(1) %>%
  rowwise() %>%
  group_split() %>%
  walk(
    function(d) {
      d <- as.list(d)
      d$FASTQ <- c(d$FASTQ, d$FASTQ2)
      d$FASTQ2 <- NULL
      sink(file.path("output", "ENA", paste0(d$NAME, ".manifest")), append = FALSE)
      iwalk(d, ~ for (x in .x) {cat(.y, x); cat("\n")})
      sink()
    }
  )


get_reads_accno <- function(sample, reads_dir = "output/ENA/reads") {
  file.path(reads_dir, sample, "submit", "receipt.xml") %>%
    xml2::read_xml() %>%
    xml2::xml_find_first("RUN") %>%
    xml2::xml_attr("accession")
}

list.files(path = "output/ENA/reads", pattern = "receipt.xml", recursive = TRUE) %>%
  stringr::str_split_fixed("/", 3) %>%
  magrittr::extract( , 1) %>%
  set_names(., .) %>%
  purrr::map_chr(get_reads_accno) %>%
  tibble::enframe(name = "sample_alias", value = "accno") %>%
  tidyr::extract(
    "sample_alias",
    into = c("sample_alias", "tech"),
    regex = "(.+)-((?:Illumina|PacBio|IonTorrent)-?[fr]?)"
  ) %>%
  tidyr::pivot_wider(names_from = "tech", values_from = "accno") %>% View()

get_reads_accno("OT-control-001-L-PacBio-r")
get_reads_md5 <- function(sample, reads_dir = "output/ENA/reads") {
  file.path(reads_dir, sample, "submit", "run.xml") %>%
    xml2::read_xml() %>%
    xml2::xml_find_first("//FILE") %>%
    xml2::xml_attr("checksum")
}

ncbi_get_taxon_summary()%>%
  stringi::stri_replace_all_fixed(
    c("Fungi", "Alveolata"),
    c("fungus", "alveolate"),
    vectorize_all = FALSE
  ) %>%
  paste("uncultured", .) %T>%
  writeLines("output/ENA/env_taxa.txt")

drake_cache()$list() %>% purrr::keep(stringr::str_detect, "physeq")
loadd(asv_table)
pooled_inext <- asv_table %>%
  dplyr::mutate_at("seq", tzara::seqhash) %>%
  tibble::column_to_rownames("seq") %>%
  as.matrix() %>%
  subtract(1) %>%
  pmax(0) %>%
  iNEXT::iNEXT()

sample_inext <-
  readd(big_seq_table_ITS2) %>%
  magrittr::extract(rowSums(.) > 100,) %>%
  t() %>%
  iNEXT::iNEXT(knots = 200)

alpha100 <- readd(big_seq_table_ITS2) %>%
  magrittr::extract(rowSums(.) > 100,) %>%
  t() %>%
  {
    dplyr::bind_cols(
      list(sample = rownames(.)),
      iNEXT::estimateD(., level = 100) %>%
        dplyr::filter(order == 0)
    )
  }


  only_obs <- function(d) {
    dplyr::filter(d$data, method == "observed")
}

library(ggplot2)
sample_inext %>%
  ggplot2::fortify() %>%
  tidyr::extract(
    site,
    into = c("seq_run", "t", "plate", "well", "region"),
    regex = "((is|pb|SH)[-_][0-9]{3,4})_([^_]+)_([A-H]1?[0-9])_(.+)",
    remove = FALSE
  ) %>%
  dplyr::filter(method != "extrapolated", seq_run != "is_057") %>%
  ggplot(aes(x = x, y = y, color = seq_run)) +
  # lines for each sample
  geom_line(aes(group = site), alpha = 0.08) +
  # use the alpha values I provided as the literal alpha values
  scale_alpha_identity() +
  # small points to represent the actual observed read depth and richness
  geom_point(data = ~dplyr::filter(., method == "observed"), alpha = 0.8, size = 1, shape = 1) +
  # seperate facet for each forest type.
  # They are side-by-side so we can compare them visually.
  # facet_wrap(facets = ~seq_run, nrow = 1) +
  xlab("reads") +
  ylab("Species richness") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

source(here::here("scripts", "singleton_estimator.R"))
