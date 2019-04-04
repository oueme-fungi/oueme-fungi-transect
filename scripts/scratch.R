library(drake)
library(dplyr)
library(tidyverse)
library(magrittr)
library(Biostrings)
library(rgbif)
# try to combine ASVs based on co-occurence
# problem : need to process co-occurence groups in order of reads, not by joining the first group first.??

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

tax_silva <- DECIPHER::IdTaxa(testLSU, classifier_silva_LSU, strand = "top", processors = 3, threshold = 30)
tax_rdp <- DECIPHER::IdTaxa(testLSU, classifier_rdp_LSU, strand = "top", processors = 3, threshold = 50)
tax_unite <- DECIPHER::IdTaxa(testITS, classifier_unite_ITS2, strand = "top", processors = 3, threshold = 30)

tax_rdp2 <- dada2::assignTaxonomy(as.character(Biostrings::DNAStringSet(testLSU)), "reference/rdp.LSU.dada2.fasta.gz")
tax_unite2 <- dada2::assignTaxonomy(as.character(Biostrings::DNAStringSet(testITS)), "reference/unite.ITS2.fasta.gz")

tax_dada_unite <- dada2::assignTaxonomy(testITS, unite_file)

inocybe_rdp <- dplyr::filter(rdp_c12n, grepl("Inocybe", c12n))
unique(inocybe_rdp$c12n)

penicillium_rdp <- dplyr::filter(rdp_c12n, grepl("Penicillium", c12n))
unique(penicillium_rdp$species)

silva_db <- unique(Biostrings::readRNAStringSet(silva_file))
silva_LR5 <- Biostrings::vmatchPattern(pattern = "CGAAGUUUCCCUCAGGA", silva_db, 
                                       max.mismatch = 2)
whichLR5 <- XVector::elementNROWS(silva_LR5) > 0
silva_LR5_db <- silva_db[whichLR5]
silva_LR5 <- silva_LR5[whichLR5]
silva_LR5_db <- Biostrings::subseq(silva_LR5_db, end = purrr::map_int(startIndex(silva_LR5), dplyr::first))

dada2:::

rdp_LR5 <- Biostrings::vmatchPattern(pattern = "CGAAGTTTCCCTCAGGA", rdp_db, 
                                       max.mismatch = 2)
whichLR5 <- XVector::elementNROWS(rdp_LR5) > 0
rdp_LR5_db <- rdp_db[whichLR5]
rdp_LR5 <- rdp_LR5[whichLR5]
rdp_LR5_db <- Biostrings::subseq(rdp_LR5_db, end = purrr::map_int(startIndex(rdp_LR5), dplyr::first))


rdp_LF402 <- Biostrings::vmatchPattern(pattern = "TGAAAT", rdp_LR5_db, 
                                     max.mismatch = 2)
which_LF402 <- XVector::elementNROWS(rdp_LF402) > 0
sum(which_LF402)
rdp_LR5_db[grepl("Inocybe", names(rdp_LR5_db))]

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

  asvtab <-  %>%
    phyloseq::otu_table(taxa_are_rows = FALSE)
}
      


st <- 

longtree <- ape::read.tree("data/pasta/pasta_raxml.tree")

phyloseq::taxa_names(longtree) <-
  tibble::tibble(hash = phyloseq::taxa_names(longtree)) %>%
  dplyr::left_join(conseq_filt) %$%
  seqhash(ITS2)

physeq <- phyloseq::phyloseq(pm, st, longtree) %>%
  phyloseq::prune_samples(samples = phyloseq::sample_data(.)$Dataset == "long-pacbio"
                          & phyloseq::sample_data(.)$Qual == ""
                          & phyloseq::sample_data(.)$sample_type == "Sample")

unif_dist <- phyloseq::distance(physeq, "unifrac")
bc_dist <- phyloseq::distance(physeq, "bray")
spatial_dist <- phyloseq::sample_data(physeq) %>%
  with(X + 30000 * as.integer(Site) + 100000 * as.integer(Year)) %>%
  dist
which(is.na(unif_dist))
which(is.na(bc_dist))
which(is.na(spatial_dist))
unif_corr <- vegan::mantel.correlog(unif_dist, spatial_dist,
                       break.pts = 0:12 + 0.5,
                       cutoff = FALSE)
bc_corr <- vegan::mantel.correlog(bc_dist, spatial_dist,
                       break.pts = 0:12 + 0.5,
                       cutoff = FALSE)
plot(unif_corr)
plot(bc_corr)

bc_fit <- lm(log(Mantel.cor) ~ class.index, data = as.data.frame(bc_corr$mantel.res))
unif_fit <- lm(log(Mantel.cor) ~ class.index, data = as.data.frame(unif_corr$mantel.res))

as.data.frame(bc_corr$mantel.res) %>%
  ggplot(aes(class.index, Mantel.cor)) +
  geom_point() +
  geom_line(aes(y = exp(predict(bc_fit)))) +
  scale_y_continuous(limits = c(0, NA), name = "Bray-Curtis -- Spatial Distance correlation") +
  xlab("distance (m)")

as.data.frame(unif_corr$mantel.res) %>%
  ggplot(aes(class.index, Mantel.cor)) +
  geom_point() +
  geom_line(aes(y = exp(predict(unif_fit)))) +
  scale_y_continuous(limits = c(0, NA), name = "Bray-Curtis -- Spatial Distance correlation") +
  xlab("distance (m)")

# phyloseq::distanceMethodList
  left_join(select(datasets, Dataset, Tech, Seq.Run, Forward, Reverse)) %>%
  mutate(Primer.Pair = paste(str_replace(Forward, "_tag.*$", ""),
                             str_replace(Reverse, "_tag.*$", ""),
                             sep = "_")) %>%
  select(Seq.Run, Plate, Well, Region, Primer.Pair, everything()) %>%
  mutate_at(c("Seq.Run", "Region", "Primer.Pair"), factor) %>%
  mutate_at("Well", factor, levels = levels(pm$Well)) %>%
  mutate_at("Plate", factor, levels = levels(pm$Plate))


%>%
  tidyr::gather(key = "ASV", value = "Reads",
         -one_of(c("Seq.Run", "Plate", "Well", "Region", "Primer.Pair", "Dataset", "Tech", "Runs", "Forward", "Reverse", "Regions", "PlateKey", "Homopolymer.Gap.Penalty", "Band.Size", "Pool"))) %>%
  mutate_at("Reads", as.integer) %>%
  filter(Reads > 0L)

guilds <- readd(guilds_table_ITS2_unite) %>%
  mutate(hash = map_chr(seq, digest::digest) %>% stringr::str_sub(end = 8)) %>%
  mutate_at(c("trophicMode", "guild"), stringr::str_split, "-")

ecm_taxa <- filter(guilds, map_lgl(guild, magrittr::is_in, x = "Ectomycorrhizal")) %$% hash

samples_to_use <- 
  right_join(pm, st) %>%
    filter(!is.na(Year), Qual == "", Dataset == "short-pacbio") %>%
    group_by(Dataset, X, Site) %>%
    summarize(Reads = sum(Reads)) %>%
  filter(Reads > 0) %>%
  select(-Reads)


testdata <- 
  right_join(pm, st) %>%
  filter(!is.na(Year), Qual == "") %>%
  semi_join(samples_to_use) %>%
  select(Site, X, ASV, Reads) %>%
  group_by(Site, X, ASV) %>%
  summarize(Reads = sum(Reads)) %>%
  unite("Sample", Site, X, remove = FALSE) %>%
  arrange(Sample)

sitecounts <- testdata %>%
  group_by(Sample, X, Site, .drop = TRUE) %>%
  summarize(Reads = sum(Reads))

library(ggplot2)
ggplot(sitecounts, aes(x = X, y = Reads)) +
  geom_point() +
  facet_wrap(~Site, ncol = 1)

countdist <- sitecounts %$%
  set_names(Reads, Sample) %>%
  log %>%
  stats::dist()

comm <- testdata %>%
  ungroup() %>%
  tidyr::spread(key = ASV, value = Reads, fill = 0) %>%
  dplyr::select(-Site, -X) %>%
  column_to_rownames("Sample") %>%
  as.matrix() %>%
  vegan::decostand("total")

ecm_comm <- comm %>%
  magrittr::extract(,ecm_taxa[ecm_taxa %in% colnames(.)]) %>%
  vegan::decostand("total")

ecodist <- vegan::vegdist(comm, "bray")
ecodist_ecm <- vegan::vegdist(ecm_comm, "bray")

spatialdist <- testdata %>%
  ungroup() %>%
  select(Sample, Site, X) %>%
  unique() %>%
  tidyr::crossing(., .) %>%
  mutate(dist = ifelse(Site == Site1,
                       abs(X - X1),
                       30000)) %>%
  select(Sample, Sample1, dist) %>%
  tidyr::spread(key = Sample1, value = dist) %>%
  column_to_rownames("Sample") %>%
  as.matrix() %>%
  as.dist()

corgram <- vegan::mantel.correlog(D.eco = ecodist, D.geo = spatialdist, break.pts = 0:12 + 0.5, cutoff = FALSE)

plot(corgram)
corgram

corgram_ecm <- vegan::mantel.correlog(D.eco = ecodist_ecm, D.geo = spatialdist, break.pts = 0:12 + 0.5, cutoff = FALSE)

plot(corgram_ecm)
corgram_ecm

corgram2 <- vegan::mantel.correlog(D.eco = countdist,
                                   D.geo = spatialdist,
                                   break.pts = 0:12 + 0.5,
                                   cutoff = FALSE)
plot(corgram2)
corgram2
