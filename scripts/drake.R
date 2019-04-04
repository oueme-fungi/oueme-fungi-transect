if (interactive()) {
  cat("Creating drake plan in interactive session...\n")
  library(here)
  base.dir <- here()
  r.dir <- here("scripts")
  data.dir <- here("data")
  lab.dir <- here("config")
  seq.dir <- here("sequences")
  trim.dir <- file.path(seq.dir, "trim")
  region.dir <- file.path(seq.dir, "regions")
  filter.dir <- file.path(seq.dir, "filter")
  pasta.dir <- file.path(data.dir, "pasta")
  plan.dir <- file.path(data.dir, "plan")
  ref.dir <- here("reference")
  rmd.dir <- here("writing")
  out.dir <- here("output")
  in.files <- list.files(path = trim.dir,
                         pattern = "\\.trim\\.fastq\\.gz$",
                         full.names = TRUE, recursive = TRUE)
  dataset.file <- file.path(lab.dir, "datasets.csv")
  regions.file <- file.path(lab.dir, "regions.csv")
  platemap.file <- file.path(lab.dir, "Brendan_soil2.xlsx")
  platemap.sheet <- "Concentration samples"
  rdp_file <- file.path(ref.dir, "rdp.LSU.fasta.gz")
  silva_file <- file.path(ref.dir, "silva.LSU.fasta.gz")
  silva_tax_file <- file.path(ref.dir, "silva_tax.txt")
  unite_file <- file.path(ref.dir, "unite.fasta.gz")
  unite_patch_file <- file.path(ref.dir, "unite_patch.csv")
  rdp_patch_file <- file.path(ref.dir, "rdp_patch.csv")
  silva_patch_file <- file.path(ref.dir, "silva_patch.csv")
  plan_file <- file.path(plan.dir, "plan.rds")
  meta1_file <- file.path(plan.dir, "meta1.rds")
  meta2_file <- file.path(plan.dir, "meta2.rds")
  meta3_file <- file.path(plan.dir, "meta3.rds")
  meta4_file <- file.path(plan.dir, "meta4.rds")
  meta4csv_file <- file.path(plan.dir, "meta4.csv")
  tid_file <- file.path(plan.dir, "tids.txt")
  drakedata.file <- file.path(plan.dir, "drake.Rdata")
  bigsplit <- 80L
} else if (exists("snakemake")) {
  cat("Creating drake plan in snakemake session...\n")
  snakemake@source(".Rprofile", echo = FALSE)
  r.dir <- snakemake@config$rdir
  dataset.file <- snakemake@input$dataset
  regions.file <- snakemake@input$regions
  platemap.file <- snakemake@input$platemap
  platemap.sheet <- snakemake@config$platemap_sheet
  target <- sub("^.", "", snakemake@output)
  trim.dir <- snakemake@config$trimdir
  region.dir <- snakemake@config$regiondir
  filter.dir <- snakemake@config$filterdir
  ref.dir <- snakemake@config$ref_root
  rdp_file <- snakemake@config$rdp_file
  rdp_patch_file <- snakemake@config$rdp_patch_file
  silva_file <- snakemake@config$silva_file
  silva_tax_file <- snakemake@config$silva_tax_file
  silva_patch_file <- snakemake@config$silva_patch_file
  unite_file <- snakemake@config$unite_file
  unite_patch_file <- snakemake@config$unite_patch_file
  rmd.dir <- snakemake@config$rmddir
  out.dir <- snakemake@config$outdir
  pasta.dir <- snakemake@config$pastadir
  plan.dir <- snakemake@config$plandir
  plan_file <- snakemake@output$plan
  meta1_file <- snakemake@output$meta1
  meta2_file <- snakemake@output$meta2
  meta3_file <- snakemake@output$meta3
  meta4_file <- snakemake@output$meta4
  drakedata_file <- snakemake@output$drakedata
  meta4csv_file <- snakemake@output$meta4csv
  tids_file <- snakemake@output$tids
  prereqs <- unlist(snakemake@input)
  in.files <- stringr::str_subset(prereqs,
                         pattern = "\\.trim\\.fastq\\.gz$")
  bigsplit <- snakemake@config$bigsplit
  plan.file <- snakemake@output[["plan"]]
  meta1.file <- snakemake@output[["meta1"]]
  meta2.file <- snakemake@output[["meta2"]]
  meta3.file <- snakemake@output[["meta3"]]
  meta4.file <- snakemake@output[["meta4"]]
  meta4.csv <- snakemake@output[["meta4.csv"]]
  pid.file <- snakemake@output[["pids"]]
  tid.file <- snakemake@output[["tids"]]
  drakedata.file <- snakemake@output[["drakedata"]]
  logfile <- file(snakemake@log[[1]], open = "at")
  sink(logfile)
  sink(logfile, type = "message")
} else {
  stop("Can't find Snakemake object in non-interactive session!")
}

library(magrittr)
library(tidyverse)
library(rlang)
library(glue)
library(drake)
library(assertr)


#### read scripts and configs ####
# load various pipeline functions
source(file.path(r.dir, "parallel_helpers.R"))
source(file.path(r.dir, "combine_derep.R"))
source(file.path(r.dir, "extract_regions.R"))
source(file.path(r.dir, "dada.R"))
source(file.path(r.dir, "taxonomy.R"))
source(file.path(r.dir, "plate_check.R"))
source(file.path(r.dir, "map_LSU.R"))

# load the dataset and region definitions
datasets <- read_csv(dataset.file)
regions <- read_csv(regions.file)

# "meta" dataframes containing definitions for each instance of mapped targets

#### meta1 ####
# meta1 has one row per demultiplexed, primer-trimmed fastq.gz file
# The "index" is FileID
meta1 <- datasets %>%
  mutate(PrimerID = paste0(str_replace(Forward, "_tag.*$", ""),
                            str_replace(Reverse, "_tag.*$", "")),
         Plate = map(Runs, ~formatC(seq(.), width = 3, flag = "0"))) %>%
  unnest(Plate) %>%
  mutate(Direction = if_else(Tech == "PacBio",
                             list(c("f", "r")),
                             list(""))) %>%
  unnest(Direction) %>%
  mutate(Well = list(tidyr::crossing(Row = LETTERS[1:8], Col = 1:12) %>%
                       transmute(Well = paste0(Row, Col)))) %>%
  unnest(Well) %>%
  mutate(Trim.File = glue("{Seq.Run}_{Plate}/{Seq.Run}_{Plate}-{Well}{Direction}.trim.fastq.gz"),
         SeqRunID = Seq.Run,
         
         FileID = glue("{Seq.Run}{Plate}{Well}{Direction}") %>%
           str_replace_all("[:punct:]", ""),
         positions = make.names(glue("positions_{FileID}")),
         join_derep = make.names(glue("join_derep_{PrimerID}")),
         itsxtrim = make.names(glue("itsxtrim_{PrimerID}"))) %>%
  filter(file.exists(file.path(trim.dir, Trim.File))) %>%
  filter(file.size(file.path(trim.dir, Trim.File)) > 40) %>%
  verify(file.path(trim.dir, Trim.File) %in% in.files) %>%
  mutate_at(c("join_derep", "itsxtrim", "FileID", "SeqRunID", "PrimerID"), syms)

#### meta2 ####
# meta2 has one row per region per input file
# The index is WellID + RegionID
# PlateID is defined because it will be the index for meta3
meta2 <- meta1 %>%
  mutate_at("Regions", str_split, ",") %>%
  unnest(Region = Regions, .preserve = ends_with("ID")) %>%
  left_join(regions, by = "Region") %>%
  mutate(WellID = glue("{Seq.Run}{Plate}{Well}{Direction}"),
         PlateID = glue("{Seq.Run}_{Plate}"), 
         RegionID = Region,
         Region.File = glue("{Seq.Run}_{Plate}-{Well}{Direction}-{Region}.trim.fastq.gz"),
         Filter.File = glue("{Seq.Run}_{Plate}-{Well}{Direction}-{Region}.qfilt.fastq.gz")) %>%
  mutate_at(c("positions", "WellID", "PlateID", "RegionID"), syms)

if (interactive()) {
  #use a subset
  meta1 %<>%
    group_by(Seq.Run, Plate) %>%
    mutate(filesize = file.size(file.path(trim.dir, Trim.File))) %>%
    arrange(desc(filesize), .by_group = TRUE) %>%
    dplyr::slice(1:4) %>%
    ungroup()
  
  meta2 %<>%
    semi_join(meta1, by = "Trim.File")
} else {
  saveRDS(meta1, meta1_file)
  saveRDS(meta2, meta2_file)
}

#### meta3 ####
# meta3 has one row per region per plate
# The index is PlateID + RegionID (Plate ID)
meta3 <- meta2 %>%
  select(Seq.Run, Plate, Region, PlateID, RegionID, SeqRunID) %>%
  unique() %>%
  arrange(Seq.Run, Plate, Region) %>%
  mutate(region_derep = glue("region_derep_{PlateID}_{RegionID}")) %>%
  left_join(datasets %>% select(-Regions), by = "Seq.Run") %>%
  left_join(regions, by = "Region") %>%
  mutate_at(c("region_derep"), syms)
if (!interactive()) {
  saveRDS(meta3, meta3_file)
}

#### meta4 ####
# meta4 has one row per region
# The index is TaxID (Taxonomy ID)
meta4 <- meta3 %>%
  mutate_if(is.list, as.character) %>%
  group_by_at(names(regions)) %>%
  summarize(PlateRegionID = paste(PlateID, RegionID, sep = "_", collapse = ",")) %>%
  ungroup() %>%
  mutate(Reference = strsplit(Reference, " *, *") %>%
                                map(unlist)) %>%
  unnest(Reference, .preserve = Region) %>%
  mutate(TaxID = glue("{Region}_{Reference}"),
         big_seq_table = glue("big_seq_table_{Region}"),
         Reference.File = Reference) %>%
  arrange(TaxID) %>%
  mutate_at(c("TaxID", "big_seq_table"), syms)
if (!interactive()) {
  saveRDS(meta4, meta4_file)
  mutate_if(meta4, is.list, as.character) %>%
  readr::write_csv(meta4csv_file)
  readr::write_lines(meta4$TaxID, tids_file)
}

#### drake plan ####

cat("\nbuilding plan...\n")
tictoc::tic()
plan <- drake_plan(
  
  datasets = readr::read_csv(file_in(!!dataset.file)),
  
  # derep1 ----
  # dereplicate input sequences before running ITSx
  # this is done independently in parallel
  derep1 = target(
    tibble::tibble(
      file = Trim.File,
      derep = list(dada2::derepFastq(file_in(!!file.path(trim.dir, Trim.File)),
                                     qualityType = "FastqQuality",
                                     n = 1e4) %>% inset2("quals", NULL))),
    transform = map(.data = !!select(meta1, Trim.File, FileID, PrimerID),
                    .tag_in = step,
                    .id = FileID)),
  
  # join_derep----
  # join the results from dereplicating into a list of unique sequences
  # and a map from the original sequences to the uniques
  join_derep = target(
    combine_derep(derep1),
    transform = combine(derep1,
                        .tag_in = step,
                        .by = PrimerID)),
  # split_fasta----
  # break the list of unique files into equal sized chunks for ITSx
  split_fasta = target(
    split(join_derep$fasta, seq_along(join_derep$fasta) %% !!bigsplit + 1),
    transform = map(join_derep,
                    .tag_in = step)),
  
  # itsx_shard----
  # find the location of LSU, 5.8S, and SSU (and this ITS1 and ITS2)
  # in the sequences from each chunk
  itsx_shard = target({
    library(Biostrings)
    library(dplyr)
    rITSx::itsx(in.file = split_fasta[[shard]],
                # nhmmer = TRUE,
                read_function = Biostrings::readDNAStringSet,
                fasta = FALSE, summary = FALSE, graphical = FALSE,
                save_regions = "none", not_found = FALSE,
                complement = FALSE, cpu = ignore(itsx_cpus))[["positions"]]
    },
    transform = cross(split_fasta,
                      shard = !!(1:bigsplit),
                      .tag_in = step,
                      .id = PrimerID)),
  
  # itsxtrim ----
  # combine the chunks into a master positions list.
  itsxtrim = target(
    dplyr::bind_rows(itsx_shard),
    transform = combine(itsx_shard,
                        .tag_in = step,
                        .by = PrimerID)),
  
  # positions----
  # find the entries in the positions list which correspond to the sequences
  # in each file
  positions = target(
    dplyr::filter(join_derep$map, file == Trim.File) %>%
      dplyr::left_join(itsxtrim, by = c("newmap" = "seq")),
    transform = map(.data = !!select(meta1, join_derep, itsxtrim, Trim.File,
                                     FileID),
                    .tag_in = step,
                    .id = FileID)),
  
  # regions----
  # cut the required regions out of each sequence
  regions = target(
    extract_region(infile = file_in(!!file.path(trim.dir, Trim.File)),
                   outfile = file_out(!!file.path(region.dir, Region.File)),
                   region = Region,
                   positions = positions),
    transform = map(.data = !!select(meta2, Seq.Run, Trim.File, Region.File, Region,
                                     positions, WellID, PlateID, RegionID, SeqRunID),
                    .tag_in = step,
                    .id = c(WellID, RegionID))),
  
  # filter----
  # quality filter the sequences
  filter = target(
    {
      outfile <- file_out(!!file.path(filter.dir, Filter.File))
      unlink(outfile)
      out <- dada2::filterAndTrim(
        fwd = file_in(!!file.path(region.dir, Region.File)),
        filt = file_out(!!file.path(filter.dir, Filter.File)),
        maxLen = MaxLength,
        minLen = MinLength,
        maxEE = MaxEE,
        qualityType = "FastqQuality")
      if (!file.exists(outfile)) {
        ShortRead::writeFastq(ShortRead::ShortReadQ(),
                              outfile)
      }
      return(out)
    },
    transform = map(.data = !!select(meta2, Region.File, Filter.File,
                                     MaxLength, MinLength, MaxEE, WellID, RegionID),
                    .tag_in = step,
                    .id = c(WellID, RegionID))),
  
  # derep2----
  # Do a second round of dereplication on the filtered regions.
  derep2 = target({
    infile <- file_in(!!file.path(filter.dir, Filter.File))
    if (file.size(infile) > 50) {
      out <- dada2::derepFastq(fls = infile, n = 1e4,
                               qualityType = "FastqQuality")
      fqs <- ShortRead::FastqStreamer(infile, n = 1e4)
      seqids <- character(0)
      while (length(fq <- ShortRead::yield(fqs, qualityType = "FastqQuality"))) {
        seqids <- c(seqids, as.character(fq@id))
      }
      close(fqs)
      names(out[["map"]]) <- seqids
    } else {
      out = NULL
    }
    result <- list()
    result[[Filter.File]] <- out
    result
  },
  transform = map(.data = !!meta2,
                  .tag_in = step,
                  .id = c(WellID, RegionID))),
  
  # region_derep ----
  # put together all the dereplicated files for each region for each plate.
  region_derep = target(
    purrr::compact(c(derep2)),
    transform = combine(derep2, .by = c(PlateID, RegionID), .tag_in = step),
  ),
  
  # err----
  # Calibrate dada error models
  err = target({
    err.fun <- ifelse(Tech == "PacBio",
                      dada2::PacBioErrfun,
                      dada2::loessErrfun)
    dada2::learnErrors(
      fls = region_derep, nbases = 1e8,
      multithread = ignore(dada_cpus), randomize = TRUE,
      errorEstimationFunction = err.fun, 
      HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(Homopolymer.Gap.Penalty)),
      BAND_SIZE = Band.Size,
      pool = eval(rlang::parse_expr(Pool)),
      verbose = TRUE,
      qualityType = "FastqQuality")},
    transform = map(.data = !!meta3,
                    .tag_in = step, .id = c(PlateID, RegionID))),
  
  # dada----
  # Run dada denoising algorithm
  dada = target(
    dada2::dada(
      derep = region_derep, err = err,
      multithread = ignore(dada_cpus), 
      HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(Homopolymer.Gap.Penalty)),
      BAND_SIZE = Band.Size,
      pool = eval(rlang::parse_expr(Pool))),
    transform = map(err, .data = !!meta3,
                    .tag_in = step,
                    .id = c(PlateID, RegionID))),
  
  # dadamap----
  # Make maps from individual sequences in source files to dada ASVs
  dada_map = target(
    dadamap(region_derep, dada),
    transform = map(dada, .data = !!meta3,
                    .tag_in = step, .id = c(PlateID, RegionID))),
  
  # seq_table----
  # Make a sample x ASV abundance matrix
  seq_table = target(
    dada2::makeSequenceTable(dada),
    transform = map(dada, .data = !!meta3, .tag_in = step, .id = c(PlateID, RegionID))),
  
  # nochim----
  # Remove likely bimeras
  nochim = target(
    dada2::removeBimeraDenovo(seq_table, method = "consensus",
                              multithread = ignore(dada_cpus)),
    transform = map(seq_table, .data = !!meta3, .tag_in = step, .id = c(PlateID, RegionID))),
  
  # big_seq_table ----
  # Join all the sequence tables for each region
  big_seq_table = target(
    dada2::mergeSequenceTables(tables = list(nochim)),
    transform = combine(nochim, .tag_in = step, .by = RegionID)),
  
  # dbprep ----
  # import the RDP, Silva, and UNITE databases, and format them for use with
  # DECIPHER.
  dbprep = target(
    formatdb(db = Reference,
             seq_file = file_in(seq_file),
             patch_file = file_in(patch_file),
             tax_file = file_in(tax_file),
             maxGroupSize = 20),
    
    transform = map(Reference = c("rdp", "silva"),
                    ReferenceID = !!rlang::syms(c("rdp", "silva")),
                    seq_file = !!c(rdp_file, silva_file),
                    patch_file = !!c(rdp_patch_file, silva_patch_file),
                    tax_file = !!c(NA, silva_tax_file),
                    .id = ReferenceID)),
  # classifier----
  # Train a DECIPHER classifier for each database.
  classifier = target(
    refine_classifier(dbprep$train,
                      dbprep$taxonomy,
                      dbprep$rank,
                      out_root = !!file.path(ref.dir, Reference)),
    transform = map(dbprep, .id = ReferenceID)
  ),
  
  # taxon ----
  # Assign taxonomy to each ASV
  taxon = target(
    taxonomy(big_seq_table,
             reference = file_in(!!file.path(ref.dir, paste0(Reference, ".fasta.gz"))),
             multithread = ignore(dada_cpus)),
    transform = map(.data = !!meta4, .tag_in = step, .id = TaxID)),
  
  # funguild_db ----
  # Download the FUNGuild database
  funguild_db = FUNGuildR::get_funguild_db(),
  
  # guilds_table ----
  # Assign ecological guilds to the ASVs based on taxonomy.
  guilds_table = target(
    FUNGuildR::funguild_assign(taxon, db = funguild_db),
    transform = map(taxon, TaxID, .tag_in = step, .id = TaxID)),
  
  # platemap ----
  # Read the map between plate locations and samples
  platemap = read_platemap(file_in(!!platemap.file), platemap.sheet),
  
  # raw ----
  # Join the raw read info (for long pacbio reads).
  raw = target(
    raw_reads(regions, 
              filenames = purrr::map_chr(rlang::exprs(regions), rlang::as_string),
              max_ee = 3),
    transform = combine(regions, Seq.Run, Region,
                        .by = c(SeqRunID, RegionID),
                        .tag_in = step,
                        .id = c(SeqRunID, RegionID))
  ),
  
  # combined ----
  # Replace raw reads which were mapped to a DADA2 ASV with the ASV,
  # and put them all in one tibble.  We'll only do this for pacbio long reads.
  combined = target(
    combine_bigmaps(list(dada_map), dplyr::bind_rows(raw)),
    transform = combine(dada_map, raw, .by = SeqRunID, .tag_in = step)
  ),
  
  # conseq ----
  # Calculate consensus sequences for full ITS and LSU for each ITS2 ASV with
  # with at least 3 sequences.
  conseq = 
    combined_pb_500 %>%
    dplyr::group_by(ITS2) %>%
    dplyr::filter(!is.na(ITS2), dplyr::n() >= 3) %>%
    dplyr::summarize(nreads = dplyr::n(),
                     long = as.character(
                       calculate_consensus(long, seq.id, ignore(dada_cpus))),
                     ITS = as.character(
                       calculate_consensus(ITS, seq.id, ignore(dada_cpus))),
                     LSU = as.character(
                       calculate_consensus(LSU, seq.id, ignore(dada_cpus)))),
  # conseq_filt ----
  # Remove consensus sequences which did not have a strong consensus.
  conseq_filt =
    conseq %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(
      qlong = Biostrings::RNAStringSet(long) %>%
        Biostrings::letterFrequency(., "MRWSYKVHDBN"),
      qITS = Biostrings::RNAStringSet(ITS) %>%
        Biostrings::letterFrequency(., "MRWSYKVHDBN"),
      qLSU = Biostrings::RNAStringSet(LSU) %>%
        Biostrings::letterFrequency(., "MRWSYKVHDBN")) %>%
    dplyr::filter(qITS < 3, qLSU < 3, qlong < 3) %>%
    dplyr::mutate(hash = seqhash(long)) %>%
    dplyr::arrange(hash),
  
  # taxon_LSUcons_rdp----
  # Assign taxonomy to the consensus LSU sequences.
  taxon_LSUcons_rdp =
    conseq_filt$LSU %>%
    unique() %>%
    stringr::str_replace_all("U", "T") %>%
    taxonomy(reference = file_in("reference/rdp.fasta.gz"),
             multithread = ignore(dada_cpus)) %>%
    dplyr::mutate_at("seq", stringr::str_replace_all, "T", "U"),
  
  # cons_tax----
  # Make names which combine the ITS2 and LSU identifications to put on the tree
  # Also make hash names, because PASTA will destroy non-alphanumeric names
  cons_tax  =
    conseq_filt %>%
    dplyr::left_join(dplyr::select(taxon_ITS2_unite, ITS2 = seq, Name),
                     by = "ITS2") %>%
    dplyr::left_join(dplyr::select(taxon_LSUcons_rdp, LSU = seq, Name),
                     suffix = c("_ITS2", "_LSU"),
                     by = "LSU") %>%
    tidyr::unite("Name", Name_ITS2, Name_LSU, sep = "/"),
  
  # long_consensus----
  # get the long amplicon consensus and convert to RNAStringSet.
  # Use the sequence hashes as the name, so that names will be robust in PASTA
  long_consensus = conseq_filt %$%
    rlang::set_names(long, hash) %>%
    Biostrings::RNAStringSet(),
  
  # lsualn----
  # Align the LSU consensus.
  # DECIPHER has a fast progressive alignment algorthm that calculates RNA
  # secondary structure
  lsualn =
    DECIPHER::AlignSeqs(long_consensus,
                        iterations = 10, refinements = 10,
                        processors = ignore(dada_cpus)),
  
  # write_lsualn----
  # write the aligned consensus for PASTA.
  # also touch a file to mark that we really need to re-run PASTA.
  write_lsualn = {
    Biostrings::writeXStringSet(lsualn,
                                file_out(ignore(longASV_file)))
    Sys.setFileTime(file.path(pasta.dir, ".dopasta"), Sys.time())
  },
  
  # read the long amplicon tree in from 
  
  # seq_count ----
  # Count the sequences in each fastq file
  seq_count = target(
    tibble::tibble(
      file = file_in(fq.file),
      reads = system(glue::glue("zcat {file} | grep -c '^@'", file = fq.file),
                     intern = TRUE) %>%
        as.integer),
    transform = map(fq.file = !!c(file.path(trim.dir, meta1$Trim.File),
                                 file.path(region.dir, meta2$Region.File),
                                 file.path(filter.dir, meta2$Filter.File)),
                    .tag_in = step)),
  
  # seq_counts----
  # merge all the sequence counts into one data.frame
  seq_counts = target(
    dplyr::bind_rows(seq_count),
    transform = combine(seq_count)),

  # qstats----
  # calculate quality stats for the different regions
  qstats = target(
    q_stats(file_in(!!file.path(region.dir, Region.File))),
    transform = map(.data = !!meta2,
                    .tag_in = step,
                    .id = c(WellID, RegionID))),
  #qstats_join----
  # join all the quality stats into one data.frame
  qstats_join = target(
    dplyr::bind_rows(qstats) %>%
      tidyr::extract(
        col = "file",
        into = c("Seq.Run", "Plate", "Well", "Direction", "Region"),
        regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.trim\\.fastq\\.gz") %>%
      dplyr::left_join(datasets),
    transform = combine(qstats)),
  # qstats_knit----
  # knit a report about the quality stats.
  qstats_knit = {
    if (!dir.exists(!!out.dir)) dir.create(!!out.dir, recursive = TRUE)
    rmarkdown::render(
      knitr_in(!!file.path(rmd.dir, "qual-check.Rmd")),
      output_file = file_out(!!file.path(out.dir, "qual-check.pdf")),
      output_dir = !!out.dir)},
  trace = TRUE
) %>%
  filter(!(step %in% c("raw", "combined")) | SeqRunID == 'pb_500')
tictoc::toc()

if (!interactive()) saveRDS(plan, plan_file)

remove(meta1)
remove(meta2)
remove(meta3)
remove(meta4)
remove(datasets)
remove(regions)

cat("\nCalculating outdated targets...\n")
tictoc::tic()
dconfig <- drake_config(plan, jobs_preprocess = local_cpus())
od <- outdated(dconfig)
tictoc::toc()

if (interactive()) {
  # plansteps <- unique(na.omit(plan[["step"]]))
  vis_drake_graph(dconfig,
                  # group = "step", clusters = plansteps,
                  targets_only = TRUE)
}

if (interactive()) {
  dconfig <- drake_config(plan)
  vis_drake_graph(dconfig)
}
remove(dconfig)
remove(snakemake)

save(list = ls(),
     file = drakedata_file)

drake_plan_file_io <- function(plan) {
  if ("file_in" %in% names(plan)) {
    plan$file_in <-
      pmap_chr(plan,
               function(file_in, ...)
                 deparse(eval(file_in, envir = list(...)))
      )
    plan$command <-
      glue_data(plan,
                "{{",
                "file_in({file_in})",
                "{command}",
                "}}",
                .sep = "\n"
      )
  }
  if ("file_out" %in% names(plan)) {
    plan$file_out <-
      pmap_chr(plan,
               function(file_out, ...)
                 deparse(eval(file_out, envir = list(...)))
      )
    plan$command <-
      glue_data(plan,
                "{{",
                "file_out({file_out})",
                "{command}",
                "}}",
                .sep = "\n"
      )
  }
  return(plan)
}

