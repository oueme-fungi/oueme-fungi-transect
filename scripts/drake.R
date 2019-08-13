# define or input file names and parameters
if (interactive()) {
  cat("Creating drake plan in interactive session...\n")
  library(here)
  base_dir <- here()
  r_dir <- here("scripts")
  data_dir <- here("data")
  lab_dir <- here("config")
  seq_dir <- here("sequences")
  trim_dir <- file.path(seq_dir, "trim")
  region_dir <- file.path(seq_dir, "regions")
  filter_dir <- file.path(seq_dir, "filter")
  cluster_dir <- file.path(data_dir, "clusters")
  pasta_dir <- file.path(data_dir, "pasta")
  plan_dir <- file.path(data_dir, "plan")
  ref_dir <- here("reference")
  rmd_dir <- here("writing")
  out_dir <- here("output")
  in_files <- list.files(path = trim_dir,
                         pattern = "\\.trim\\.fastq\\.gz$",
                         full.names = TRUE, recursive = TRUE)
  dataset_file <- file.path(lab_dir, "datasets.csv")
  regions_file <- file.path(lab_dir, "regions.csv")
  platemap_file <- file.path(lab_dir, "Brendan_soil2.xlsx")
  platemap_sheet <- "Concentration samples"
  rdp_file <- file.path(ref_dir, "rdp.LSU.fasta.gz")
  unite_file <- file.path(ref_dir, "unite.fasta.gz")
  unite_patch_file <- file.path(ref_dir, "unite_patch.csv")
  rdp_patch_file <- file.path(ref_dir, "rdp_patch.csv")
  plan_file <- file.path(plan_dir, "plan.rds")
  itsx_meta_file <- file.path(plan_dir, "itsx_meta.rds")
  predada_meta_file <- file.path(plan_dir, "predada_meta.rds")
  dada_meta_file <- file.path(plan_dir, "dada_meta.rds")
  taxonomy_meta_file <- file.path(plan_dir, "taxonomy_meta.rds")
  taxonomy_meta_csv_file <- file.path(plan_dir, "taxonomy_meta.csv")
  tid_file <- file.path(plan_dir, "tids.txt")
  drakedata_file <- file.path(plan_dir, "drake.Rdata")
  longtree_file <- file.path(pasta_dir, "pasta_raxml.tree")
  bigsplit <- 80L
} else if (exists("snakemake")) {
  cat("Creating drake plan in snakemake session...\n")
  snakemake@source(".Rprofile", echo = FALSE)
  r_dir <- snakemake@config$rdir
  dataset_file <- snakemake@input$dataset
  regions_file <- snakemake@input$regions
  platemap_file <- snakemake@input$platemap
  platemap_sheet <- snakemake@config$platemap_sheet
  target <- sub("^.", "", snakemake@output)
  trim_dir <- snakemake@config$trimdir
  region_dir <- snakemake@config$regiondir
  filter_dir <- snakemake@config$filterdir
  ref_dir <- snakemake@config$ref_root
  rdp_file <- snakemake@config$rdp_file
  rdp_patch_file <- snakemake@config$rdp_patch_file
  unite_file <- snakemake@config$unite_file
  unite_patch_file <- snakemake@config$unite_patch_file
  rmd_dir <- snakemake@config$rmddir
  out_dir <- snakemake@config$outdir
  cluster_dir <- snakemake@config$clusterdir
  pasta_dir <- snakemake@config$pastadir
  plan_dir <- snakemake@config$plandir
  plan_file <- snakemake@output$plan
  itsx_meta_file <- snakemake@output$itsx_meta
  predada_meta_file <- snakemake@output$predada_meta
  dada_meta_file <- snakemake@output$dada_meta
  taxonomy_meta_file <- snakemake@output$taxonomy_meta
  drakedata_file <- snakemake@output$drakedata
  taxonomy_meta_csv_file <- snakemake@output$taxonomy_meta_csv
  tids_file <- snakemake@output$tids
  prereqs <- unlist(snakemake@input)
  in_files <- stringr::str_subset(prereqs,
                         pattern = "\\.trim\\.fastq\\.gz$")
  bigsplit <- snakemake@config$bigsplit
  plan_file <- snakemake@output[["plan"]]
  itsx_meta_file <- snakemake@output[["itsx_meta"]]
  predada_meta_file <- snakemake@output[["predada_meta"]]
  dada_meta_file <- snakemake@output[["dada_meta"]]
  taxonomy_meta_file <- snakemake@output[["taxonomy_meta"]]
  pid_file <- snakemake@output[["pids"]]
  tid_file <- snakemake@output[["tids"]]
  drakedata_file <- snakemake@output[["drakedata"]]
  longtree_file <- snakemake@output[["pasta"]][["tree"]]
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
source(file.path(r_dir, "parallel_helpers.R"))
source(file.path(r_dir, "extract_regions.R"))
source(file.path(r_dir, "dada.R"))
source(file.path(r_dir, "taxonomy.R"))
source(file.path(r_dir, "plate_check.R"))
source(file.path(r_dir, "map_LSU.R"))
source(file.path(r_dir, "mantel.R"))

# load the dataset and region definitions
cat("Loading datasets file.\n")
datasets <- read_csv(dataset_file, col_types = "cccicccccic")
cat("Loading regions file.\n")
regions <- read_csv(regions_file, col_types = "cccciiiiic") %>%
  mutate_at("range", replace_na, "") %>%
  mutate(min_length_post = coalesce(min_length_post, min_length),
         max_length_post = coalesce(max_length_post, max_length))

# "meta" dataframes containing definitions for each instance of mapped targets

#### itsx_meta ####
# itsx_meta has one row per demultiplexed, primer-trimmed fastq.gz file
# used in targets derep1 and positions
# The "index" is file_ID
cat("Making itsx_meta.\n")
itsx_meta <- datasets %>%
  mutate(primer_ID = paste0(str_replace(forward, "_tag.*$", ""),
                            str_replace(reverse, "_tag.*$", "")),
         plate = map(runs, ~formatC(seq(.), width = 3, flag = "0"))) %T>%
  {cat("finished first mutate...\n")} %>%
  unnest(plate) %>%
  mutate(direction = if_else(tech == "PacBio",
                             list(c("f", "r")),
                             list(""))) %T>%
  {cat("finished second mutate...\n")} %>%
  unnest(direction) %>%
  mutate(well = list(tidyr::crossing(Row = LETTERS[1:8], Col = 1:12) %>%
                       transmute(well = paste0(Row, Col)))) %T>%
  {cat("finished third mutate...\n")} %>%
  unnest(well) %>%
  mutate(trim_file = glue("{seq_run}_{plate}/{seq_run}_{plate}-{well}{direction}.trim.fastq.gz"),
         seq_run_ID = seq_run,
         
         file_ID = glue("{seq_run}{plate}{well}{direction}") %>%
           str_replace_all("[:punct:]", ""),
         positions = make.names(glue("positions_{file_ID}")),
         derep = make.names(glue("derep1_{file_ID}")),
         join_derep = make.names(glue("join_derep_{primer_ID}")),
         join_derep_map = make.names(glue("join_derep_map_{primer_ID}")),
         itsxtrim = make.names(glue("itsxtrim_{primer_ID}"))) %T>%
  {cat("finished fourth mutate...\n")} %>%
  filter(file.exists(file.path(trim_dir, trim_file))) %>%
  filter(normalizePath(file.path(trim_dir, trim_file)) %in% normalizePath(in_files)) %>%
  filter(file.size(file.path(trim_dir, trim_file)) > 40) %>%
  verify(file.path(trim_dir, trim_file) %in% in_files) %>%
  mutate_at(c("join_derep", "join_derep_map", "itsxtrim", "file_ID", "seq_run_ID", "primer_ID"), syms)

#### predada_meta ####
# predada_meta has one row region per input file
# it is used in targets regions, filter, and derep2
# The index is well_ID + region_ID
# plate_ID is defined because it will be the index for dada_meta
cat("Making predada_meta.\n")
predada_meta <- itsx_meta %>%
  mutate_at("regions", str_split, ",") %>%
  unnest(region = regions, .preserve = ends_with("_ID")) %>%
  left_join(regions, by = "region") %>%
  mutate(well_ID = glue("{seq_run}{plate}{well}{direction}"),
         plate_ID = glue("{seq_run}_{plate}"), 
         region_ID = region,
         range_ID = glue("{region}{range}"),
         region_file = glue("{seq_run}_{plate}-{well}{direction}-{region}.trim.fastq.gz"),
         filter_file = glue("{seq_run}_{plate}-{well}{direction}-{region}{range}.qfilt.fastq.gz")) %>%
  mutate_at(c("positions", "well_ID", "plate_ID", "region_ID", "range_ID"), syms)

if (interactive()) {
  #use a subset
  itsx_meta %<>%
    group_by(seq_run, plate) %>%
    mutate(filesize = file.size(file.path(trim_dir, trim_file))) %>%
    arrange(desc(filesize), .by_group = TRUE) %>%
    dplyr::slice(1:4) %>%
    ungroup()
  
  predada_meta %<>%
    semi_join(itsx_meta, by = "trim_file")
} else {
  saveRDS(itsx_meta, itsx_meta_file)
  saveRDS(predada_meta, predada_meta_file)
}

#### dada_meta ####
# dada_meta has one row per region per plate
# it is used for targets err, dada, dadamap, seq_table, and nochim
# The index is plate_ID + region_ID (plate_ID)
cat("Making dada_meta.\n")
dada_meta <- predada_meta %>%
  select(seq_run, plate, region, range, plate_ID, region_ID, range_ID, seq_run_ID) %>%
  unique() %>%
  arrange(seq_run, plate, region) %>%
  mutate(region_derep = glue("region_derep_{plate_ID}_{range_ID}")) %>%
  left_join(datasets %>% select(-regions), by = "seq_run") %>%
  left_join(regions, by = c("region", "range")) %>%
  mutate_at(c("region_derep"), syms)
if (!interactive()) {
  saveRDS(dada_meta, dada_meta_file)
}

#### region_meta ####
# region_meta has one row per region
# it is used in target big_fasta.
# the index is region_ID
cat("Making region_meta.\n")
region_meta <- regions %>%
  select(region_ID = region) %>%
  mutate(big_seq_table = glue("big_seq_table_{region_ID}"),
         big_fasta_file = glue("{cluster_dir}/{region_ID}.fasta.gz")) %>%
  mutate_at(c("region_ID", "big_seq_table"), syms) %>%
  unique()

#### taxonomy_meta ####
# taxonomy_meta has one row per region and reference DB
# it is used in target taxon,
# and mapped to target guilds_table
# The index is tax_ID (Taxonomy ID)
cat("Making taxonomy_meta.\n")
taxonomy_meta <- dada_meta %>%
  mutate_if(is.list, as.character) %>%
  group_by_at(names(regions)) %>%
  summarize(plate_region_ID = paste(plate_ID, region_ID, sep = "_", collapse = ",")) %>%
  ungroup() %>%
  mutate(reference_file = strsplit(reference, " *, *") %>%
                                map(unlist)) %>%
  unnest(reference_file, .preserve = region) %>%
  mutate(reference = str_split_fixed(reference_file, "\\.", 2)[,1],
         reference_file = glue("{ref_dir}/{reference_file}.vsearch.fasta.gz"),
         tax_ID = glue("{region}_{reference}"),
         big_seq_table = glue("big_seq_table_{region}")) %>%
  arrange(tax_ID) %>%
  mutate_at(c("tax_ID", "big_seq_table"), syms)
if (!interactive()) {
  saveRDS(taxonomy_meta, taxonomy_meta_file)
  mutate_if(taxonomy_meta, is.list, as.character) %>%
  readr::write_csv(taxonomy_meta_csv_file)
  readr::write_lines(taxonomy_meta$tax_ID, tids_file)
}

#### ref_meta ####
# ref_meta has one row per (region specific) reference database.
# it is used for target dbprep
# and mapped to target classifier
cat("Making ref_meta.\n")
ref_meta <- select(taxonomy_meta, reference, reference_file) %>%
  unique() %>%
  mutate(ref_ID = rlang::syms(str_replace(reference_file, "\\.", "_")),
         seq_file = glue("{ref_dir}/{reference_file}.fasta.gz"),
         patch_file = glue("{ref_dir}/{reference}_patch.csv"),
         tax_file = ifelse(reference == "silva",
                           glue("{ref_dir}/{reference}_tax.txt"),
                           NA_character_),
         problem_root = glue("{ref_dir}/{reference_file}_problems"))

#### drake plan ####

cat("\nbuilding plan...\n")
tictoc::tic()
plan <- drake_plan(
  
  datasets = readr::read_csv(file_in(!!dataset_file)),
  
  # derep1 ----
  # dereplicate input sequences before running ITSx
  # this is done independently in parallel
  derep1 = target(
    dada2::derepFastq(file_in(!!file.path(trim_dir, trim_file)),
                      qualityType = "FastqQuality",
                      verbose = TRUE,
                      n = 1e4) %>%
      # add the names of the sequences
      tzara::add_derep_names(file_in(!!file.path(trim_dir, trim_file))) %>%
      # we will not run dada on these, so we don't need the quality
      # information.  Removing it drastically reduces the size in memory.
      magrittr::inset2("quals", NULL),
    transform = map(.data = !!select(itsx_meta, trim_file, file_ID, primer_ID),
                    .tag_in = step,
                    .id = file_ID)),
  
  # join_derep----
  # join the results from dereplicating into a list of unique sequences
  # and a map from the original sequences to the uniques
  join_derep = target(
    tzara::combine_derep(drake_combine(derep1)),
    transform = combine(derep1,
                        .tag_in = step,
                        .by = primer_ID)),
  
  # Get the two parts separately
  join_derep_fasta = target(
    join_derep$fasta,
    transform = map(join_derep, .id = primer_ID)),
  join_derep_map = target(
    join_derep$map,
    transform = map(join_derep, .id = primer_ID)),
  
  # split_fasta----
  # break the list of unique files into equal sized chunks for ITSx
  split_fasta = target(
    split(join_derep$fasta, seq_along(join_derep$fasta) %% !!bigsplit + 1),
    transform = map(join_derep,
                    .tag_in = step)),
  
  # itsx_shard----
  # find the location of LSU, 5.8S, and SSU (and thus ITS1 and ITS2)
  # in the sequences from each chunk
  itsx_shard = target({
    library(Biostrings)
    library(dplyr)
    rITSx::itsx(in.file = split_fasta[[shard]],
                nhmmer = TRUE,
                read_function = Biostrings::readDNAStringSet,
                fasta = FALSE, summary = FALSE, graphical = FALSE,
                save_regions = "none", not_found = FALSE,
                complement = FALSE, cpu = ignore(itsx_cpus))[["positions"]]
    },
    transform = cross(split_fasta,
                      shard = !!(1:bigsplit),
                      .tag_in = step,
                      .id = primer_ID)),
  
  # itsxtrim ----
  # combine the chunks into a master positions list.
  itsxtrim = target(
    dplyr::bind_rows(itsx_shard),
    transform = combine(itsx_shard,
                        .tag_in = step,
                        .by = primer_ID)),
  
  # positions----
  # find the entries in the positions list which correspond to the sequences
  # in each file
  positions = target(
    dplyr::filter(join_derep_map, name == derep) %>% #here
      dplyr::left_join(itsxtrim, by = c("map" = "seq")) %>%
      dplyr::rename(seq = seq.id),
    transform = map(.data = !!select(itsx_meta, derep, join_derep_map, itsxtrim, trim_file,
                                     file_ID),
                    .tag_in = step,
                    .id = file_ID)),
                  
  # regions----
  # cut the required regions out of each sequence
  regions = target(
    tzara::extract_region(seq = file_in(!!file.path(trim_dir, trim_file)),
                          outfile = file_out(!!file.path(region_dir, region_file)),
                          region = region_start,
                          region2 = region_end,
                          positions = positions),
    transform = map(.data = !!unique(select(predada_meta, seq_run, trim_file,
                                     region_file, region, region_start, region_end,
                                     positions,
                                     well_ID, plate_ID, region_ID, seq_run_ID)),
                    .tag_in = step,
                    .id = c(well_ID, region_ID))),
  
  # filter----
  # quality filter the sequences
  filter = target(
    {
      outfile <- file_out(!!file.path(filter_dir, filter_file))
      unlink(outfile)
      out <- dada2::filterAndTrim(
        fwd = file_in(!!file.path(region_dir, region_file)),
        filt = outfile,
        maxLen = max_length,
        minLen = min_length,
        maxEE = max_ee,
        qualityType = "FastqQuality")
      if (!file.exists(outfile)) {
        ShortRead::writeFastq(ShortRead::ShortReadQ(),
                              outfile)
      }
      return(out)
    },
    transform = map(.data = !!select(predada_meta, region_file, filter_file,
                                     max_length, min_length, max_ee, well_ID,
                                     region_ID, range_ID),
                    .tag_in = step,
                    .id = c(well_ID, range_ID))),
  
  # derep2----
  # Do a second round of dereplication on the filtered regions.
  derep2 = target({
    infile <- file_in(!!file.path(filter_dir, filter_file))
    if (file.size(infile) > 50) {
      out <- dada2::derepFastq(fls = infile, n = 1e4,
                               qualityType = "FastqQuality",
                               verbose = TRUE) %>%
        tzara::add_derep_names(infile)
    } else {
      out = NULL
    }
    result <- list()
    result[[filter_file]] <- out
    result
  },
  transform = map(.data = !!predada_meta,
                  .tag_in = step,
                  .id = c(well_ID, range_ID))),
  
  # region_derep ----
  # put together all the dereplicated files for each region for each plate.
  region_derep = target(
    purrr::compact(c(derep2)),
    transform = combine(derep2, .by = c(plate_ID, range_ID), .tag_in = step),
  ),
  
  # err----
  # Calibrate dada error models (on different size ranges)
  err = target({
    err.fun <- ifelse(tech == "PacBio",
                      dada2::PacBioErrfun,
                      dada2::loessErrfun)
    dada2::learnErrors(
      fls = region_derep,
      nbases = 1e9,
      multithread = ignore(dada_cpus),
      randomize = TRUE,
      errorEstimationFunction = err.fun, 
      HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
      BAND_SIZE = band_size,
      pool = eval(rlang::parse_expr(pool)),
      verbose = TRUE,
      qualityType = "FastqQuality")},
    transform = map(.data = !!dada_meta,
                    .tag_in = step, .id = c(plate_ID, range_ID))),
  
  # dada----
  # Run dada denoising algorithm (on different size ranges)
  dada = target(
    dada2::dada(
      derep = region_derep, err = err,
      multithread = ignore(dada_cpus), 
      HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
      BAND_SIZE = band_size,
      pool = eval(rlang::parse_expr(pool))),
    transform = map(err, .data = !!dada_meta,
                    .tag_in = step,
                    .id = c(plate_ID, range_ID))),
  
  # dadamap----
  # Make maps from individual sequences in source files to dada ASVs
  # (for each size range, if multiple)
  dada_map = target(
    tzara::dadamap(region_derep, dada),
    transform = map(dada, .data = !!dada_meta,
                    .tag_in = step, .id = c(plate_ID, range_ID))),
  
  # seq_table----
  # Make a sample x ASV abundance matrix
  # (for each size range, if multiple)
  seq_table = target(
    dada2::makeSequenceTable(dada) %>%
      magrittr::extract(,nchar(colnames(.)) >= min_length_post &
                          nchar(colnames(.)) <= max_length_post),
    transform = map(dada, .data = !!dada_meta, .tag_in = step, .id = c(plate_ID, range_ID))),
  
  # nochim----
  # Remove likely bimeras
  nochim = target(
    dada2::removeBimeraDenovo(
      list(seq_table) %>%
        purrr::map(tibble::as_tibble, rownames = "file") %>%
        purrr::reduce(dplyr::full_join, by = "file") %>%
        tibble::column_to_rownames("file") %>%
        as.matrix,
      method = "consensus",
      multithread = ignore(dada_cpus)),
    transform = combine(seq_table, .by = c(plate_ID, region_ID))),
  
  # big_seq_table ----
  # Join all the sequence tables for each region
  big_seq_table = target(
    dada2::mergeSequenceTables(tables = list(nochim)),
    transform = combine(nochim, .tag_in = step, .by = region_ID)),
  
  # big_fasta ----
  # write the ITS2 big_seq_table as a fasta file so that it can be clustered by
  # VSEARCH.
  big_fasta = target(
    write_big_fasta(big_seq_table,
                    file_out(!!big_fasta_file)),
    transform = map(.data = !!region_meta, .id = region_ID)),
  
  # dbprep ----
  # import the RDP, Silva, and UNITE databases, and format them for use with
  # DECIPHER.
  dbprep = target(
    formatdb(db = reference,
             seq_file = file_in(seq_file),
             patch_file = file_in(patch_file),
             tax_file = file_in(tax_file),
             maxGroupSize = 20),
    
    transform = map(.data = !!ref_meta,
                    .id = ref_ID)),
  # classifier----
  # Train a DECIPHER classifier for each database.
  classifier = target(
    refine_classifier(dbprep$train,
                      dbprep$taxonomy,
                      dbprep$rank,
                      out_root = problem_root),
    transform = map(dbprep, .id = ref_ID)
  ),
  
  # taxon ----
  # Assign taxonomy to each ASV
  taxon = target(
    colnames(big_seq_table) %>%
      sintax(db = file_in(!!file.path(ref_dir, paste0(reference_file, ".vsearch.fasta.gz"))),
             sintax_cutoff = 0.9,
             multithread = ignore(dada_cpus)) %>%
      sintax_format(),
    transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID)),
  
  # funguild_db ----
  # Download the FUNGuild database
  funguild_db = FUNGuildR::get_funguild_db(),
  
  # guilds_table ----
  # Assign ecological guilds to the ASVs based on taxonomy.
  guilds_table = target(
    FUNGuildR::funguild_assign(taxon, db = funguild_db),
    transform = map(taxon, tax_ID, .tag_in = step, .id = tax_ID)),
  
  # platemap ----
  # Read the map between plate locations and samples
  platemap = read_platemap(file_in(!!platemap_file), platemap_sheet),
  
  # raw ----
  # Join the raw read info (for long pacbio reads).
  raw = target(
    tzara::summarize_sread(list(regions),
                           name = !!symbols_to_values(filter_file),
                           max_ee = !!eval(parse(text = unique(symbols_to_values(max_ee))))),
    transform = combine(regions, filter_file, max_ee, seq_run, region,
                        .by = c(seq_run_ID, region_ID),
                        .tag_in = step,
                        .id = c(seq_run_ID, region_ID))
  ),
  
  # combined ----
  # Replace raw reads which were mapped to a DADA2 ASV with the ASV,
  # and put them all in one tibble.  We'll only do this for pacbio long reads.
  combined = target(
    tzara::combine_bigmaps(purrr::imap_dfr(c(dada_map), ~mutate(.x, name = .y)),
                           dplyr::bind_rows(raw)),
    transform = combine(dada_map, raw, .by = seq_run_ID, .tag_in = step)
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
    sintax(db = file_in("reference/rdp_train.LSU.fasta.gz"),
           sintax_cutoff = 0.9,
           multithread = ignore(dada_cpus)) %>%
    dplyr::mutate_at("seq", stringr::str_replace_all, "T", "U"),
  
  # cons_tax----
  # Make names which combine the ITS2 and LSU identifications to put on the tree
  # Also make hash names, because PASTA will destroy non-alphanumeric names
  cons_tax  =
    conseq_filt %>%
    dplyr::left_join(dplyr::select(taxon_ITS2_unite, ITS2 = seq, name),
                     by = "ITS2") %>%
    dplyr::left_join(dplyr::select(taxon_LSUcons_rdp, LSU = seq, name),
                     suffix = c("_ITS2", "_LSU"),
                     by = "LSU") %>%
    tidyr::unite("name", name_ITS2, name_LSU, sep = "/"),
  
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
    Sys.setFileTime(file.path(pasta_dir, ".dopasta"), Sys.time())
  },
  
  # longtree----
  # read the long amplicon tree in
  longtree = {
    lt = ape::read.tree(file_in(longtree_file))
    phyloseq::taxa_names(lt) <-
      tibble::tibble(hash = phyloseq::taxa_names(lt)) %>%
      dplyr::left_join(conseq_filt) %$%
      seqhash(ITS2)
    lt},
  
  physeq = assemble_physeq(platemap, datasets, relabel_seqtable(big_seq_table_ITS2), longtree),
  
  longphyseq = physeq %>%
    phyloseq::prune_samples(samples = phyloseq::sample_data(.)$Dataset == "long-pacbio"
                            & phyloseq::sample_data(.)$sample_type == "Sample"
                            & rowSums(phyloseq::otu_table(.)) > 0),
  
  unif_dist = phyloseq::distance(longphyseq, "unifrac"),
  bc_dist = phyloseq::distance(longphyseq, "bray"),
  spatial_dist = phyloseq::sample_data(longphyseq) %>%
    with(X + 30000 * as.integer(Site) + 100000 * as.integer(Year)) %>%
    dist,
  spatial_dist2 = spatial_dist + ifelse(spatial_dist > 50000, -100000, 100000),
  unif_corr = vegan::mantel.correlog(unif_dist, spatial_dist,
                                      break.pts = 0:13 - 0.5,
                                      cutoff = FALSE),
  bc_corr = vegan::mantel.correlog(bc_dist, spatial_dist,
                                    break.pts = 0:13 - 0.5,
                                    cutoff = FALSE),
  
  
  unif_corr2 = vegan::mantel.correlog(unif_dist, spatial_dist2,
                                       break.pts = 0:13 - 0.5,
                                       cutoff = FALSE),
  bc_corr2 = vegan::mantel.correlog(bc_dist, spatial_dist2,
                                     break.pts = 0:13 - 0.5,
                                     cutoff = FALSE),
  
  # seq_count ----
  # Count the sequences in each fastq file
  seq_count = target(
    tibble::tibble(
      file = file_in(fq_file),
      reads = system(glue::glue("zcat {file} | grep -c '^@'", file = fq_file),
                     intern = TRUE) %>%
        as.integer),
    transform = map(fq_file = !!c(file.path(trim_dir, itsx_meta$trim_file),
                                 file.path(region_dir, predada_meta$region_file),
                                 file.path(filter_dir, predada_meta$filter_file)),
                    .tag_in = step)),
  
  # seq_counts----
  # merge all the sequence counts into one data.frame
  seq_counts = target(
    dplyr::bind_rows(seq_count),
    transform = combine(seq_count)),

  # qstats----
  # calculate quality stats for the different regions
  qstats = target(
    q_stats(file_in(!!file.path(region_dir, region_file))),
    transform = map(.data = !!predada_meta,
                    .tag_in = step,
                    .id = c(well_ID, region_ID))),
  #qstats_join----
  # join all the quality stats into one data.frame
  qstats_join = target(
    dplyr::bind_rows(qstats) %>%
      tidyr::extract(
        col = "file",
        into = c("seq_run", "plate", "well", "direction", "region"),
        regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.trim\\.fastq\\.gz") %>%
      dplyr::left_join(datasets),
    transform = combine(qstats)),
  # qstats_knit----
  # knit a report about the quality stats.
  qstats_knit = {
    if (!dir.exists(!!out_dir)) dir.create(!!out_dir, recursive = TRUE)
    rmarkdown::render(
      knitr_in(!!file.path(rmd_dir, "qual-check.Rmd")),
      output_file = file_out(!!file.path(out_dir, "qual-check.pdf")),
      output_dir = !!out_dir)},
  trace = TRUE
) %>%
  filter(!(step %in% c("raw", "combined")) | seq_run_ID == 'pb_500')
tictoc::toc()

if (!interactive()) saveRDS(plan, plan_file)

remove(itsx_meta)
remove(predada_meta)
remove(dada_meta)
remove(taxonomy_meta)
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

