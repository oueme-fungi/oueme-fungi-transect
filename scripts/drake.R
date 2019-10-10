library(futile.logger)

# define or input file names and parameters
if (interactive()) {
  flog.info("Creating drake plan in interactive session...")
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
  locarna_dir <- file.path(data_dir, "mlocarna")
  plan_dir <- file.path(data_dir, "plan")
  ref_dir <- here("reference")
  rmd_dir <- here("writing")
  out_dir <- here("output")
  in_files <- list.files(path = trim_dir,
                         pattern = "\\.trim\\.fastq\\.gz$",
                         full.names = TRUE, recursive = TRUE)
  in_files <- in_files[!grepl("is_057", in_files)]
  dataset_file <- file.path(lab_dir, "datasets.csv")
  regions_file <- file.path(lab_dir, "regions.csv")
  methods_file <- file.path(lab_dir, "taxonomy_methods.csv")
  platemap_file <- file.path(lab_dir, "Brendan_soil2.xlsx")
  platemap_sheet <- "Concentration samples"
  rdp_file <- file.path(ref_dir, "rdp.LSU.fasta.gz")
  unite_file <- file.path(ref_dir, "unite.fasta.gz")
  unite_patch_file <- file.path(ref_dir, "unite_patch.csv")
  rdp_patch_file <- file.path(ref_dir, "rdp_patch.csv")
  cm_5.8S <- file.path(ref_dir, "RF00002.cm")
  cm_32S <- file.path(ref_dir, "fungi_32S_LR5.cm")
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
  flog.info("Creating drake plan in snakemake session...")
  snakemake@source(".Rprofile", echo = FALSE)
  r_dir <- snakemake@config$rdir
  dataset_file <- snakemake@input$dataset
  regions_file <- snakemake@input$regions
  methods_file <- snakemake@input$methods
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
  cm_5.8S <- snakemake@config$cm_5_8S
  cm_32S <- snakemake@config$cm_32S
  rmd_dir <- snakemake@config$rmddir
  out_dir <- snakemake@config$outdir
  cluster_dir <- snakemake@config$clusterdir
  pasta_dir <- snakemake@config$pastadir
  locarna_dir <- snakemake@config$locarnadir
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
source(file.path(r_dir, "inferrnal.R"))

setup_log("drakeplan")

# load the dataset and region definitions
flog.info("Loading datasets file.")
datasets <- read_csv(dataset_file, col_types = "cccicccccic")
flog.info("Loading regions file.\n")
regions <- read_csv(regions_file, col_types = "cccciiiiic") %>%
  mutate_at("range", replace_na, "") %>%
  mutate(min_length_post = coalesce(min_length_post, min_length),
         max_length_post = coalesce(max_length_post, max_length))
flog.info("Loading taxonomic assignment methods file.")
methods <- read_csv(methods_file, col_types = "ci")

# "meta" dataframes containing definitions for each instance of mapped targets

#### itsx_meta ####
# itsx_meta has one row per demultiplexed, primer-trimmed fastq.gz file
# used in targets derep1 and positions
# The "index" is file_ID
flog.info("Making itsx_meta.")
itsx_meta <- datasets %>%
  mutate(primer_ID = paste0(str_replace(forward, "_tag.*$", ""),
                            str_replace(reverse, "_tag.*$", "")),
         plate = map(runs, ~formatC(seq(.), width = 3, flag = "0"))) %T>%
         {flog.info("finished first mutate...")} %>%
  unnest(plate) %>%
  mutate(direction = if_else(tech == "PacBio",
                             list(c("f", "r")),
                             list(""))) %T>%
                             {flog.info("finished second mutate...")} %>%
  unnest(direction) %>%
  mutate(well = list(tidyr::crossing(Row = LETTERS[1:8], Col = 1:12) %>%
                       transmute(well = paste0(Row, Col)))) %T>%
                       {flog.info("finished third mutate...")} %>%
  unnest(well) %>%
  mutate(trim_file = glue("{seq_run}_{plate}/{seq_run}_{plate}-{well}{direction}.trim.fastq.gz"),
         file_ID = glue("{seq_run}{plate}{well}{direction}") %>%
           str_replace_all("[:punct:]", ""),
         positions = make.names(glue("positions_{file_ID}")),
         derep = make.names(glue("derep1_{file_ID}")),
         join_derep = make.names(glue("join_derep_{primer_ID}")),
         join_derep_map = make.names(glue("join_derep_map_{primer_ID}")),
         itsxtrim = make.names(glue("itsxtrim_{primer_ID}")),
         lsux_pos = make.names(glue("lsux_pos_{primer_ID}")),
         search5_8S = make.names(glue("search5_8S_{primer_ID}"))) %T>%
         {flog.info("finished fourth mutate...")} %>%
  filter(file.exists(file.path(trim_dir, trim_file))) %>%
  filter(normalizePath(file.path(trim_dir, trim_file)) %in% normalizePath(in_files)) %>%
  filter(file.size(file.path(trim_dir, trim_file)) > 40) %>%
  verify(file.path(trim_dir, trim_file) %in% in_files) %>%
  mutate_at(c("join_derep", "join_derep_map", "itsxtrim", "lsux_pos", "search5_8S"), syms)

#### predada_meta ####
# predada_meta has one row per region per well
# it is used in targets regionx, filter, and derep2
# The index is well_ID + region
# plate_ID is defined because it will be the index for dada_meta
flog.info("Making predada_meta.")
predada_meta <- itsx_meta %>%
  group_by_at(c("seq_run", "regions", "plate", "well")) %>%
  summarize(trim_file = list(c(trim_file)),
            positions = list(c(positions))) %>%
  ungroup() %>%
  mutate_at("regions", str_split, ",") %>%
  unnest(region = regions, .preserve = c("trim_file", "positions")) %>%
  left_join(regions, by = "region") %>%
  mutate(well_ID = glue("{seq_run}{plate}{well}"),
         plate_ID = glue("{seq_run}_{plate}"),
         range_ID = glue("{region}{range}"),
         region_file = glue("{seq_run}_{plate}-{well}-{region}.trim.fastq.gz"),
         filter_file = glue("{seq_run}_{plate}-{well}-{region}{range}.qfilt.fastq.gz")) %>%
  mutate_at("positions", lapply, syms)

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
# dada_meta has one row per region per sequencing run
# it is used for targets err, dada, dadamap, seq_table, and nochim
# The index is seq_run + region
flog.info("Making dada_meta.")
dada_meta <- predada_meta %>%
  select(seq_run, region, range, range_ID) %>%
  unique() %>%
  arrange(seq_run, region) %>%
  mutate(region_derep = glue("region_derep_{seq_run}_{range_ID}"),
         error_model = glue("err_{seq_run}_5_8S")) %>%
  left_join(datasets %>% select(-regions), by = "seq_run") %>%
  left_join(regions, by = c("region", "range")) %>%
  mutate_at(c("region_derep", "error_model"), syms)
if (!interactive()) {
  saveRDS(dada_meta, dada_meta_file)
}

#### region_meta ####
# region_meta has one row per region
# it is used in target big_fasta.
# the index is region
flog.info("Making region_meta.")
region_meta <- regions %>%
  select(region) %>%
  mutate(big_seq_table = glue("big_seq_table_{region}"),
         big_fasta_file = glue("{cluster_dir}/{region}.fasta.gz")) %>%
  mutate_at("big_seq_table", syms) %>%
  unique()

#### taxonomy_meta ####
# taxonomy_meta has one row per region and reference DB
# it is used in target taxon,
# and mapped to target guilds_table
# The index is tax_ID (Taxonomy ID)
flog.info("Making taxonomy_meta.")
taxonomy_meta <- dada_meta %>%
  mutate_if(is.list, as.character) %>%
  select_at(names(regions)) %>%
  separate_rows(reference, sep = " *, *") %>%
  separate(reference, c("reference", "refregion"), sep = "\\.") %>%
  filter(!is.na(reference)) %>%
  crossing(methods) %>%
  mutate(methodfile = ifelse(method == "idtaxa", "sintax", method),
         reference_file = glue("{ref_dir}/{reference}.{refregion}.{methodfile}.fasta.gz"),
         refdb = glue("refdb_{method}_{reference}_{refregion}"),
         tax_ID = glue("{region}_{reference}_{refregion}_{method}"),
         big_seq_table = glue("big_seq_table_{region}")) %>%
  arrange(tax_ID) %>%
  mutate_at(c("big_seq_table", "refdb"), syms)
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
flog.info("Making ref_meta.")
ref_meta <- select(taxonomy_meta, "reference", "refregion", "method", "reference_file") %>%
  unique() %>%
  mutate(ref_ID = glue("{reference}_{refregion}"))

#### drake plan ####

flog.info("\nbuilding plan...")
tictoc::tic()
plan <- drake_plan(
  
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
    transform = map(.data = !!select(itsx_meta, "trim_file", "file_ID", "primer_ID"),
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
    transform = map(join_derep, .id = primer_ID),
    format = "fst"),
  
  # split_fasta----
  # break the list of unique files into equal sized chunks for ITSx
  split_fasta = target(
    split(join_derep_fasta, seq_along(join_derep_fasta) %% !!bigsplit + 1),
    transform = map(join_derep_fasta,
                    .id = primer_ID,
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
                    .id = primer_ID),
  format = "fst"),
  
  # itsxtrim ----
  # combine the chunks into a master positions list.
  itsxtrim = target(
    dplyr::bind_rows(itsx_shard),
    transform = combine(itsx_shard,
                        .tag_in = step,
                        .by = primer_ID),
    format = "fst"),
  
  lsux_shard = target(
    LSUx(seq = split_fasta[[shard]],
         cm_5.8S = file_in(!!cm_5.8S),
         cm_32S = file_in(!!cm_32S),
         glocal = TRUE,
         ITS1 = TRUE,
         cpu = ignore(itsx_cpus)) %>%
      dplyr::mutate_at("seq_name", as.integer),
    transform = cross(split_fasta,
                      shard = !!(1:bigsplit),
                      .id = primer_ID),
    format = "fst"
  ),
  
  lsux_pos = target(
    dplyr::bind_rows(lsux_shard) %>%
      gather_regions(),
    transform = combine(lsux_shard,
                        .by = primer_ID),
    format = "fst"
  ),
  
  # positions----
  # find the entries in the positions list which correspond to the sequences
  # in each file
  positions = target(
    lsux_pos %>%
      dplyr::left_join(dplyr::filter(join_derep_map, name == derep), .,
                       by = c("map" = "seq_name")) %>%
      dplyr::rename(seq = seq.id),
    transform = map(.data = !!select(itsx_meta, "derep", "join_derep_map",
                                     "lsux_pos", "trim_file", "file_ID"),
                    .tag_in = step,
                    .id = file_ID),
    format = "fst"),
  
  # regionx----
  # cut the required regions out of each sequence
  regionx = target(
    tzara::extract_region(seq = file_in(!!file.path(trim_dir, trim_file)),
                          outfile = file_out(!!file.path(region_dir, region_file)),
                          region = region_start,
                          region2 = region_end,
                          positions = positions),
    transform = map(.data = !!predada_meta,
                    .tag_in = step,
                    .id = c(well_ID, region))),
  
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
    transform = map(.data = !!predada_meta,
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
    transform = combine(derep2, .by = c(seq_run, range_ID), .tag_in = step),
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
    transform = map(.data = !!filter(dada_meta, range_ID == "5_8S"),
                    .tag_in = step, .id = c(seq_run, range_ID))),
  
  # dada----
  # Run dada denoising algorithm (on different size ranges)
  dada = target(
    dada2::dada(
      derep = region_derep, err = error_model,
      multithread = ignore(dada_cpus), 
      HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
      BAND_SIZE = band_size,
      pool = eval(rlang::parse_expr(pool))),
    transform = map(.data = !!dada_meta,
                    .tag_in = step,
                    .id = c(seq_run, range_ID))),
  
  # dadamap----
  # Make maps from individual sequences in source files to dada ASVs
  # (for each size range, if multiple)
  dada_map = target(
    tzara::dadamap(region_derep, dada),
    transform = map(dada,
                    .tag_in = step, .id = c(seq_run, range_ID))),
  
  # seq_table----
  # Make a sample x ASV abundance matrix
  # (for each size range, if multiple)
  seq_table = target(
    dada2::makeSequenceTable(dada) %>%
      magrittr::extract(,nchar(colnames(.)) >= min_length_post &
                          nchar(colnames(.)) <= max_length_post),
    transform = map(dada, .tag_in = step, .id = c(seq_run, range_ID))),
  
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
    transform = combine(seq_table, .by = c(seq_run, region))),
  
  # big_seq_table ----
  # Join all the sequence tables for each region
  big_seq_table = target(
    dada2::mergeSequenceTables(tables = list(nochim)),
    transform = combine(nochim, .tag_in = step, .by = region)),
  
  # big_fasta ----
  # write the ITS2 big_seq_table as a fasta file so that it can be clustered by
  # VSEARCH.
  big_fasta = target(
    write_big_fasta(big_seq_table,
                    file_out(!!big_fasta_file)),
    transform = map(.data = !!region_meta, .id = region)),
  
  # raw ----
  # Join the raw read info (for long pacbio reads).
  raw = target(
    tzara::summarize_sread(list(regionx),
                           name = !!symbols_to_values(filter_file),
                           max_ee = 1.5 * !!eval(parse(text = unique(symbols_to_values(max_ee))))),
    transform = combine(regionx, filter_file, max_ee, seq_run, region,
                        .by = c(seq_run, region),
                        .tag_in = step,
                        .id = c(seq_run, region)),
    format = "fst"
  ),
  
  # combined ----
  # Replace raw reads which were mapped to a DADA2 ASV with the ASV,
  # and put them all in one tibble.  We'll only do this for pacbio long reads.
  combined = target(
    combine_bigmaps(dplyr::bind_rows(dada_map),
                    dplyr::bind_rows(raw)),
    transform = combine(dada_map, raw, .by = seq_run, .tag_in = step),
    format = "fst"
  ),
  
  # reconstructed ----
  # 
  reconstructed = target(
    tzara::reconstruct(seqtabs = list(dada_map),
                       rawtabs = list(raw),
                       regions = !!symbols_to_values(region),
                       order = c("ITS1", "5_8S", "ITS2", "LSU1", "D1", "LSU2",
                                 "D2", "LSU3", "D3", "LSU4"),
                       output = "long",
                       read_column = "seq.id",
                       asv_column = "dada.seq",
                       raw_column = "seq",
                       sample_column = "name",
                       sample_regex = paste0("(pb|is)_\\d{3}_\\d{3}-[A-H]1?\\d"),
                       ncpus = ignore(dada_cores)) %>%
      dplyr::mutate(`32S` = stringr::str_c(`5_8S`, ITS2, LSU1, D1, LSU2, D2, LSU3,
                                           D3, LSU4),
                    ITS = stringr::str_c(ITS1, `5_8S`, ITS2),
                    LSU = stringr::str_c(LSU1, D1, LSU2, D2, LSU3, D3, LSU4),
                    hash = tzara::seqhash(long)),
    transform = combine(dada_map, raw, region, .by = seq_run, .tag_in = step),
    format = "fst"
  ),
  
  # conseq ----
  # Calculate consensus sequences for full ITS and LSU for each ITS2 ASV with
  # with at least 3 sequences.
  preconseq = target(
    combined_pb_500 %>%
      tidyr::extract(name, c("seq_run", "plate", "well", "direction", "region"), "([pi][sb]_\\d{3})_(\\d{3})-([A-H]\\d{1,2})([fr])-(.+)\\.qfilt\\.fastq\\.gz") %>%
      tidyr::spread(key = "region", value = "seq") %>%
      dplyr::group_by(ITS2) %>%
      dplyr::filter(!is.na(ITS2), dplyr::n() >= 3)),
  
  conseq = target(
    preconseq %>%
      dplyr::summarize(
        nreads = dplyr::n(),
        !! region := as.character(
          tzara::cluster_consensus(.data[[region]], seq.id, ignore(dada_cpus)))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(stringr::str_count(.[[region]], "[MRWSYKVHDBN]") < 3,
                    !is.na(.[[region]]),
                    nchar(.[[region]]) >= min_length),
    transform = map(.data = !!filter(regions, region %in% c("long", "ITS", "ITS1", "LSU", "LSU1", "D1", "LSU2", "D2", "LSU3", "D3", "LSU4", "32S", "5_8S")), .id = region),
    format = "fst"),
  
  # allseqs ----
  # make a data frame of all consensus sequences and ASVs; each row is an ITS2 ASV,
  # or if 
  allseqs = target(
    make_allseq_table(list(conseq),
                      drake_combine(big_seq_table)),
    transform = combine(conseq, big_seq_table),
    format = "fst"
  ),
  
  # dbprep ----
  # import the RDP, Silva, and UNITE databases, and format them for use with
  # DECIPHER.
  refdb_dada2 = target(
    file_in(reference_file),
    transform = map(.data = !!filter(ref_meta, method == "dada2"),
                    .tag_out = refdb,
                    .id = ref_ID)),
  
  
  refdb_sintax = target(
    file_in(reference_file),
    transform = map(.data = !!filter(ref_meta, method == "sintax"),
                    .tag_out = refdb,
                    .id = ref_ID)),
  
  refdb_idtaxa = target(
    fit_idtaxa(file_in(reference_file)),
    transform = map(.data = !!filter(ref_meta, method == "idtaxa"),
                    .tag_out = refdb,
                    .id = ref_ID)),
  
  # taxon ----
  # Assign taxonomy to each ASV
  taxon = target(
    dplyr::select(reconstructed, region, "hash") %>%
      dplyr::filter(complete.cases(.)) %>%
      {set_names(.[[region]], .$hash)} %>%
      taxonomy(reference = refdb,
               method = method,
               multithread = ignore(dada_cpus)) %>%
      taxtable(),
    transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID),
    format = "fst"),
  
  # funguild_db ----
  # Download the FUNGuild database
  funguild_db = FUNGuildR::get_funguild_db(),
  
  # guilds_table ----
  # Assign ecological guilds to the ASVs based on taxonomy.
  guilds_table = target(
    FUNGuildR::funguild_assign(taxon, db = funguild_db),
    transform = map(taxon, tax_ID, .tag_in = step, .id = tax_ID),
    format = "fst"),
  
  # platemap ----
  # Read the map between plate locations and samples
  platemap = target(
    read_platemap(file_in(!!platemap_file), platemap_sheet),
    format = "fst"),
  
  
  # taxon_LSUcons_rdp----
  # Assign taxonomy to the consensus LSU sequences.
  taxon_LSUcons_rdp = target(
    allseqs$LSU %>%
      unique() %>%
      stringr::str_replace_all("U", "T") %>%
      sintax(db = file_in("reference/rdp_train.LSU.fasta.gz"),
             sintax_cutoff = 0.9,
             multithread = ignore(dada_cpus)) %>%
      dplyr::mutate_at("seq", stringr::str_replace_all, "T", "U"),
    format = "fst"),
  
  # cons_tax----
  # Make names which combine the ITS2 and LSU identifications to put on the tree
  # Also make hash names, because PASTA will destroy non-alphanumeric names
  cons_tax  = target(
    allseqs %>%
      dplyr::left_join(dplyr::select(taxon_ITS2_unite, ITS2 = seq, name),
                       by = "ITS2") %>%
      dplyr::left_join(dplyr::select(taxon_LSUcons_rdp, LSU = seq, name),
                       suffix = c("_ITS2", "_LSU"),
                       by = "LSU") %>%
      tidyr::unite("name", name_ITS2, name_LSU, sep = "/"),
    format = "fst"),
  
  # long_consensus----
  # get the long amplicon consensus and convert to RNAStringSet.
  # Use the sequence hashes as the name, so that names will be robust in PASTA
  reconst_32S = reconstructed %>%
    dplyr::select("32S", ITS1, hash) %>%
    dplyr::filter(complete.cases(.)) %$%
    rlang::set_names(`32S`, hash) %>%
    chartr(old = "T", new = "U") %>% 
    Biostrings::RNAStringSet(),
  reconst_LSU = reconstructed %>%
    dplyr::select(LSU, ITS, hash) %>%
    dplyr::filter(complete.cases(.)) %$%
    rlang::set_names(LSU, hash) %>%
    chartr(old = "T", new = "U") %>%
    Biostrings::RNAStringSet(),
  reconst_long = reconstructed %>%
    dplyr::select(long, ITS, hash) %>%
    dplyr::filter(complete.cases(.)) %$%
    rlang::set_names(long, hash) %>%
    chartr(old = "T", new = "U") %>%
    Biostrings::RNAStringSet(),
  
  # cmaln_32S
  # Align conserved positions and annotate conserved secondary structure of the
  # 32S consensus using Infernal.
  cmaln_32S =
    cmalign(cmfile = file_in(!!cm_32S),
            seq = reconst_32S,
            glocal = TRUE,
            cpu = ignore(dada_cpus)),
  
  # cmaln_long
  # Concatenate ITS1 (padded with gaps) to the beginning of the Infernal 32S
  # alignment to make a preliminary long amplicon alignment with conserved
  # 5.8S/LSU secondary structure annotations.
  cmaln_long =
    dplyr::inner_join(
      reconstructed %>%
        dplyr::select(hash, ITS1) %>%
        dplyr::filter(complete.cases(.)) %>%
        unique(),
      tibble::tibble(
        hash = cmaln_32S$names,
        aln = as.character(cmaln_32S$alignment@unmasked)
      ),
      by = "hash"
    ) %>%
    dplyr::mutate(
      ITS1 = chartr("T", "U", ITS1),
      ITS1 = stringr::str_pad(ITS1, max(nchar(ITS1)), "right", "-"),
      aln = stringr::str_c(ITS1, aln)
    ) %>% {
      write_clustalw_ss(
        aln = magrittr::set_names(.$aln, .$hash) %>%
          Biostrings::RNAStringSet(),
        sec_str = paste0(
          # The first four bases after the ITS1 primer are paired to other
          # parts of SSU, so should not be paired in the rest of the amplicon.
          # Otherwise, we have no structure information for the ITS1 region (it
          # does not have eukaryote-wide or fungi-wide conserved structure)
          stringr::str_pad(
            "xxxx",
            max(nchar(reconstruct_long_aln$ITS1)),
            "right",
            "-"),
          cmaln_32S$SS_cons
        ),
        # Don't take the alignment in the variable regions as
        # conserved.
        ref = stringr::str_pad(chartr("v", ".", cmaln_32S$RF),
                               max(nchar(.$aln)), "left", "-"),
        seq_names = .$hash,
        file = file_out(cmaln_file_long)
      )
    },
  
  # guidetree_32S ----
  # Take only the conserved, gap-free positions in the 32S alignment
  aln_32S_conserv =
    remove_nonconsensus_nongaps(aln_32S),
  
  # Generate a distance matrix from the conserved 32S alignment.
  dist_32S_conserv =
    DECIPHER::DistanceMatrix(
      type = "dist",
      processors = ignore(dada_cpus)
    ),
  
  # Make a UPGMA guide tree for mlocarna based on the aligned positions in
  # the Infernal alignment for 32S.
  guidetree_32S =
    DECIPHER::IdClusters(
      myDistMatrix = dist_32S_conserv,
      collapse = -1, #break polytomies
      type = "dendrogram",
      processors = ignore(dada_cpus)
    ) %T>%
    DECIPHER::WriteDendrogram(
      file = file_out(guide_tree_file)
    ),
  
  # DECIPHER has a fast progressive alignment algorithm that calculates RNA
  # secondary structure
  aln_LSU =
    DECIPHER::AlignSeqs(consensus_LSU,
                        iterations = 10, refinements = 10,
                        processors = ignore(dada_cpus)),
  
  
  
  # write_lsualn----
  # write the aligned consensus for PASTA.
  # also touch a file to mark that we really do need to re-run PASTA.
  write_aln_LSU = {
    Biostrings::writeXStringSet(aln_LSU,
                                file_out(aln_file_LSU))
    Sys.setFileTime(file.path(pasta_dir, ".dopasta"), Sys.time())
  },

  write_aln_32S = 
    write_clustalw_ss(aln = aln_32S$alignment@unmasked,
                      sec_str = aln_32S$SS_cons,
                      file = file_out(aln_file_32S)),
  
  # longtree----
  # read the long amplicon tree in
  longtree = {
    lt = ape::read.tree(file_in(longtree_file))
    phyloseq::taxa_names(lt) <-
      tibble::tibble(hash = phyloseq::taxa_names(lt)) %>%
      dplyr::left_join(allseqs) %$%
      hash
    lt},
  
  physeq = assemble_physeq(platemap, datasets, relabel_seqtable(big_seq_table_ITS2)),#, longtree),
  
  longphyseq = physeq %>%
    phyloseq::prune_samples(samples = phyloseq::sample_data(.)$dataset == "long-pacbio"
                            & phyloseq::sample_data(.)$sample_type == "Sample"
                            & (phyloseq::sample_data(.)$qual == "X" | phyloseq::sample_data(.)$year == "2015")
                            & rowSums(phyloseq::otu_table(.)) > 0),
  
  shortphyseq = physeq %>%
    phyloseq::prune_samples(samples = phyloseq::sample_data(.)$dataset == "short-pacbio"
                            & phyloseq::sample_data(.)$sample_type == "Sample"
                            & (phyloseq::sample_data(.)$qual == "X" | phyloseq::sample_data(.)$year == "2015")
                            & rowSums(phyloseq::otu_table(.)) > 0),
  
  unif_dist = phyloseq::distance(longphyseq, "unifrac"),
  bc_dist_long = phyloseq::distance(longphyseq, "bray"),
  bc_dist_short = phyloseq::distance(shortphyseq, "bray"),
  spatial_dist_long = phyloseq::sample_data(longphyseq) %>%
    with(x + 30000 * as.integer(site) + 100000 * as.integer(year)) %>%
    dist,
  spatial_dist_long2 = spatial_dist_long + ifelse(spatial_dist_long > 50000, -100000, 100000),
  spatial_dist_short = phyloseq::sample_data(shortphyseq) %>%
    with(x + 30000 * as.integer(site) + 100000 * as.integer(year)) %>%
    dist,
  spatial_dist_short2 = spatial_dist_short + ifelse(spatial_dist_short > 50000, -100000, 100000),
  unif_corr = vegan::mantel.correlog(unif_dist, spatial_dist,
                                     break.pts = 0:13 - 0.5,
                                     cutoff = FALSE),
  bc_corr_long = vegan::mantel.correlog(bc_dist_long, spatial_dist_long,
                                        break.pts = 0:13 - 0.5,
                                        cutoff = FALSE),
  bc_corr_short = vegan::mantel.correlog(bc_dist_short, spatial_dist_short,
                                         break.pts = 0:13 - 0.5,
                                         cutoff = FALSE),
  
  
  unif_corr2 = vegan::mantel.correlog(unif_dist, spatial_dist2,
                                      break.pts = 0:13 - 0.5,
                                      cutoff = FALSE),
  bc_corr_long2 = vegan::mantel.correlog(bc_dist_long, spatial_dist_long2,
                                         break.pts = 0:13 - 0.5,
                                         cutoff = FALSE),
  bc_corr_short2 = vegan::mantel.correlog(bc_dist_short, spatial_dist_short2,
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
                    .tag_in = step),
    format = "fst"),
  
  # seq_counts----
  # merge all the sequence counts into one data.frame
  seq_counts = target(
    dplyr::bind_rows(seq_count),
    transform = combine(seq_count),
    format = "fst"),
  
  # qstats----
  # calculate quality stats for the different regions
  qstats = target(
    q_stats(file_in(!!file.path(region_dir, region_file))),
    transform = map(.data = !!predada_meta,
                    .tag_in = step,
                    .id = c(well_ID, region)),
    format = "fst"),
  #qstats_join----
  # join all the quality stats into one data.frame
  qstats_join = target(
    dplyr::bind_rows(qstats) %>%
      tidyr::extract(
        col = "file",
        into = c("seq_run", "plate", "well", "direction", "region"),
        regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.trim\\.fastq\\.gz") %>%
      dplyr::left_join(datasets),
    transform = combine(qstats),
    format = "fst"),
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
  filter(!(step %in% c("raw", "combined")) | seq_run == '"pb_500"')
tictoc::toc()

if (!interactive()) saveRDS(plan, plan_file)

remove(itsx_meta)
remove(predada_meta)
remove(dada_meta)
remove(taxonomy_meta)
remove(regions)

flog.info("\nCalculating outdated targets...")
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

