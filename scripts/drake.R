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
  cmaln_file_long <- file.path(locarna_dir, "long_cmalign.aln")
  guide_tree_file <- file.path(locarna_dir, "32S_guide.tree")
  mlocarna_pp_dir <- file.path(locarna_dir, "consensus_pp")
  mlocarna_result_dir <- file.path(locarna_dir, "consensus")
  makelocarna <- file.path(r_dir, "snakemakelocarna.sh")
  makelocarna_profile <- file.path(lab_dir, "snakemakelocarnaUPPMAX")
  makelocarna_conda <- file.path(lab_dir, "conda", "snakemakelocarna.yaml")
  mlocarna_aln_file <- file.path(mlocarna_result_dir, "results", "result.stk")
  raxml_locarna_out_dir <- file.path(data_dir, "raxml_locarna")
  raxml_decipher_out_dir <- file.path(data_dir, "raxml_decipher")
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
  #longtree_file <- file.path(pasta_dir, "pasta_raxml.tree")
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
  cmaln_file_long <- snakemake@config$cmaln_long
  guide_tree_file <- snakemake@config$guide_tree
  mlocarna_aln_file <- snakemake@config$mlocarna_aln
  mlocarna_pp_dir <- snakemake@config$mlocarna_pp_dir
  mlocarna_result_dir <- snakemake@config$mlocarna_dir
  makelocarna <- snakemake@config$makelocarna
  makelocarna_profile <- snakemake@config$makelocarna_profile
  makelocarna_conda <- snakemake@config$makelocarna_conda
  raxml_locarna_out_dir <- snakemake@config$raxml_locarna_dir
  raxml_decipher_out_dir <- snakemake@config$raxml_decipher_dir
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
  #longtree_file <- snakemake@output[["pasta"]][["tree"]]
  logfile <- file(snakemake@log[[1]], open = "at")
} else {
  stop("Can't find Snakemake object in non-interactive session!")
}

library(magrittr)
library(tidyverse)
library(rlang)
library(glue)
library(drake)
library(assertr)
library(disk.frame)


#### read scripts and configs ####
# load various pipeline functions
source(file.path(r_dir, "dada.R"))
source(file.path(r_dir, "extract_regions.R"))
source(file.path(r_dir, "inferrnal.R"))
source(file.path(r_dir, "locarrna.R"))
source(file.path(r_dir, "lsux.R"))
source(file.path(r_dir, "mantel.R"))
source(file.path(r_dir, "parallel_helpers.R"))
source(file.path(r_dir, "plate_check.R"))
source(file.path(r_dir, "raxml.R"))
source(file.path(r_dir, "taxonomy.R"))
source(file.path(r_dir, "utils.R"))

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
         fastq = make.names(glue("fastq_{file_ID}")),
         join_derep = make.names(glue("join_derep_{primer_ID}")),
         join_derep_map = make.names(glue("join_derep_map_{primer_ID}")),
         itsxtrim = make.names(glue("itsxtrim_{primer_ID}")),
         lsux_pos = make.names(glue("lsux_pos_{primer_ID}"))) %T>%
         {flog.info("finished fourth mutate...")} %>%
  filter(file.exists(file.path(trim_dir, trim_file))) %>%
  filter(normalizePath(file.path(trim_dir, trim_file)) %in% normalizePath(in_files)) %>%
  filter(file.size(file.path(trim_dir, trim_file)) > 40) %>%
  mutate(md5 = tools::md5sum(file.path(trim_dir, trim_file))) %>%
  verify(file.path(trim_dir, trim_file) %in% in_files) %>%
  mutate_at(c("join_derep", "join_derep_map", "itsxtrim", "lsux_pos", "fastq"), syms)

#### predada_meta ####
# predada_meta has one row per region per well
# it is used in targets regionx, filter, and derep2
# The index is well_ID + region
# plate_ID is defined because it will be the index for dada_meta
flog.info("Making predada_meta.")
predada_meta <- select(itsx_meta, seq_run, regions, primer_ID) %>%
  unique() %>%
  mutate_at("regions", str_split, ",") %>%
  unnest(regions) %>%
  dplyr::rename(region = regions) %>%
  mutate(positions = glue("positions_{primer_ID}")) %>%
  mutate_at("positions", syms) %>%
  left_join(regions, by = "region")

if (!interactive()) {
  saveRDS(itsx_meta, itsx_meta_file)
  saveRDS(predada_meta, predada_meta_file)
}

# #### dada_meta ####
# # dada_meta has one row per region per sequencing run
# # it is used for targets err, dada, dadamap, seq_table, and nochim
# # The index is seq_run + region
# flog.info("Making dada_meta.")
# dada_meta <- predada_meta %>%
#   select(seq_run, region, range, range_ID, primer_ID) %>%
#   unique() %>%
#   arrange(seq_run, region) %>%
#   mutate(region_derep = glue("region_derep_{seq_run}_{range_ID}"),
#          error_model = glue("err_{seq_run}_5_8S"),
#          preconseq = glue("preconseq_{primer_ID}")) %>%
#   left_join(datasets %>% select(-regions), by = "seq_run") %>%
#   left_join(regions, by = c("region", "range")) %>%
#   mutate_at(c("region_derep", "error_model", "preconseq"), syms)
# if (!interactive()) {
#   saveRDS(dada_meta, dada_meta_file)
# }
# 
# #### region_meta ####
# # region_meta has one row per region
# # it is used in target big_fasta.
# # the index is region
# flog.info("Making region_meta.")
# region_meta <- regions %>%
#   select(region) %>%
#   mutate(big_seq_table = glue("big_seq_table_{region}"),
#          big_fasta_file = glue("{cluster_dir}/{region}.fasta.gz")) %>%
#   mutate_at("big_seq_table", syms) %>%
#   unique()
# 
# #### taxonomy_meta ####
# # taxonomy_meta has one row per region and reference DB
# # it is used in target taxon,
# # and mapped to target guilds_table
# # The index is tax_ID (Taxonomy ID)
# flog.info("Making taxonomy_meta.")
# taxonomy_meta <- dada_meta %>%
#   mutate_if(is.list, as.character) %>%
#   select_at(names(regions)) %>%
#   separate_rows(reference, sep = " *, *") %>%
#   separate(reference, c("reference", "refregion"), sep = "\\.") %>%
#   filter(!is.na(reference)) %>%
#   crossing(methods) %>%
#   mutate(methodfile = ifelse(method == "idtaxa", "sintax", method),
#          reference_file = glue("{ref_dir}/{reference}.{refregion}.{methodfile}.fasta.gz"),
#          refdb = glue("refdb_{method}_{reference}_{refregion}"),
#          tax_ID = glue("{region}_{reference}_{refregion}_{method}"),
#          big_seq_table = glue("big_seq_table_{region}")) %>%
#   arrange(tax_ID) %>%
#   mutate_at(c("big_seq_table", "refdb"), syms)
# if (!interactive()) {
#   saveRDS(taxonomy_meta, taxonomy_meta_file)
#   mutate_if(taxonomy_meta, is.list, as.character) %>%
#     readr::write_csv(taxonomy_meta_csv_file)
#   readr::write_lines(taxonomy_meta$tax_ID, tids_file)
# }
# 
# #### ref_meta ####
# # ref_meta has one row per (region specific) reference database.
# # it is used for target dbprep
# # and mapped to target classifier
# flog.info("Making ref_meta.")
# ref_meta <- select(taxonomy_meta, "reference", "refregion", "method", "reference_file") %>%
#   unique() %>%
#   mutate(ref_ID = glue("{reference}_{refregion}"))

#### drake plan ####
shard <- 1L:bigsplit

flog.info("\nbuilding plan...")
tictoc::tic()
plan <- drake_plan(
  file_meta = target(
    dplyr::filter(itsx_meta, .data[["primer_ID"]] == primer_ID),
    transform = map(primer_ID = !!unique(itsx_meta$primer_ID))
  ),
  derep1 = target(
    dada2::derepFastq(
      file.path(trim_dir, file_meta$trim_file),
      qualityType = "FastqQuality",
      verbose = TRUE,
      n = 1e4
    ) %>%
      # add the names of the sequences
      tzara::add_derep_names(file.path(trim_dir, file_meta$trim_file)) %>%
      # we will not run dada on these, so we don't need the quality
      # information.  Removing it drastically reduces the size in memory.
      magrittr::inset2("quals", NULL),
    transform = map(file_meta, .id = primer_ID),
    dynamic = map(file_meta, .trace = file_meta)
  ),
  trace1 = target(
    get_trace(paste0("file_meta_", primer_ID), derep1),
    transform = map(derep1, .id = primer_ID)
  ),
  
  join_derep = target(
    tzara::combine_derep(
      derep1,
      .data = dplyr::bind_rows(trace1)[ , c("seq_run", "plate", "well",
                                            "direction", "trim_file")]
    ),
    transform = map(derep1, trace1, .id = primer_ID),
    dynamic = group(derep1)
  ),
  
  join_derep_map = target(
    join_derep$map,
    transform = map(join_derep, .id = primer_ID),
    dynamic = map(join_derep)
  ),
  
  split_derep = target(
    drake_slice(join_derep$fasta, index = shard, slices = bigsplit),
    dynamic = cross(join_derep, shard),
    transform = map(join_derep, .id = primer_ID)
  ),
  
  lsux = target(
    LSUx(
      seq = split_derep,
      cm_5.8S = file_in(!!cm_5.8S),
      cm_32S = file_in(!!cm_32S),
      glocal = TRUE,
      ITS1 = TRUE,
      cpu = ignore(itsx_cpus)
    ) %>%
      dplyr::mutate_at("seq_name", as.integer) %>%
      as.data.frame(),
    transform = map(split_derep, .id = primer_ID),
    dynamic = map(split_derep),
    format = "fst"
  ),
  
  lsux_pos = target(
    combine_dynamic_diskframe(lsux),
    transform = map(lsux, .id = primer_ID),
    format = "diskframe"
  ),
  
  join_derep_by = target(
    join_derep_map %>%
      dplyr::group_by("seq_run", "plate", "well") %>%
      dplyr::group_indices(),
    transform = map(join_derep_map, .id = primer_ID),
    dynamic = map(join_derep_map)
  ),
  
  join_derep_groups = target(
    join_derep_map %>%
      dplyr::group_by("seq_run", "plate", "well") %>%
      dplyr::group_keys(),
    transform = map(join_derep_map, .id = primer_ID),
    dynamic = map(join_derep_map)
  ),
  
  positions = target(
    lsux_pos %>%
     dplyr::filter(seq_name %in% join_derep_map[[1]]$map) %>%
      as.data.frame(row.names = NULL) %>%
    dplyr::left_join(
      join_derep_map[[1]],
      .,
      by = c("map" = "seq_name")) %>%
      dplyr::rename(seq = seq.id) %>%
      gather_regions(),
    transform = map(lsux_pos, join_derep_map, join_derep_by, .id = primer_ID),
    dynamic = group(join_derep_map, .by = join_derep_by)
  ),
  
  # regionx----
  # cut the required regions out of each sequence
  regionx = target(
    if (nrow(positions) > 0) {
      pos <- dplyr::group_by(positions, trim_file)
      tzara::extract_region(
        seq = file.path(trim_dir, dplyr::group_keys(pos)$trim_file),
        region = region_start,
        region2 = region_end,
        positions = dplyr::group_split(pos)
      )
    } else {
      ShortRead::ShortReadQ()
    },
    transform = map(
      .data = !!predada_meta,
      .tag_in = step,
      .id = c(seq_run, region)
    ),
    dynamic = map(positions)
  ),
  
  # filter----
  # quality filter the sequences
  filter = target(
    filterReads(
      regionx, 
      maxLen = max_length,
      minLen = min_length,
      maxEE = max_ee
    ),
    transform = map(regionx, max_length, min_length, max_ee,
                    .id = c(seq_run, region)),
    dynamic = map(regionx)
  ),
  trace = TRUE
)

#   # raw ----
#   # Join the raw read info
#   raw = target(
#     tzara::summarize_sread(
#       list(regionx),
#       name = !!symbols_to_values(filter_file),
#       max_ee = 1.5 * !!eval(parse(text = unique(symbols_to_values(max_ee))))
#     ),
#     transform = combine(
#       regionx, filter_file, max_ee, seq_run, region, primer_ID,
#       .by = c(primer_ID, region),
#       .tag_in = step,
#       .id = c(primer_ID, region)
#     ),
#     format = "fst"
#   ),
#   
#   # derep2----
#   # Do a second round of dereplication on the filtered regions.
#   derep2 = target({
#       infile <- file_in(!!file.path(filter_dir, filter_file))
#       if (file.size(infile) > 50) {
#         out <- dada2::derepFastq(
#           fls = infile, n = 1e4,
#           qualityType = "FastqQuality",
#           verbose = TRUE
#         ) %>%
#           tzara::add_derep_names(infile)
#       } else {
#         out = NULL
#       }
#       result <- list()
#       result[[filter_file]] <- out
#       result
#     },
#     transform = map(
#       .data = !!predada_meta,
#       .tag_in = step,
#       .id = c(well_ID, range_ID)
#     )
#   ),
#   
#   # region_derep ----
#   # put together all the dereplicated files for each region for each plate.
#   region_derep = target(
#     purrr::compact(c(derep2)),
#     transform = combine(derep2, .by = c(seq_run, range_ID), .tag_in = step),
#   ),
#   
#   # err----
#   # Calibrate dada error models (on different size ranges)
#   err = target({
#     err.fun <- if (tech == "PacBio" ) dada2::PacBioErrfun else
#       dada2::loessErrfun
#     dada2::learnErrors(
#       fls = region_derep,
#       nbases = 1e9,
#       multithread = ignore(dada_cpus),
#       randomize = TRUE,
#       errorEstimationFunction = err.fun, 
#       HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
#       BAND_SIZE = band_size,
#       pool = eval(rlang::parse_expr(pool)),
#       verbose = TRUE,
#       qualityType = "FastqQuality"
#     )},
#     transform = map(
#       .data = !!filter(dada_meta, range_ID == "5_8S"),
#       .tag_in = step,
#       .id = c(seq_run, range_ID)
#     )
#   ),
#   
#   # dada----
#   # Run dada denoising algorithm (on different size ranges)
#   dada = target(
#     dada2::dada(
#       derep = region_derep, err = error_model,
#       multithread = ignore(dada_cpus), 
#       HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
#       BAND_SIZE = band_size,
#       pool = eval(rlang::parse_expr(pool))),
#     transform = map(
#       .data = !!dada_meta,
#       .tag_in = step,
#       .id = c(seq_run, range_ID)
#     )
#   ),
#   
#   # dadamap----
#   # Make maps from individual sequences in source files to dada ASVs
#   # (for each size range, if multiple)
#   dada_map = target(
#     tzara::dadamap(region_derep, dada),
#     transform = map(dada,
#                     .tag_in = step, .id = c(seq_run, range_ID))),
#   
#   # seq_table----
#   # Make a sample x ASV abundance matrix
#   # (for each size range, if multiple)
#   seq_table = target(
#     dada2::makeSequenceTable(dada) %>%
#       magrittr::extract(,nchar(colnames(.)) >= min_length_post &
#                           nchar(colnames(.)) <= max_length_post),
#     transform = map(dada, .tag_in = step, .id = c(seq_run, range_ID))),
#   
#   # nochim----
#   # Remove likely bimeras
#   nochim = target(
#     dada2::removeBimeraDenovo(
#       list(seq_table) %>%
#         purrr::map(tibble::as_tibble, rownames = "file") %>%
#         purrr::reduce(dplyr::full_join, by = "file") %>%
#         tibble::column_to_rownames("file") %>%
#         as.matrix,
#       method = "consensus",
#       multithread = ignore(dada_cpus)),
#     transform = combine(seq_table, .by = c(seq_run, region))),
#   
#   # big_seq_table ----
#   # Join all the sequence tables for each region
#   big_seq_table = target(
#     dada2::mergeSequenceTables(tables = list(nochim)),
#     transform = combine(nochim, .tag_in = step, .by = region)),
#   
#   # big_fasta ----
#   # write the big_seq_table as a fasta files so that they can be clustered by
#   # VSEARCH.
#   big_fasta = target(
#     write_big_fasta(big_seq_table,
#                     file_out(!!big_fasta_file)),
#     transform = map(.data = !!region_meta, .id = region)),
#   
#   # combined ----
#   # Replace raw reads which were mapped to a DADA2 ASV with the ASV,
#   # and put them all in one tibble.
#   combined = target(
#     tzara::combine_bigmaps(dplyr::bind_rows(dada_map),
#                     dplyr::bind_rows(raw)),
#     transform = combine(dada_map, raw, .by = primer_ID, .tag_in = step),
#     format = "fst"
#   ),
#   
#   # conseq ----
#   # Calculate consensus sequences for full ITS and LSU for each ITS2 ASV with
#   # with at least 3 sequences.
#   preconseq = target(
#     combined %>%
#       tidyr::extract(
#         name,
#         c("seq_run", "plate", "well", "region"),
#         "([pi][sb]_\\d{3})_(\\d{3})-([A-H]\\d{1,2})-(.+)\\.qfilt\\.fastq\\.gz"
#       ) %>%
#       tidyr::spread(key = "region", value = "seq") %>%
#       dplyr::group_by(ITS2) %>%
#       dplyr::filter(!is.na(ITS2), dplyr::n() >= 3) %>%
#       dplyr::mutate(primer_ID = primer_ID) %>%
#       region_concat("LSU", c("LSU1", "D1", "LSU2", "D2", "LSU3", "D3", "LSU4")) %>%
#       region_concat("ITS", c("ITS1", "5_8S", "ITS2")) %>%
#       region_concat("long", c("ITS", "LSU")) %>%
#       region_concat("short", c("5_8S", "ITS2", "LSU1")) %>%
#       region_concat("32S", c("5_8S", "ITS2", "LSU")) ,
#     transform = map(combined, .id = primer_ID)
#   ),
#   
#   conseq = target(
#     preconseq %>%
#       # for LSU1 and 5.8S, use only long reads if there are enough
#       # if there are not enough, use short reads
#       dplyr::mutate(
#         !! region := if (region %in% c("5_8S", "LSU1")) {
#         if (sum(seq_run == "pb_500" && !is.na(.data[[region]])) >= 3) {
#           ifelse(seq_run == "pb_500", .data[[region]], NA_character_)
#         } else {
#           ifelse(seq_run == "pb_500", NA_character_, .data[[region]])
#         }
#       } else {
#         .data[[region]]
#       }  ) %>%
#       dplyr::filter(sum(!is.na(.data[[region]])) >= 3) %>%
#       dplyr::summarize(
#         nreads = dplyr::n(),
#         !! region := as.character(
#           tzara::cluster_consensus(.data[[region]], seq.id, ignore(dada_cpus)))
#       ) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(
#         stringr::str_count(.[[region]], "[MRWSYKVHDBN]") < 3,
#         !is.na(.[[region]]),
#         nchar(.[[region]]) >= min_length
#       ),
#     transform = map(
#       .data = !!filter(
#         dada_meta,
#         region %in% c("ITS", "LSU", "32S", "long", "short")
#       ),
#       .id = c(primer_ID, region)
#     ),
#     format = "fst"
#   ),
#   
#   # allseqs ----
#   # make a data frame of all consensus sequences and ASVs; each row is an ITS2 ASV
#   allseqs = target(
#     make_allseq_table(list(conseq),
#                       drake_combine(big_seq_table)) %>%
#       dplyr::filter(!is.na(ITS2)) %>%
#       dplyr::mutate(full = dplyr::coalesce(long, short)),
#     transform = combine(conseq, big_seq_table),
#     format = "fst"
#   ),
#   
#   # dbprep ----
#   # import the RDP, Silva, and UNITE databases, and format them for use with
#   # DECIPHER.
#   refdb_dada2 = target(
#     file_in(reference_file),
#     transform = map(.data = !!filter(ref_meta, method == "dada2"),
#                     .tag_out = refdb,
#                     .id = ref_ID)),
#   
#   
#   refdb_sintax = target(
#     file_in(reference_file),
#     transform = map(.data = !!filter(ref_meta, method == "sintax"),
#                     .tag_out = refdb,
#                     .id = ref_ID)),
#   
#   refdb_idtaxa = target(
#     fit_idtaxa(file_in(reference_file)),
#     transform = map(.data = !!filter(ref_meta, method == "idtaxa"),
#                     .tag_out = refdb,
#                     .id = ref_ID)),
#   
#   # taxon ----
#   # Assign taxonomy to each ASV
#   taxon = target({
#     seqs <- dplyr::select(allseqs, region, "hash") %>%
#       dplyr::filter(complete.cases(.)) %>%
#       unique() %>%
#       {set_names(.[[region]], .$hash)}
#     taxonomy(seq = seqs,
#              reference = refdb,
#              method = method,
#              multithread = ignore(dada_cpus)) %>%
#       taxtable(names = names(seqs))
#     },
#     transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID),
#     format = "fst"),
#   
#   taxon_table = target(
#     drake_combine(taxon) %>%
#      combine_taxon_tables(allseqs),
#     transform = combine(taxon)
#   ),
#   
#   taxon_labels = target(
#     make_taxon_labels(taxon_table)
#   ),
#   
#   # funguild_db ----
#   # Download the FUNGuild database
#   funguild_db = FUNGuildR::get_funguild_db(),
#   
#   # guilds_table ----
#   # Assign ecological guilds to the ASVs based on taxonomy.
#   guilds_table = target(
#     FUNGuildR::funguild_assign(taxon, db = funguild_db),
#     transform = map(taxon, tax_ID, .tag_in = step, .id = tax_ID),
#     format = "fst"),
#   
#   # platemap ----
#   # Read the map between plate locations and samples
#   platemap = target(
#     read_platemap(file_in(!!platemap_file), platemap_sheet),
#     format = "fst"),
#   
#   # long_consensus----
#   # get the long amplicon consensus and convert to RNAStringSet.
#   # Use the sequence hashes as the name, so that names will be robust
#   # make sure that we only use sequences that actually match the long sequence.
#   cons_32S = allseqs %>%
#     dplyr::select("32S", ITS1, long, hash) %>%
#     dplyr::filter(
#       complete.cases(.),
#       endsWith(long, `32S`)
#     ) %>%
#     dplyr::select("32S", hash) %>%
#     unique() %>%
#     assertr::assert(assertr::is_uniq, hash) %$%
#     rlang::set_names(`32S`, hash) %>%
#     chartr(old = "T", new = "U") %>% 
#     Biostrings::RNAStringSet(),
#   cons_LSU = allseqs %>%
#     dplyr::select(LSU, ITS1, hash) %>%
#     dplyr::filter(
#       complete.cases(.),
#       endsWith(long, LSU)
#     ) %>%
#     dplyr::select(LSU, hash) %>%
#     unique() %$%
#     rlang::set_names(LSU, hash) %>%
#     chartr(old = "T", new = "U") %>%
#     Biostrings::RNAStringSet(),
#   cons_long = allseqs %>%
#     dplyr::select(long, ITS1, hash) %>%
#     dplyr::filter(complete.cases(.)) %>%
#     dplyr::select(long, hash) %>%
#     unique() %$%
#     rlang::set_names(long, hash) %>%
#     chartr(old = "T", new = "U") %>%
#     Biostrings::RNAStringSet(),
#   
#   # cmaln_32S
#   # Align conserved positions and annotate conserved secondary structure of the
#   # 32S consensus using Infernal.
#   cmaln_32S =
#     cmalign(cmfile = file_in(!!cm_32S),
#             seq = cons_32S,
#             glocal = TRUE,
#             cpu = ignore(dada_cpus)),
#   
#   # cmaln_long
#   # Concatenate ITS1 (padded with gaps) to the beginning of the Infernal 32S
#   # alignment to make a preliminary long amplicon alignment with conserved
#   # 5.8S/LSU secondary structure annotations.
#   cmaln_long =
#     dplyr::inner_join(
#       allseqs %>%
#         dplyr::select(hash, long, ITS1) %>%
#         dplyr::filter(complete.cases(.), startsWith(long, ITS1)) %>%
#         unique(),
#       tibble::tibble(
#         hash = cmaln_32S$names,
#         aln = as.character(cmaln_32S$alignment@unmasked)
#       ),
#       by = "hash"
#     ) %>%
#     dplyr::mutate(
#       ITS1 = chartr("T", "U", ITS1),
#       ITS1 = stringr::str_pad(ITS1, max(nchar(ITS1)), "right", "-"),
#       aln = stringr::str_c(ITS1, aln)
#     ) %>% {
#       write_clustalw_ss(
#         aln = magrittr::set_names(.$aln, .$hash) %>%
#           Biostrings::RNAStringSet(),
#         sec_str = paste0(
#           # The first four bases after the ITS1 primer are paired to other
#           # parts of SSU, so should not be paired in the rest of the amplicon.
#           # Otherwise, we have no structure information for the ITS1 region (it
#           # does not have eukaryote-wide or fungi-wide conserved structure)
#           stringr::str_pad(
#             "xxxx",
#             max(nchar(.$ITS1)),
#             "right",
#             "-"),
#           cmaln_32S$SS_cons
#         ),
#         # Don't take the alignment in the variable regions as
#         # conserved.
#         ref = stringr::str_pad(chartr("v", ".", cmaln_32S$RF),
#                                max(nchar(.$aln)), "left", "-"),
#         seq_names = .$hash,
#         file = file_out(!!cmaln_file_long)
#       )
#     },
#   
#   # guidetree_32S ----
#   # Take only the conserved, gap-free positions in the 32S alignment
#   aln_32S_conserv =
#     remove_nonconsensus_nongaps(cmaln_32S),
#   
#   # Generate a distance matrix from the conserved 32S alignment.
#   dist_32S_conserv =
#     DECIPHER::DistanceMatrix(
#       aln_32S_conserv,
#       type = "dist",
#       processors = ignore(dada_cpus)
#     ),
#   
#   # Make a UPGMA guide tree for mlocarna based on the aligned positions in
#   # the Infernal alignment for 32S.
#   guidetree_32S =
#     DECIPHER::IdClusters(
#       myDistMatrix = dist_32S_conserv,
#       collapse = -1, #break polytomies
#       type = "dendrogram",
#       processors = ignore(dada_cpus)
#     ) %T>%
#     DECIPHER::WriteDendrogram(
#       file = file_out(!!guide_tree_file)
#     ),
#   
#   # create the pair probability files for mlocarna
#   # this only needs to be redone if the sequences change, and it can be
#   # efficiently done in parallel
#   # It is made in a mirror directory, because mlocarna will write a file to the
#   # same directory.
#   mlocarna_pp_long = {
#     mlocarna_realign(
#       alignment = file_in(!!cmaln_file_long),
#       target_dir = file_out(!!mlocarna_pp_dir),
#       cpus = ignore(raxml_cpus),
#       only_dps = TRUE
#     )
#     mirror_dir(!!mlocarna_pp_dir, !!mlocarna_result_dir)
#   },
#   
#   # realign the consensus sequences using mlocarna
#   # this is a progressive alignment, where each locarna alignment is single
#   # threaded, and the ones closer to the root of the tree take a lot longer.
#   # it is efficient to do this in parallel near the beginning, but only one or
#   # a few cores can be utilized at the end; mlocarna itself does not execute
#   # them in parallel at all.
#   # snakemakelocarna.sh is a wrapper for locarna which writes empty files for
#   # every locarna call (thus spoofing mlocarna into thinking that it is working)
#   # while collecting the dependency structure of the calls.  On the last call, it
#   # submits all the calls as separate single-threaded jobs on SLURM using
#   # Snakemake.
#   # This requires a lot of configuration to work properly on a given cluster
#   # system; see files in makelocarna_profile.
#   # (default "config/snakemakelocarnaUPPMAX")
#   # This will always fail on the first run, because it does not wait for all of
#   # the locarna calls to complete.  Check the cluster queue and make sure
#   # all of the jobs have completed, then this rule should pass on the next run
#   # of drake.
#   realign_long = {
#     mlocarna_pp_long
#     file_out(!!mlocarna_aln_file)
#     mlocarna_realign(
#       alignment = file_in(!!cmaln_file_long),
#       guide_tree = file_in(!!guide_tree_file),
#       target_dir = file_out(!!mlocarna_result_dir),
#       stockholm = TRUE,
#       consensus_structure = "alifold",
#       cpus = 1,
#       skip_pp = TRUE,
#       pw_aligner = normalizePath(file_in(!!makelocarna)),
#       pw_aligner_options = paste(
#         "--profile", normalizePath(!!makelocarna_profile),
#         "--conda", normalizePath(!!makelocarna_conda)
#       )
#     )
#   },
#   
#   # read the mlocarna output
#   aln_locarna_long =
#     read_stockholm_msa(file_in(!!mlocarna_aln_file)),
#   
#   # make a tree based on the realigned consensus using RAxML
#   raxml_locarna_long = {
#       raxml_RNA(
#         RNAaln = aln_locarna_long$alignment,
#         S = aln_locarna_long$SS_cons,
#         m = "GTRGAMMA",
#         A = "S16",
#         f = "a",
#         N = "autoMRE_IGN",
#         p = 12345,
#         x = 827,
#         k = TRUE,
#         file = "locarna_long",
#         dir = !!raxml_locarna_out_dir,
#         exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
#         threads = ignore(raxml_cpus)
#       )
#     setwd(wd)
#     result
#   },
#   
#   aln_decipher_long =
#     allseqs %>%
#     dplyr::select(hash, long, ITS1) %>%
#     dplyr::filter(complete.cases(.), startsWith(long, ITS1)) %>%
#     unique() %$%
#     set_names(long, hash) %>%
#     chartr("T", "U", .) %>%
#     Biostrings::RNAStringSet() %>%
#     DECIPHER::AlignSeqs(iterations = 10,
#                         refinements = 10,
#                         processors = ignore(raxml_cpus)),
#   
#   aln_decipher_LSU =
#     allseqs %>%
#     dplyr::select(hash, long, LSU, ITS1) %>%
#     dplyr::filter(complete.cases(.), startsWith(long, ITS1)) %>%
#     unique() %$%
#     set_names(LSU, hash) %>%
#     dplyr::filter(!duplicated(LSU)) %>%
#     chartr("T", "U", .) %>%
#     Biostrings::RNAStringSet() %>%
#     DECIPHER::AlignSeqs(iterations = 10,
#                         refinements = 10,
#                         processors = ignore(raxml_cpus)),
#   
#   aln_decipher_LSU_trim = trim_LSU_intron(aln_decipher_LSU),
#   
#   raxml_decipher_LSU = {
#     if (!dir.exists(!!raxml_decipher_out_dir))
#       dir.create(!!raxml_decipher_out_dir, recursive = TRUE)
#     wd <- setwd(!!raxml_decipher_out_dir)
#     result <- 
#       aln_decipher_LSU_trim %>%
#       Biostrings::DNAStringSet() %>%
#       ape::as.DNAbin() %>%
#       as.matrix() %>%
#       ips::raxml(
#         DNAbin = .,
#         m = "GTRGAMMA",
#         f = "a",
#         N = "autoMRE_IGN",
#         p = 12345,
#         x = 827,
#         k = TRUE,
#         file = "decipher_LSU",
#         exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
#         threads = ignore(raxml_cpus))
#     setwd(wd)
#     result
#   },
#   
#   raxml_decipher_long = {
#     if (!dir.exists(!!raxml_decipher_out_dir))
#       dir.create(!!raxml_decipher_out_dir, recursive = TRUE)
#     wd <- setwd(!!raxml_decipher_out_dir)
#     result <- 
#       aln_decipher_long %>%
#       Biostrings::DNAStringSet() %>%
#       ape::as.DNAbin() %>%
#       as.matrix() %>%
#       ips::raxml(
#         DNAbin = .,
#         m = "GTRGAMMA",
#         f = "a",
#         N = "autoMRE_IGN",
#         p = 12345,
#         x = 827,
#         k = TRUE,
#         file = "decipher_long",
#         exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
#         threads = ignore(raxml_cpus))
#     setwd(wd)
#     result
#     },
#   
#   labeled_tree_decipher_long =
#     relabel_tree(
#       tree = raxml_decipher_long$bipartitions,
#       old = taxon_labels$label,
#       new = taxon_labels$tip_label
#     ) %T>%
#     castor::write_tree(file_out("data/trees/decipher_long.tree")),
#   
#   labeled_tree_mlocarna_long =
#     relabel_tree(
#       tree = raxml_mlocarna_long$bipartitions,
#       old = taxon_labels$label,
#       new = taxon_labels$tip_label
#     ) %T>%
#     castor::write_tree(file_out("data/trees/mlocarna_long.tree")),
#   
#   
#   # Tulasnella is on a long branch because of a high divergence rate.
#   # it ends up with Metazoa in the tree.  Exclude it.
#   tulasnella =
#     taxon_table %>%
#     dplyr::filter(
#       rank == "family",
#       taxon == "Tulasnellaceae"
#     ) %$%
#     label %>%
#     unique(),
#   
#   # find the sequences that were confidently identified as fungi
#   # (and placed in a phylum).  We extract the common ancestor of these
#   surefungi =
#     taxon_table %>%
#     dplyr::group_by(label) %>%
#     dplyr::filter(
#       rank == "phylum",
#       n_tot >= 6,
#       n_diff == 1,
#       grepl("mycota", taxon)
#     ) %$%
#     label %>%
#     unique() %>%
#     setdiff(tulasnella),
#   
#   # find the most abundant sequence which we are confident is a plant
#   bestplant =
#     taxon_table %>%
#     dplyr::ungroup() %>%
#     dplyr::filter(rank == "phylum",
#                   taxon == "Chlorophyta",
#                   n_tot == 6,
#                   n_diff == 1) %>%
#     dplyr::filter(n_reads == max(n_reads)) %$%
#     label[1],
#   
#   # root the tree within plants.  This isn't accurate, but it should guarantee
#   # that it is not rooted within fungi.
#   tree_decipher_long = 
#     ape::root(raxml_decipher_long$bipartitions,
#               bestplant),
#   
#   # Extract the minimally inclusive clade including all confidently identified
#   # fungi (except Tulasnella).
#   fungi_tree_decipher_long =
#     ape::getMRCA(
#       phy = tree_decipher_long,
#       tip = intersect(surefungi, tree_decipher_long$tip.label)
#     ) %>%
#     ape::extract.clade(phy = tree_decipher_long),
#   
#   taxon_phy = phylotax(
#     tree = fungi_tree_decipher_long,
#     taxa = taxon_table
#   ),
#   
#   physeq = assemble_physeq(
#     platemap,
#     datasets,
#     relabel_seqtable(big_seq_table_ITS2),
#     fungi_tree_decipher_long),
# 
#   longphyseq = physeq %>%
#     phyloseq::prune_samples(
#       samples = phyloseq::sample_data(.)$dataset == "long-pacbio"
#       & phyloseq::sample_data(.)$sample_type == "Sample"
#       & (phyloseq::sample_data(.)$qual == "X" | phyloseq::sample_data(.)$year == "2015")
#       & rowSums(phyloseq::otu_table(.)) > 0
#       ),
#   
#   shortphyseq = physeq %>%
#     phyloseq::prune_samples(
#       samples = phyloseq::sample_data(.)$dataset == "short-pacbio"
#       & phyloseq::sample_data(.)$sample_type == "Sample"
#       & (phyloseq::sample_data(.)$qual == "X" | phyloseq::sample_data(.)$year == "2015")
#       & rowSums(phyloseq::otu_table(.)) > 0),
#   
#   unif_dist = phyloseq::distance(longphyseq, "unifrac"),
#   bc_dist_long = phyloseq::distance(longphyseq, "bray"),
#   bc_dist_short = phyloseq::distance(shortphyseq, "bray"),
#   spatial_dist_long = phyloseq::sample_data(longphyseq) %>%
#     with(x + 30000 * as.integer(site) + 100000 * as.integer(year)) %>%
#     dist,
#   spatial_dist_long2 = spatial_dist_long + ifelse(spatial_dist_long > 50000, -100000, 100000),
#   spatial_dist_short = phyloseq::sample_data(shortphyseq) %>%
#     with(x + 30000 * as.integer(site) + 100000 * as.integer(year)) %>%
#     dist,
#   spatial_dist_short2 = spatial_dist_short + ifelse(spatial_dist_short > 50000, -100000, 100000),
#   unif_corr = vegan::mantel.correlog(unif_dist, spatial_dist,
#                                      break.pts = 0:13 - 0.5,
#                                      cutoff = FALSE),
#   bc_corr_long = vegan::mantel.correlog(bc_dist_long, spatial_dist_long,
#                                         break.pts = 0:13 - 0.5,
#                                         cutoff = FALSE),
#   bc_corr_short = vegan::mantel.correlog(bc_dist_short, spatial_dist_short,
#                                          break.pts = 0:13 - 0.5,
#                                          cutoff = FALSE),
#   
#   
#   unif_corr2 = vegan::mantel.correlog(unif_dist, spatial_dist2,
#                                       break.pts = 0:13 - 0.5,
#                                       cutoff = FALSE),
#   bc_corr_long2 = vegan::mantel.correlog(bc_dist_long, spatial_dist_long2,
#                                          break.pts = 0:13 - 0.5,
#                                          cutoff = FALSE),
#   bc_corr_short2 = vegan::mantel.correlog(bc_dist_short, spatial_dist_short2,
#                                           break.pts = 0:13 - 0.5,
#                                           cutoff = FALSE),
#   
#   # seq_count ----
#   # Count the sequences in each fastq file
#   seq_count = target(
#     tibble::tibble(
#       file = file_in(fq_file),
#       reads = system(glue::glue("zcat {file} | grep -c '^@'", file = fq_file),
#                      intern = TRUE) %>%
#         as.integer),
#     transform = map(fq_file = !!c(file.path(trim_dir, itsx_meta$trim_file),
#                                   file.path(region_dir, predada_meta$region_file),
#                                   file.path(filter_dir, predada_meta$filter_file)),
#                     .tag_in = step),
#     format = "fst"),
#   
#   # seq_counts----
#   # merge all the sequence counts into one data.frame
#   seq_counts = target(
#     dplyr::bind_rows(seq_count),
#     transform = combine(seq_count),
#     format = "fst"),
#   
#   # qstats----
#   # calculate quality stats for the different regions
#   qstats = target(
#     q_stats(file_in(!!file.path(region_dir, region_file))),
#     transform = map(.data = !!predada_meta,
#                     .tag_in = step,
#                     .id = c(well_ID, region)),
#     format = "fst"),
#   #qstats_join----
#   # join all the quality stats into one data.frame
#   qstats_join = target(
#     dplyr::bind_rows(qstats) %>%
#       tidyr::extract(
#         col = "file",
#         into = c("seq_run", "plate", "well", "direction", "region"),
#         regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.trim\\.fastq\\.gz") %>%
#       dplyr::left_join(datasets),
#     transform = combine(qstats),
#     format = "fst"),
#   # qstats_knit----
#   # knit a report about the quality stats.
#   qstats_knit = {
#     if (!dir.exists(!!out_dir)) dir.create(!!out_dir, recursive = TRUE)
#     rmarkdown::render(
#       knitr_in(!!file.path(rmd_dir, "qual-check.Rmd")),
#       output_file = file_out(!!file.path(out_dir, "qual-check.pdf")),
#       output_dir = !!out_dir)},
#   # max_expand = if (interactive()) 9 else NULL,
#   trace = TRUE
# )
tictoc::toc()

if (!interactive()) saveRDS(plan, plan_file)

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

remove(snakemake)

save(list = ls(),
     file = drakedata_file)
