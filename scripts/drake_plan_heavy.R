suppressPackageStartupMessages({
  library(futile.logger)
  library(magrittr)
  library(tidyverse)
  library(rlang)
  library(glue)
  library(drake)
  library(assertr)
  library(disk.frame)
})

if (interactive()) {
  flog.info("Creating drake plan in interactive session...")
  config <- yaml::read_yaml("config/config.yaml")
  config$bigsplit <- config$smallsplit
} else if (exists("snakemake")) {
  flog.info("Creating drake plan in snakemake session...")
  snakemake@source(".Rprofile", echo = FALSE)
  config <- snakemake@config
} else {
  stop("Can't find Snakemake object in non-interactive session!")
}


#### read scripts and configs ####
# load various pipeline functions

source(file.path(config$rdir, "dada.R"))
source(file.path(config$rdir, "qstats.R"))
source(file.path(config$rdir, "inferrnal.R"))
source(file.path(config$rdir, "locarrna.R"))
source(file.path(config$rdir, "lsux.R"))
source(file.path(config$rdir, "mantel.R"))
source(file.path(config$rdir, "parallel_helpers.R"))
source(file.path(config$rdir, "plate_check.R"))
source(file.path(config$rdir, "raxml.R"))
source(file.path(config$rdir, "taxonomy.R"))
source(file.path(config$rdir, "utils.R"))
source(file.path(config$rdir, "epa-ng.R"))
source(file.path(config$rdir, "AlignSeqs.R"))
source(file.path(config$rdir, "DECIPHER-utils.R"))
source(file.path(config$rdir, "epa_iterate.R"))
source(file.path(config$rdir, "mafft.R"))

setup_log("drakeplan")

# load the dataset and region definitions
flog.info("Loading datasets file.")
datasets <- read_csv(config$dataset, col_types = "cccccccicccccccccccicc")
flog.info("Loading regions file.\n")
regions <- read_csv(config$regions, col_types = "cccciiiiic") %>%
  mutate_at("range", replace_na, "") %>%
  mutate(min_length_post = coalesce(min_length_post, min_length),
         max_length_post = coalesce(max_length_post, max_length))
flog.info("Loading taxonomic assignment methods file.")
methods <- read_csv(config$methods, col_types = "ci")

flog.info("Loading table of demultiplexed files.")
demux_meta <- read_delim(
  config$demux_file,
  delim = " ",
  col_names = c("md5", "trim_file"),
  col_types = "cc",
  trim_ws = TRUE
)


# "meta" dataframes containing definitions for each instance of mapped targets

#### itsx_meta ####
# itsx_meta has one row per demultiplexed, primer-trimmed fastq.gz file
# Illumina files need to be treated differently, so they are not included.
flog.info("Making itsx_meta.")
itsx_meta <- datasets %>%
  filter(tech != "Illumina") %>%
  mutate(primer_ID = paste0(str_replace(forward, "_tag.*$", ""),
                            str_replace(reverse, "_tag.*$", "")),
         plate = map(runs, ~formatC(seq(.), width = 3, flag = "0"))) %>%
  unnest(plate) %>%
  mutate(direction = if_else(tech == "PacBio",
                             list(c("f", "r")),
                             list(""))) %>%
  unnest(direction) %>%
  mutate(well = list(tidyr::crossing(Row = LETTERS[1:8], Col = 1:12) %>%
                       transmute(well = paste0(Row, Col)))) %>%
  unnest(well) %>%
  mutate(trim_file = glue("{config$trimdir}/{seq_run}_{plate}/{seq_run}_{plate}-{well}{direction}.trim.fastq.gz") %>%
           unclass()) %>%
  inner_join(demux_meta, by = "trim_file") %>%
  filter(file.size(trim_file) > 40)

#### positions_meta ####
flog.info("Making positions_meta.")
positions_meta <- select(itsx_meta, "seq_run", "primer_ID") %>%
  unique() %>%
  mutate(lsux_pos = glue("lsux_pos_{primer_ID}"),
         derep_map = glue("derep_map_{primer_ID}")) %>%
  mutate_at(c("lsux_pos", "derep_map"), syms)

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
  mutate(positions = glue("positions_{seq_run}")) %>%
  mutate_at("positions", syms) %>%
  left_join(regions, by = "region")

#### dada_meta ####
# dada_meta has one row per region per sequencing run
# it is used for targets err, dada, dadamap, seq_table, and nochim
# The index is seq_run + region
# Illumina runs are treated separately.
flog.info("Making dada_meta.")
dada_meta <- predada_meta %>%
  select(seq_run, region, reference, primer_ID, positions, min_length_post,
         max_length_post) %>%
  unique() %>%
  arrange(seq_run, region) %>%
  left_join(datasets %>% select(-regions), by = "seq_run") %>%
  mutate(derep2 = glue("derep2_{seq_run}_{region}"),
         derep_groups = glue("derep_groups_{seq_run}"),
         error_model = glue("err_{seq_run}_{err_region}"),
         preconseq = glue("preconseq_{primer_ID}"),
         position_map = glue("position_map_{seq_run}")) %>%
  mutate_at(
    c("derep2", "derep_groups", "error_model", "preconseq", "position_map"),
    syms
    )

illumina_meta <- datasets %>%
  filter(tech == "Illumina") %>%
  mutate(
    primer_ID = paste0(str_replace(forward, "_tag.*$", ""),
                            str_replace(reverse, "_tag.*$", "")),
         plate = map(runs, ~formatC(seq(.), width = 3, flag = "0")),
         region = err_region
         ) %>%
  unnest(plate) %>%
  tidyr::expand_grid(
    read = c("R1", "R2"),
    row = LETTERS[1:8],
    col = 1:12,
    direction = c("f", "r")
    ) %>%
  mutate(
    well = paste0(row, col),
    trim_file = glue("{config$trimdir}/{seq_run}_{plate}/{seq_run}_{plate}-{well}_{read}{direction}.trim.fastq.gz") %>%
           unclass()) %>%
  inner_join(demux_meta, by = "trim_file") %>%
  pivot_wider(names_from = "read", values_from = c("trim_file", "md5")) %>%
  mutate(ID = glue::glue("{seq_run}_{plate}_{well}_{direction}_{region}"))

illumina_err_meta <- illumina_meta %>%
  select(seq_run, primer_ID, region, tech, hgp, band_size, pool) %>%
  mutate(
    derep_illumina = glue::glue("derep_illumina_{seq_run}"),
    illumina_id = glue::glue("illumina_id_{seq_run}")
  ) %>%
  mutate_at(c("derep_illumina", "illumina_id"), purrr::compose(syms, make.names)) %>% 
  unique() %>%
  expand_grid(read = c("R1", "R2"))

merge_meta <- illumina_err_meta %>%
  mutate(
    derep = make.names(glue::glue("derep_illumina_{seq_run}")),
    dada = make.names(glue::glue("dada_illumina_{seq_run}_{read}"))
  ) %>%
  pivot_wider(names_from = "read", values_from = c("dada"), names_prefix = "dada_") %>%
  mutate_at(c("derep", "dada_R1", "dada_R2"), syms) %>%
  dplyr::select(-region)

illumina_region_meta <- datasets %>%
  dplyr::filter(tech == "Illumina") %>%
  mutate(
    primer_ID = paste0(str_replace(forward, "_tag.*$", ""),
                       str_replace(reverse, "_tag.*$", ""))
  ) %>%
  tidyr::separate_rows(regions) %>%
  dplyr::rename(region = regions) %>%
  left_join(regions, by = "region") %>%
  mutate(
    seq_table_illumina = glue::glue("seq_table_illumina_{seq_run}"),
    lsux_illumina = glue::glue("lsux_illumina_{seq_run}")
  ) %>%
  mutate_at(c("seq_table_illumina", "lsux_illumina"), purrr::compose(syms, make.names))
  
flog.info("Making err_meta.")
err_meta <- dada_meta %>%
  filter(region == err_region) %>%
  select(seq_run, region, derep2, tech, hgp, band_size, pool)

flog.info("Making conseq_meta.")
conseq_meta <- dada_meta %>%
  filter(region %in% c("ITS1", "ITS", "LSU", "32S", "long", "short")) %>%
  select(region, primer_ID, preconseq) %>%
  left_join(regions %>% select(region, min_length), by = "region") %>%
  unique() %>%
  left_join(
      select(illumina_meta, seq_run_illumina = seq_run, primer_ID) %>%
        unique(),
      by = "primer_ID"
  ) %>%
  mutate(preconseq_illumina = ifelse(is.na(seq_run_illumina), list(NULL),
                                     glue("preconseq_illumina_{seq_run_illumina}") %>%
                                       make.names %>%
                                       syms))

#### region_meta ####
# region_meta has one row per region
# it is used in target big_fasta.
# the index is region
flog.info("Making region_meta.")
region_meta <- dada_meta %>%
  select(region) %>%
  unique() %>%
  mutate(big_seq_table = glue("big_seq_table_{region}"),
         big_fasta_file = glue("{config$clusterdir}/{region}.fasta.gz")) %>%
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
  select(region, reference) %>%
  separate_rows(reference, sep = " *, *") %>%
  separate(reference, c("reference", "refregion"), sep = "\\.") %>%
  filter(!is.na(reference)) %>%
  crossing(methods) %>%
  mutate(methodfile = ifelse(method == "idtaxa", "sintax", method),
         reference_file = glue("{config$ref_root}/{reference}.{refregion}.{methodfile}.fasta.gz"),
         refdb = glue("refdb_{method}_{reference}_{refregion}"),
         tax_ID = glue("{region}_{reference}_{refregion}_{method}"),
         big_seq_table = glue("big_seq_table_{region}")) %>%
  arrange(tax_ID) %>%
  mutate_at(c("big_seq_table", "refdb"), syms)

#### ref_meta ####
# ref_meta has one row per (region specific) reference database.
# it is used for target dbprep
# and mapped to target classifier
flog.info("Making ref_meta.")
ref_meta <- select(taxonomy_meta, "reference", "refregion", "method", "reference_file") %>%
  unique() %>%
  mutate(ref_ID = glue("{reference}_{refregion}"))

#### drake plan ####
flog.info("Assembling plan...")
tictoc::tic()
plan <- drake_plan(
  
  #### Plan section 1: Split reads into regions ####
  shard = 1L:!!config$bigsplit,
  file_meta = target(
    dplyr::filter(itsx_meta, .data[["primer_ID"]] == primer_ID) %>%
      dplyr::sample_n(nrow(.)),
    transform = map(primer_ID = !!unique(itsx_meta$primer_ID))
  ),
  derep1 = target(
    dada2::derepFastq(
      file_meta$trim_file,
      qualityType = "FastqQuality",
      verbose = TRUE,
      n = 1e4
    ) %>%
      # add the names of the sequences
      tzara::add_derep_names(file_meta$trim_file) %>%
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
      readd(derep1, cache = ignore(drake::drake_cache(cache_dir))),
      .data = dplyr::bind_rows(trace1)[ , c("seq_run", "plate", "well",
                                            "direction", "trim_file")]
    ),
    transform = map(derep1, trace1, .id = primer_ID)
  ),
  
  derep_map = target(
    join_derep$map,
    transform = map(join_derep, .id = primer_ID)
  ),
  
  derep_submap = target(
    dplyr::filter(derep_map, .data[["seq_run"]] == seq_run),
    transform = map(.data = !!positions_meta, .id = seq_run)
  ),
  
  lsux = target(
    drake_slice(join_derep$fasta, index = shard, slices = !!config$bigsplit) %>%
    LSUx(
      cm_5.8S = file_in(!!config$cm_5_8S),
      cm_32S = file_in(!!config$cm_32S),
      glocal = TRUE,
      ITS1 = TRUE,
      cpu = ignore(ncpus)
    ) %>%
      dplyr::mutate_at("seq_name", as.integer) %>%
      as.data.frame(),
    transform = map(join_derep, .id = primer_ID),
    dynamic = map(shard),
    format = "fst"
  ),
  
  lsux_pos = target(
    combine_dynamic_diskframe(lsux, cache = ignore(cache_dir)),
    transform = map(lsux, .id = primer_ID),
    format = "diskframe"
  ),
  
  derep_by = target(
    derep_submap %>%
      dplyr::group_by_at(c("seq_run", "plate", "well")) %>%
      dplyr::group_indices(),
    transform = map(derep_submap, .id = seq_run)
  ),
  
  derep_groups = target(
    derep_submap %>%
      dplyr::group_by_at(c("seq_run", "plate", "well", "trim_file")) %>%
      dplyr::group_keys(),
    transform = map(derep_submap, .id = seq_run)
  ),
  
  positions = target({
    ids <- unique(derep_submap$map)
    lsux_pos %>%
      dplyr::filter(seq_name %in% ids) %>%
      as.data.frame(row.names = NULL) %>%
      dplyr::left_join(
        derep_submap,
        .,
        by = c("map" = "seq_name")) %>%
      dplyr::rename(seq = seq.id) %>%
      gather_regions()},
    transform = map(lsux_pos, derep_submap, derep_by, .id = seq_run),
    dynamic = group(derep_submap, .by = derep_by, .trace = derep_by)
  ),
  
  position_map = target(
    dplyr::select(positions, "seq_run", "plate", "well") %>%
      unique(),
    transform = map(positions, .id = seq_run),
    dynamic = map(positions)
  ),
  
  # extract regions, do quality filtering, and dereplicate
  # this is all one step because the extraction and filtering are fast,
  # but saving the outputs would be large.
  derep2 = target(
    extract_and_derep(
      positions = positions,
      trim_file = unique(positions$trim_file),
      region = region,
      region_start = region_start,
      region_end = region_end,
      max_length = max_length,
      min_length = min_length,
      max_ee = max_ee
    ),
    transform = map(
      .data = !!select(predada_meta, positions, region_start, region_end,
                       max_length, min_length, max_ee, seq_run, region, primer_ID),
      .tag_in = step,
      .id = c(seq_run, region)
    ),
    dynamic = map(positions)
  ),

  derep2_length = target(
    length(derep2),
    transform = map(derep2, primer_ID, .tag_in = step, .id = c(seq_run, region)),
    dynamic = map(derep2)
  ),
  
  illumina_group = target(
    dplyr::filter(
      illumina_meta,
      .data[["seq_run"]] == seq_run
    ),
    transform = map(
      seq_run = !!unique(illumina_meta$seq_run),
      .tag_in = step,
      .id = seq_run
    )
  ),
  
  illumina_id = target(
    illumina_group,
    transform = map(illumina_group, .id = seq_run),
    dynamic = map(illumina_group)
  ),
  
  derep_illumina = target(
    filter_and_derep_pairs(
      trim_file_1 = illumina_group$trim_file_R1,
      trim_file_2 = illumina_group$trim_file_R2,
      trimR = 0,
      truncQ = 10,
      max_length = 2999,
      min_length = 50,
      max_ee = 3,
      ID = illumina_group$ID
    ),
    transform = map(
      illumina_group,
      .tag_in = step,
      .tag_out = derep_all,
      .id = seq_run
    ),
    dynamic = map(illumina_group)
  ),
  
  #### Plan section 2: Denoise using DADA ####

  # Calibrate dada error models on 5.8S
  err = target({
    err.fun <- if (tech == "PacBio" ) dada2::PacBioErrfun else
      dada2::loessErrfun
    dereps <- readd(derep2, cache = ignore(drake::drake_cache(cache_dir))) %>%
      purrr::compact()
    dada2::learnErrors(
      fls = dereps,
      nbases = 1e9,
      multithread = ignore(ncpus),
      randomize = TRUE,
      errorEstimationFunction = err.fun,
      HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
      BAND_SIZE = band_size,
      MAX_CONSIST = 50,
      pool = eval(rlang::parse_expr(pool)),
      verbose = TRUE,
      qualityType = "FastqQuality"
    )},
    transform = map(
      .data = !!err_meta,
      .tag_in = step,
      .id = c(seq_run, region)
    )
  ),
  
  err_illumina = target({
    dereps <- readd(derep_illumina, cache = ignore(drake::drake_cache(cache_dir))) %>%
      purrr::map(read) %>%
      purrr::compact()
    dada2::learnErrors(
      fls = dereps,
      nbases = 1e9,
      multithread = ignore(ncpus),
      randomize = TRUE,
      errorEstimationFunction = dada2::loessErrfun,
      HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
      BAND_SIZE = band_size,
      MAX_CONSIST = 50,
      pool = eval(rlang::parse_expr(pool)),
      verbose = TRUE,
      qualityType = "FastqQuality"
    )},
    transform = map(
      .data = !!illumina_err_meta,
      .tag_in = step,
      .id = c(seq_run, read)
    )
  ),
  
  # Run dada denoising algorithm (on different regions)
  dada = target(
    readd(derep2, cache = ignore(drake::drake_cache(cache_dir))) %>%
      set_names(
        readd(position_map, cache = ignore(drake::drake_cache(cache_dir))) %>%
          dplyr::bind_rows() %>%
          dplyr::mutate(region = region) %>%
          glue::glue_data("{seq_run}_{plate}_{well}_{region}")
      ) %>%
      robust_dada(
        derep = .,
        err = error_model,
        multithread = ignore(ncpus),
        HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
        BAND_SIZE = band_size,
        pool = eval(rlang::parse_expr(pool))),
    transform = map(
      .data = !!select(dada_meta, derep2, derep_groups, position_map, region, hgp, band_size,
                       pool, seq_run, primer_ID, error_model, min_length_post,
                       max_length_post),
      .tag_in = step,
      .id = c(seq_run, region)
    )
  ),
  
  dada_illumina = target({
    derep <- readd(derep_illumina, cache = ignore(drake::drake_cache(cache_dir))) %>%
      purrr::map(read) %>%
      set_names(
        dplyr::bind_rows(readd(illumina_id, cache = ignore(drake::drake_cache(cache_dir)))) %>%
          dplyr::mutate(read = read) %>%
          glue::glue_data("{seq_run}_{plate}_{well}_{read}_{direction}") %>%
          unique()
      )
      gc()
      robust_dada(
        derep = derep,
        err = err_illumina,
        multithread = ignore(ncpus),
        HOMOPOLYMER_GAP_PENALTY = eval(rlang::parse_expr(hgp)),
        BAND_SIZE = band_size,
        pool = eval(rlang::parse_expr(pool)))
    },
    transform = map(
      derep_illumina,
      err_illumina,
      illumina_id,
      hgp,
      band_size,
      pool,
      read,
      .tag_in = step,
      .id = c(seq_run, read)
    )
  ),
  
  merge = target(
    dada2::mergePairs(
      dadaF = purrr::compact(dada_R1),
      derepF = set_names(purrr::map(readd(derep, cache = ignore(drake::drake_cache(cache_dir))), "R1"), names(dada_R1)) %>% purrr::compact(),
      dadaR = purrr::compact(dada_R2),
      derepR = set_names(purrr::map(readd(derep, cache = ignore(drake::drake_cache(cache_dir))), "R2"), names(dada_R2)) %>% purrr::compact()
    ),
    transform = map(.data = !!merge_meta, .id = seq_run)
  ),

  # Make a sample x ASV abundance matrix
  # (for each size range, if multiple)
  seq_table1 = target(
    dada2::makeSequenceTable(purrr::compact(dada)) %>%
      magrittr::extract(,nchar(colnames(.)) >= min_length_post &
                          nchar(colnames(.)) <= max_length_post,
                        drop = FALSE),
    transform = map(dada, .tag_in = step, .tag_out = seq_table, .id = c(seq_run, region))
  ),
  
  seq_table_illumina = target({
    tab <- dada2::makeSequenceTable(merge)
    tabF <- tab[endsWith(rownames(tab), "_f"),]
    tabF <- tabF[, colSums(tabF) > 0]
    rownames(tabF) <- sub("_R._f$", "", rownames(tabF))
    tabF <- tibble::as_tibble(tabF, rownames = "name") %>%
      tidyr::pivot_longer(cols = -"name", names_to = "seq", values_to = "nread") %>%
      dplyr::filter(nread > 0)
    tabR <- tab[endsWith(rownames(tab), "_r"),]
    tabR <- tabR[, colSums(tabR) > 0]
    rownames(tabR) <- sub("_R._r$", "", rownames(tabR))
    colnames(tabR) <- dada2:::rc(colnames(tabR))
    tabR <- tibble::as_tibble(tabR, rownames = "name") %>%
      tidyr::pivot_longer(cols = -"name", names_to = "seq", values_to = "nread") %>%
      dplyr::filter(nread > 0)
    dplyr::bind_rows(tabF, tabR) %>%
      dplyr::group_by(name, seq) %>%
      dplyr::summarise(nread = sum(nread)) %>%
      tidyr::pivot_wider(names_from = "seq", values_from = "nread", values_fill = list(nread = 0L)) %>%
      tibble::column_to_rownames("name") %>%
      as.matrix()
  },
  transform = map(merge, .id = seq_run)
  ),
  
  lsux_illumina = target(
    LSUx(
      colnames(seq_table_illumina) %>% set_names(seq_along(.)),
      cm_5.8S = file_in(!!config$cm_5_8S),
      cm_32S = file_in(!!config$cm_32S),
      glocal = TRUE,
      ITS1 = FALSE,
      cpu = ignore(ncpus)
    ) %>%
      dplyr::rename(seq = seq_name) %>%
      gather_regions(),
    transform = map(seq_table_illumina, .id = seq_run)
  ),
  
  region_illumina = target(
    colnames(seq_table_illumina) %>%
      Biostrings::DNAStringSet() %>%
      ShortRead::ShortRead(Biostrings::BStringSet(as.character(seq_along(.)))) %>%
      tzara::extract_region(
        positions = lsux_illumina,
        region = region_start,
        region2 = region_end
      ),
    transform = map(.data = !!illumina_region_meta, .id = c(seq_run, region))
  ),
  
  region_table_illumina = target(
    dplyr::left_join(
      seq_table_illumina %>%
        t() %>%
        tibble::as_tibble(rownames = NULL) %>%
        dplyr::mutate(seq = as.character(seq.int(nrow(.)))) %>%
        tidyr::pivot_longer(cols = -"seq", names_to = "name", values_to = "nread") %>%
        dplyr::filter(nread > 0),
      region_illumina %>% {
        tibble::tibble(value = as.character(.@sread), seq = as.character(.@id))
      },
      by = "seq"
    ) %>%
      dplyr::group_by(value, name) %>%
      dplyr::summarize(nread = sum(nread)) %>%
      tidyr::pivot_wider(names_from = "value", values_from = "nread", values_fill = list(nread = 0L)) %>%
      dplyr::mutate_at("name", paste, region, sep = "_") %>%
      tibble::column_to_rownames("name") %>%
      as.matrix,
    transform = map(seq_table_illumina, region_illumina, region, .tag_out = seq_table, .id = c(seq_run, region))
  ),
  
  preconseq_illumina = target(
    drake_combine(region_illumina) %>%
      purrr::map(~tibble::tibble(value = as.character(.@sread), seq = as.character(.@id))) %>%
      tibble::enframe(name = "region") %>%
      dplyr::mutate_at("region", sub, pattern = "region_illumina_[[:alnum:]]+[._]\\d+_", replacement = "") %>%
      tidyr::unnest("value") %>%
      dplyr::mutate_at("value", chartr, old = "T", new = "U") %>%
      tidyr::pivot_wider(names_from = "region", values_from = "value") %>%
      dplyr::left_join(
        seq_table_illumina %>%
          t() %>%
          tibble::as_tibble(rownames = NULL) %>%
          dplyr::mutate(seq = as.character(seq.int(nrow(.)))) %>%
          tidyr::pivot_longer(cols = -"seq", names_to = "name", values_to = "nread") %>%
          dplyr::filter(nread > 0),
        by = "seq"
      ) %>%
      dplyr::select(-seq, -name) %>%
      dplyr::group_by_at(dplyr::vars(-nread)) %>%
      dplyr::summarize(nread = sum(nread)) %>%
      dplyr::mutate(primer_ID = symbols_to_values(primer_ID)) %>%
      dplyr::group_by(ITS2),
    transform = combine(region_illumina, seq_table_illumina, seq_run, primer_ID, .by = seq_run)
  ),

  # Find likely bimeras
  chimeras = target(
    if (length(seq_table) == 0) {
      logical(0)
    } else if (is.matrix(seq_table)) {
      dada2::isBimeraDenovoTable(
        seq_table,
        multithread = ignore(ncpus)
      )
    } else {
      dada2::isBimeraDenovo(
        seq_table,
        multithread = ignore(ncpus)
      )
    },
    transform = map(seq_table, .tag_in = step, .id = c(seq_run, region))
  ),

  allchimeras = target(
    names(which(c(chimeras))) %>%
      chartr("T", "U", .) %>%
      tzara::seqhash() %>%
      unname(),
    transform = combine(chimeras, .tag_in = step, .by = region)
  ),

  # big_seq_table
  # Join all the sequence tables for each region
  big_seq_table = target(
    dada2::mergeSequenceTables(tables = list(seq_table)),
    transform = combine(seq_table, .tag_in = step, .by = region)),

  # big_fasta
  # write the big_seq_table as a fasta files so that they can be clustered by
  # VSEARCH.
  big_fasta = target(
    write_big_fasta(big_seq_table,
                    file_out(big_fasta_file)),
    transform = map(.data = !!region_meta, .id = region)
  ),

  #### Plan section 3: Generate consensus sequences ####
  # Join the raw read info
  preconseq = target({
    # dada is not dynamic, so load the actual data
    dadalist <- readd_list(dada, cache = ignore(cache_dir))
    names(dadalist) <- gsub("dada_", "", names(dadalist))
    # derep2 is dynamic, so this will only load a list of hash vectors
    # This helps to keep RAM requirements down.
    dereplist <- gett_list(derep2, cache = ignore(cache_dir))
    names(dereplist) <- gsub("derep2_", "", names(dereplist))
    # lengthlist is also dynamic, but it is small so load the actual data
    lengthlist <- readd_list(derep2_length, cache = ignore(cache_dir))
    names(lengthlist) <- gsub("derep2_length_", "", names(lengthlist))
    # Make sure they are in the same order (not guaranteed by the combine transform)
    dereplist <- dereplist[names(dadalist)]
    lengthlist <- lengthlist[names(dadalist)]
    # Join all the elements together to make one list
    dadalist <- do.call(c, c(dadalist, use.names = FALSE))
    dereplist <- do.call(c, c(dereplist, use.names = FALSE))
    lengthlist <- do.call(c, c(lengthlist, use.names = FALSE))
    
    dereplist[lengthlist == 0] <- list(NULL)
    names(dereplist) <- names(dadalist)
    
    # Combine into a master list
    dada_key <- tibble::tibble(dadalist, dereplist, lengthlist, name = names(dadalist)) %>%
      tidyr::extract(
        name,
        c("seq_run", "plate", "well", "region"),
        "([pi][sb]_\\d{3})_(\\d{3})_([A-H]\\d{1,2})_(.+)",
        remove = FALSE
      ) %>%
      dplyr::filter(lengthlist > 0) %>%
      dplyr::select(-lengthlist)
    dada_regions <- unique(dada_key[["region"]])
    
    # Remove references to the original lists so the memory can be freed.
    remove(dadalist)
    remove(dereplist)
    remove(lengthlist)
    
    # process each sample seperately, removing data as it gets processed
    dadakey <- dada_key %>%
      dplyr::group_by(seq_run, plate, well) %>%
      dplyr::group_split()
    out <- tibble::tibble()
    for (i in rev(seq_along(dadakey))) {
      out <- dplyr::bind_rows(
        out,
        dplyr::mutate_at(
          dadakey[[i]],
          "dereplist",
          purrr::map,
          ignore(drake::drake_cache(cache_dir))$get_value
        ) %>%
          do.call(what = multidada)
      )# %>%
        # dplyr::group_by_at(dada_regions) %>%
        # dplyr::summarize_at("nread", sum)
      dadakey[[i]] <- NULL
      # gc()
      futile.logger::flog.info("dadakey: %d elements, size = %f", length(dadakey), object.size(dadakey))
      futile.logger::flog.info("multidada: %d rows, size = %f", nrow(out), object.size(out))
    }
    out %>%
      dplyr::group_by_at(dada_regions) %>%
      dplyr::summarize_at("nread", sum) %>%
      dplyr::filter(!is.na(ITS2)) %>%
      dplyr::mutate(primer_ID = symbols_to_values(primer_ID)) %>%
      region_concat("LSU", c("LSU1", "D1", "LSU2", "D2", "LSU3", "D3", "LSU4")) %>%
      region_concat("ITS", c("ITS1", "5_8S", "ITS2"), "ITS2") %>%
      region_concat("long", c("ITS", "LSU"), "ITS2") %>%
      region_concat("short", c("5_8S", "ITS2", "LSU1"), "ITS2") %>%
      region_concat("32S", c("5_8S", "ITS2", "LSU"), "ITS2") %>%
      region_concat("conserv", c("5_8S", "LSU1", "LSU2", "LSU3", "LSU4")) %>%
      dplyr::group_by(ITS2) %>%
      dplyr::filter(sum(nread) >= 3)
  },
  transform = combine(
    dada, derep2, derep2_length,
    primer_ID,
    .tag_in = step,
    .by = primer_ID
  ),
  memory_strategy = "unload"
  ),

  conseq = target(
    dplyr::bind_rows(preconseq, preconseq_illumina) %>%
      # for LSU1 and 5.8S, use only long reads if there are enough
      # if there are not enough, use short reads
      dplyr::mutate(
        !! region := if (region %in% c("5_8S", "LSU1")) {
        if (sum(nread[seq_run == "pb_500" && !is.na(.data[[region]])] >= 3)) {
          ifelse(seq_run == "pb_500", .data[[region]], NA_character_)
        } else {
          ifelse(seq_run == "pb_500", NA_character_, .data[[region]])
        }
      } else {
        .data[[region]]
      }  ) %>%
      dplyr::filter(sum(nread[!is.na(.data[[region]])]) >= 3) %>%
      dplyr::group_by_at(c("ITS2", !!region)) %>%
      dplyr::summarize(nread = sum(nread)) %>%
      dplyr::group_by(ITS2) %>%
      dplyr::summarize(
        !! region := as.character(
          tzara::cluster_consensus(
            seq = .data[[region]],
            nread = nread,
            names = tzara::seqhash(.data[[region]]),
            ncpus = ignore(ncpus))
          ),
        nread = sum(nread)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(
        stringr::str_count(.[[region]], "[MRWSYKVHDBN]") < 3,
        !is.na(.[[region]]),
        nchar(.[[region]]) >= min_length
      ) %>%
      as.data.frame(),
    transform = map(
      .data = !!conseq_meta,
      .id = c(primer_ID, region)
    ),
    format = "fst"
  ),

  # make a data frame of all consensus sequences and ASVs; each row is an ITS2 ASV
  allseqs = target(
    make_allseq_table(list(conseq),
                      drake_combine(big_seq_table)) %>%
      dplyr::filter(!is.na(ITS2)) %>%
      dplyr::mutate(full = dplyr::coalesce(long, short)),
    transform = combine(conseq, big_seq_table)
  ),

  #### Plan section 4: Assign taxonomy to consensus sequences ####
  # import the RDP, Silva, and UNITE databases, and format them for use with
  # DECIPHER.
  refdb_dada2 = target(
    file_in(reference_file),
    transform = map(
      .data = !!filter(ref_meta, method == "dada2"),
      .tag_out = refdb,
      .id = ref_ID)
  ),


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

  # Assign taxonomy to each ASV
  taxon = target({
    seqs <- dplyr::select(allseqs, region, "hash") %>%
      dplyr::filter(complete.cases(.)) %>%
      unique() %>%
      {set_names(.[[region]], .$hash)}
    taxonomy(seq = seqs,
             reference = refdb,
             method = method,
             multithread = ignore(ncpus)) %>%
      taxtable(names = names(seqs))
    },
    transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID)
  ),

  taxon_table = target(
    drake_combine(taxon) %>%
     combine_taxon_tables(allseqs),
    transform = combine(taxon)
  ),

  #### Plan Section 5: Build phylogenetic trees. ####
  
  # Get the LSU sequences ready to align.
  prealn_decipher_LSU = 
    allseqs %>%
    dplyr::select(hash, long, LSU, ITS1, ITS2) %>%
    dplyr::filter(
      complete.cases(.),
      startsWith(long, ITS1),
      endsWith(long, LSU),
      stringi::stri_detect_fixed(long, chartr("T", "U", ITS2)),
      !hash %in% allchimeras_ITS2
    ) %>%
    unique() %>%
    dplyr::arrange(hash) %>%
    dplyr::filter(!duplicated(long)) %>%
    dplyr::filter(!duplicated(LSU)) %$%
    set_names(LSU, hash) %>%
    chartr("T", "U", .) %>%
    Biostrings::RNAStringSet(),
  
  # Align the consensus LSU sequences
  aln_decipher_LSU = DECIPHER::AlignSeqs(
    prealn_decipher_LSU,
    iterations = 10,
    refinements = 10,
    processors = ignore(ncpus)
  ),
  
  # Trim the intron-containing regions
  aln_decipher_LSU_trim = trim_LSU_intron(aln_decipher_LSU),

  # Build a tree based on LSU
  raxml_decipher_LSU = {
    if (!dir.exists(!!config$raxml_dir))
      dir.create(!!config$raxml_dir, recursive = TRUE)
    wd <- setwd(!!config$raxml_dir)
    result <-
      aln_decipher_LSU_trim %>%
      Biostrings::DNAStringSet() %>%
      ape::as.DNAbin() %>%
      as.matrix() %>%
      ips::raxml(
        DNAbin = .,
        m = "GTRGAMMA",
        f = "a",
        N = "autoMRE_IGN",
        p = 12345,
        x = 827,
        k = TRUE,
        file = "decipher_LSU",
        exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
        threads = ignore(ncpus))
    setwd(wd)
    result
  },
  
  # Get the LSU sequences ready to align.
  prealn_decipher_long = allseqs %>%
    dplyr::select(hash, long, ITS1, ITS2) %>%
    dplyr::filter(
      complete.cases(.),
      startsWith(long, ITS1),
      stringi::stri_detect_fixed(long, ITS2),
      !hash %in% allchimeras_ITS2
    ) %>%
    unique() %>%
    dplyr::arrange(hash) %>%
    dplyr::filter(!duplicated(long)) %$%
    set_names(long, hash) %>%
    chartr("T", "U", .) %>%
    Biostrings::RNAStringSet(),
  
  # Align the consensus long reads
  aln_decipher_long =
    DECIPHER::AlignSeqs(
      prealn_decipher_long,
      iterations = 10,
      refinements = 10,
      processors = ignore(ncpus)
    ),
  
  # Trim the intron-containing region
  aln_decipher_long_trim = trim_LSU_intron(aln_decipher_long),
  
  # Build a tree based on the long reads, constrained by the LSU tree
  raxml_decipher_long = {
    if (!dir.exists(!!config$raxml_dir))
      dir.create(!!config$raxml_dir, recursive = TRUE)
    wd <- setwd(!!config$raxml_dir)
    result <-
      aln_decipher_long_trim %>%
      Biostrings::DNAStringSet() %>%
      ape::as.DNAbin() %>%
      as.matrix() %>%
      ips::raxml(
        DNAbin = .,
        m = "GTRGAMMA",
        f = "a",
        N = "autoMRE_IGN",
        p = 12345,
        x = 827,
        k = TRUE,
        backbone = raxml_decipher_LSU$bestTree,
        file = "decipher_long",
        exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
        threads = ignore(ncpus))
    setwd(wd)
    result
  },
  
  # Align the long reads to the CM
  # Don't include any positions which are in the variable regions or which
  # are insertions with respect to the CM.
  aln_infernal_32S = 
    allseqs %>%
    dplyr::select(hash, `32S`, long, ITS1, ITS2) %>%
    dplyr::filter(
      complete.cases(.),
      startsWith(long, ITS1),
      stringi::stri_detect_fixed(`32S`, ITS2),
      !hash %in% allchimeras_ITS2
    ) %>%
    unique() %>%
    dplyr::arrange(hash) %$%
    set_names(`32S`, hash) %>%
    chartr("T", "U", .) %>%
    Biostrings::RNAStringSet() %>%
    cmalign(cmfile = "reference/fungi_32S_LR5.cm") %>%
    remove_nonconsensus_nongaps(noncons_chars = c(".", "v", "p")),
  
  # Make a tree based on the conserved parts of the long reads
  raxml_infernal_32S = {
    if (!dir.exists(!!config$raxml_dir))
      dir.create(!!config$raxml_dir, recursive = TRUE)
    wd <- setwd(!!config$raxml_dir)
    result <-
      aln_infernal_32S %>%
      Biostrings::DNAStringSet() %>%
      ape::as.DNAbin() %>%
      as.matrix() %>%
      ips::raxml(
        DNAbin = .,
        m = "GTRGAMMA",
        f = "a",
        N = "autoMRE_IGN",
        p = 12345,
        x = 827,
        k = TRUE,
        file = "decipher_32S",
        exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
        threads = ignore(ncpus))
    setwd(wd)
    result
  },
  
  # Get the short and long reads ready to align with MAFFT.
  prealn_mafft_full =
    allseqs %>%
    dplyr::select(hash, full) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::arrange(dplyr::desc(nchar(full))) %>%
    dplyr::filter(!duplicated(full)) %>%
    dplyr::filter(!duplicated(hash)) %>%
    dplyr::filter(!hash %in% names(aln_decipher_long_trim)) %$%
    set_names(full, hash) %>%
    chartr("U", "T", .) %>%
    Biostrings::DNAStringSet(),
  
  # Add the short reads to the long reads alignment using MAFFT
  aln_mafft_full =
    mafft_add(
      x = Biostrings::DNAStringSet(aln_decipher_long_trim),
      y = prealn_mafft_full,
      add = "add",
      method = "10merpair",
      maxiterate = "2",
      thread = ignore(ncpus),
      quiet = FALSE,
      exec = Sys.which("mafft")
    ),
  
  # Place the short reads on the long-read tree using EPA
  epa_full = 
    epa_ng(
      ref_msa = aln_mafft_full[names(aln_mafft_full) %in% names(aln_decipher_long)],
      tree = raxml_decipher_long$bestTree,
      query = aln_mafft_full[!names(aln_mafft_full) %in% names(aln_decipher_long)],
      model = raxml_decipher_long$info,
      threads = ignore(ncpus),
      exec = Sys.which("epa-ng"),
      redo = TRUE
    ),
  
  # Graft the EPA tree
  graft_full =
    gappa_graft(
      jplace = epa_full,
      threads = ignore(ncpus),
      verbose = TRUE
    ),
  
  # Convert the grafted tree to a guide tree
  guidetree_full =
    grafts_to_polytomies(graft_full, raxml_decipher_long$bestTree),
  
  # Resolve the EPA tree
  raxml_epa_full = {
    if (!dir.exists(!!config$raxml_dir))
      dir.create(!!config$raxml_dir, recursive = TRUE)
    wd <- setwd(!!config$raxml_dir)
    result <-
      aln_mafft_full %>%
      ape::as.DNAbin() %>%
      as.matrix() %>%
      ips::raxml(
        DNAbin = .,
        m = "GTRGAMMA",
        f = "d",
        N = 1,
        p = 12345,
        backbone = guidetree_full,
        file = "epa_long",
        exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
        threads = ignore(ncpus))
    setwd(wd)
    result
  },
  
  # Build a tree based on the long alignment, unconstrained by the LSU tree
  raxml_decipher_unconst_long = {
    if (!dir.exists(!!config$raxml_dir))
      dir.create(!!config$raxml_dir, recursive = TRUE)
    wd <- setwd(!!config$raxml_dir)
    result <-
      aln_decipher_long_trim %>%
      Biostrings::DNAStringSet() %>%
      ape::as.DNAbin() %>%
      as.matrix() %>%
      ips::raxml(
        DNAbin = .,
        m = "GTRGAMMA",
        f = "a",
        N = "autoMRE_IGN",
        p = 12345,
        x = 827,
        k = TRUE,
        file = "decipher_unconst_long",
        exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
        threads = ignore(ncpus))
    setwd(wd)
    result
  },
  
  prealn_decipher_short =
    allseqs %>%
    dplyr::select(hash, short, ITS2) %>%
    dplyr::filter(
      complete.cases(.),
      stringi::stri_detect_fixed(short, chartr("T", "U", ITS2)),
      !hash %in% allchimeras_ITS2
    ) %>%
    unique() %>%
    dplyr::arrange(hash) %>%
    dplyr::filter(
      !duplicated(short),
      !hash %in% names(aln_decipher_long_trim)
    ) %$%
    set_names(short, hash) %>%
    chartr("T", "U", .) %>%
    Biostrings::RNAStringSet(),
  
  # align the short reads
  aln_decipher_short =
    DECIPHER::AlignSeqs(
      prealn_decipher_short,
      iterations = 20, # give extra iterations here
      refinements = 10,
      processors = ignore(ncpus)
    ),
  
  # merge the short read alignment with the long read alignment
  aln_decipher_full =
    DECIPHER::AlignProfiles(
      pattern = aln_decipher_long_trim,
      subject = aln_decipher_short,
      terminalGap = 0,
      processors = ignore(ncpus)
    ),
  
  # Place short reads on the unconstrained long tree using the DECIPHER alignment
  epa_decipher_unconst_full = 
    epa_ng(
      ref_msa = aln_decipher_full[names(aln_decipher_full) %in% names(aln_decipher_long_trim)],
      tree = raxml_decipher_unconst_long$bestTree,
      query = aln_decipher_full[!names(aln_decipher_full) %in% names(aln_decipher_long_trim)],
      model = raxml_decipher_unconst_long$info,
      threads = ignore(ncpus),
      exec = Sys.which("epa-ng"),
      redo = TRUE,
    ),
  
  # Graft the EPA tree
  graft_decipher_unconst_full =
    gappa_graft(
      jplace = epa_decipher_unconst_full,
      threads = ignore(ncpus),
      verbose = TRUE
    ),
  
  # Convert the grafted tree to a guide tree
  guidetree_decipher_unconst_full =
    grafts_to_polytomies(
      graft_decipher_unconst_full,
      raxml_decipher_unconst_long$bestTree
    ),
  
  # Resolve the EPA tree
  raxml_decipher_unconst_full = {
    if (!dir.exists(!!config$raxml_dir))
      dir.create(!!config$raxml_dir, recursive = TRUE)
    withr::with_dir(
      !!config$raxml_dir,
      aln_decipher_full %>%
        Biostrings::DNAStringSet() %>%
        ape::as.DNAbin() %>%
        as.matrix() %>%
        ips::raxml(
          DNAbin = .,
          m = "GTRGAMMA",
          f = "d",
          N = 1,
          p = 12345,
          backbone = guidetree_decipher_unconst_full,
          file = "epa_unconst_long",
          exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
          threads = ignore(ncpus)
        )
    )
  },
  
  epa_iterate_full = target(
    epa_iterate(
      subject = aln_decipher_full[names(aln_decipher_full) %in% names(aln_decipher_long_trim)],
      query = aln_decipher_full[!names(aln_decipher_full) %in% names(aln_decipher_long_trim)],
      subject_tree = raxml_decipher_long$bestTree,
      subject_model = raxml_decipher_long$info,
      iterations = 50,
      threads = ignore(ncpus)
    )
  ),
  
  # Graft the EPA tree
  graft_decipher_iterate_full =
    gappa_graft(
      jplace = epa_iterate_full$jplace,
      threads = ignore(ncpus),
      verbose = TRUE
    ),
  
  # Convert the grafted tree to a guide tree
  guidetree_decipher_iterate_full =
    grafts_to_polytomies(
      graft_decipher_iterate_full,
      raxml_decipher_unconst_long$bestTree
    ),
  
  # Resolve the iterated EPA tree
  raxml_decipher_iterate_full = {
    if (!dir.exists(!!config$raxml_dir))
      dir.create(!!config$raxml_dir, recursive = TRUE)
    withr::with_dir(
      !!config$raxml_dir,
      c(
        epa_iterate_full$subject,
        epa_iterate_full$query
      ) %>%
        Biostrings::DNAStringSet() %>%
        ape::as.DNAbin() %>%
        as.matrix() %>%
        ips::raxml(
          DNAbin = .,
          m = "GTRGAMMA",
          f = "d",
          N = 1,
          p = 12345,
          backbone = guidetree_decipher_unconst_full,
          file = "epa_iter_long",
          exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
          threads = ignore(ncpus)
        )
    )
  },
  
  #### Plan section 6: Gather stats ####
  
  qstats_derep2 = target(
    attr(derep2, "qstats") %>%
      mutate(region = region) %>%
      as.data.frame(),
    transform = map(derep2, region, .id = FALSE),
    dynamic = map(derep2),
    format = "fst"
  ),
  
  qstats_illumina = target(
    attr(derep_illumina, "qstats") %>% as.data.frame(),
    transform = map(derep_illumina, .id = FALSE),
    dynamic = map(derep_illumina),
    format = "fst"
  ),
  
  raw_fastq = c(
    list.files(
      file_in(!!config$fastqdir),
      pattern = ".fastq.gz",
      full.names = TRUE
    ),
    list.files(
      file_in(!!file.path(config$rawdir, datasets$dataset[datasets$tech == "Illumina"])),
      pattern = ".fastq.gz$",
      full.names = TRUE,
      recursive = TRUE
    ),
    list.files(
      file_in(!!file.path(config$rawdir, datasets$dataset[datasets$tech == "Ion Torrent"])),
      pattern = "^rawlib.basecaller.bam",
      full.names = TRUE,
      recursive = TRUE
    )
  ),
    
  qstats_raw = target(
    q_stats(raw_fastq, step = "raw") %>%
      dplyr::mutate(
        file = ifelse(basename(file) == "rawlib.basecaller.bam", "is_057-001", file)
      ) %>%
      as.data.frame(),
    dynamic = map(raw_fastq),
    format = "fst"
  ),
  
  qstats_demux = target(
    q_stats(file_meta$trim_file, step = "demux") %>% as.data.frame(),
    transform = map(file_meta, .id = FALSE),
    dynamic = map(file_meta),
    format = "fst"
  ),
  
  qstats_demux_illumina = target(
    dplyr::bind_rows(
      q_stats(illumina_group$trim_file_R1),
      q_stats(illumina_group$trim_file_R2)
    ) %>%
      as.data.frame(),
    transform = map(illumina_group, .id = FALSE),
    dynamic = map(illumina_group),
    format = "fst"
  ),
  
  qstats = target(
    combine_dynamic_diskframe(
      c(
        qstats_raw,
        qstats_demux,
        qstats_derep2,
        qstats_illumina,
        qstats_demux_illumina
      ),
      cache = ignore(cache_dir)
    ) %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      dplyr::group_by_at(dplyr::vars(-nreads)) %>%
      dplyr::summarize_at("nreads", sum) %>%
      dplyr::ungroup(),
    transform = combine(
      qstats_derep2,
      qstats_demux,
      qstats_illumina,
      qstats_demux_illumina
    )
  ),
  
  trace = TRUE
)
tictoc::toc()

cache_dir <- ".drake_heavy"
cache <- drake::drake_cache(cache_dir)
if (is.null(cache)) cache <- drake::new_cache(".drake_heavy")


tictoc::tic()
flog.info("Configuring plan...")
dconfig <- drake_config(plan, jobs_preprocess = local_cpus(),
                        cache = cache)
tictoc::toc()
flog.info("Calculating outdated targets...")
tictoc::tic()
od <- outdated(dconfig)
tictoc::toc()

if (interactive()) {
  # plansteps <- unique(na.omit(plan[["step"]]))
  vis_drake_graph(dconfig,
                  # group = "step", clusters = plansteps,
                  targets_only = TRUE)
  ncpus <- local_cpus()
  make(plan, target = "raxml_epa_full", max_expand = config$smallsplit,
       cache = cache)
}

remove(snakemake)
remove(dconfig)
remove(cache)

save(list = ls(),
     file = config$drakedata)
