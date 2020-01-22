library(futile.logger)
library(magrittr)
library(tidyverse)
library(rlang)
library(glue)
library(drake)
library(assertr)
library(disk.frame)

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
source(file.path(config$rdir, "extract_regions.R"))
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
datasets <- read_csv(config$dataset, col_types = "cccicccccicc")
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
flog.info("Making itsx_meta.")
itsx_meta <- datasets %>%
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

flog.info("Making err_meta.")
err_meta <- dada_meta %>%
  filter(region == err_region) %>%
  select(seq_run, region, derep2, tech, hgp, band_size, pool)

flog.info("Making conseq_meta.")
conseq_meta <- dada_meta %>%
  filter(region %in% c("ITS1", "ITS", "LSU", "32S", "long", "short")) %>%
  select(region, primer_ID, preconseq) %>%
  left_join(regions %>% select(region, min_length), by = "region") %>%
  unique()

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
      readd(derep1),
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
    combine_dynamic_diskframe(lsux),
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
                       max_length, min_length, max_ee, seq_run, region),
      .tag_in = step,
      .id = c(seq_run, region)
    ),
    dynamic = map(positions)
  ),
  
  #### Plan section 2: Denoise using DADA ####

  # Calibrate dada error models on 5.8S
  err = target({
    err.fun <- if (tech == "PacBio" ) dada2::PacBioErrfun else
      dada2::loessErrfun
    dereps <- readd(derep2) %>%
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

  # Run dada denoising algorithm (on different regions)
  dada = target(
    readd(derep2) %>%
      set_names(
        readd(position_map) %>%
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

  # Make a sample x ASV abundance matrix
  # (for each size range, if multiple)
  seq_table = target(
    dada2::makeSequenceTable(purrr::compact(dada)) %>%
      magrittr::extract(,nchar(colnames(.)) >= min_length_post &
                          nchar(colnames(.)) <= max_length_post,
                        drop = FALSE),
    transform = map(dada, .tag_in = step, .id = c(seq_run, region))
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
    dadalist <- drake_combine(dada)
    names(dadalist) <- gsub("dada_", "", names(dadalist))
    dereplist <- readd_list(derep2)
    names(dereplist) <- gsub("derep2_", "", names(dereplist))
    dereplist <- dereplist[names(dadalist)]
    dadalist <- do.call(c, c(dadalist, use.names = FALSE))
    dereplist <- do.call(c, c(dereplist, use.names = FALSE))
    dereplist[vapply(dereplist, length, 1L) == 0] <- list(NULL)
    names(dereplist) <- names(dadalist)
    dada_map <- tzara::dadamap(dereplist, dadalist) %>%
      tidyr::extract(
        name,
        c("seq_run", "plate", "well", "region"),
        "([pi][sb]_\\d{3})_(\\d{3})_([A-H]\\d{1,2})_(.+)"
      )
    regions <- unique(dada_map[["region"]])
    dada_map %>%
      dplyr::select(seq.id, seq_run, plate, well, region, dada.seq) %>%
      dplyr::mutate_at("dada.seq", chartr, old = "T", new = "U") %>%
      tidyr::spread(key = "region", value = "dada.seq") %>%
      dplyr::group_by(ITS2) %>%
      dplyr::filter(!is.na(ITS2), dplyr::n() >= 3) %>%
      dplyr::group_by_at(regions) %>%
      dplyr::summarize(nread = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(primer_ID = symbols_to_values(primer_ID)) %>%
      region_concat("LSU", c("LSU1", "D1", "LSU2", "D2", "LSU3", "D3", "LSU4")) %>%
      region_concat("ITS", c("ITS1", "5_8S", "ITS2"), "ITS2") %>%
      region_concat("long", c("ITS", "LSU"), "ITS2") %>%
      region_concat("short", c("5_8S", "ITS2", "LSU1"), "ITS2") %>%
      region_concat("32S", c("5_8S", "ITS2", "LSU"), "ITS2") %>%
      dplyr::group_by(ITS2)
    },
    transform = combine(
      dada, derep2,
       primer_ID,
      .tag_in = step,
      .by = primer_ID
    )
  ),

  conseq = target(
    preconseq %>%
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
  
  # Align the consensus LSU sequences
  aln_decipher_LSU =
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
    Biostrings::RNAStringSet() %>%
    DECIPHER::AlignSeqs(iterations = 10,
                        refinements = 10,
                        processors = ignore(ncpus)),

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

  # Align the consensus long reads
  aln_decipher_long =
    allseqs %>%
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
    Biostrings::RNAStringSet() %>%
    DECIPHER::AlignSeqs(iterations = 10,
                        refinements = 10,
                        processors = ignore(ncpus)),
  
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
  
  # Add the short reads to the long reads alignment using MAFFT
  aln_mafft_full =
    allseqs %>%
    dplyr::select(hash, full) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::arrange(dplyr::desc(nchar(full))) %>%
    dplyr::filter(!duplicated(full)) %>%
    dplyr::filter(!duplicated(hash)) %>%
    dplyr::filter(!hash %in% names(aln_decipher_long_trim)) %$%
    set_names(full, hash) %>%
    chartr("U", "T", .) %>%
    Biostrings::DNAStringSet() %>%
    mafft_add(
      x = Biostrings::DNAStringSet(aln_decipher_long_trim),
      y = .,
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
  
  # align the short reads
  aln_decipher_short =
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
    Biostrings::RNAStringSet() %>%
    DECIPHER::AlignSeqs(iterations = 20, # give extra iterations here
                        refinements = 10,
                        processors = ignore(ncpus)),
  
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
    attr(derep2, "qstats") %>% as.data.frame(),
    transform = map(derep2, .id = FALSE),
    dynamic = map(derep2),
    format = "fst"
  ),
  
  raw_fastq = list.files(
    file_in(!!config$fastqdir),
    pattern = ".fastq.gz",
    full.names = TRUE
  ),
  
  qstats_raw = target(
    q_stats(raw_fastq) %>% as.data.frame(),
    dynamic = map(raw_fastq),
    format = "fst"
  ),
  
  qstats_demux = target(
    q_stats(file_meta$trim_file) %>% as.data.frame(),
    transform = map(file_meta, .id = FALSE),
    dynamic = map(file_meta),
    format = "fst"
  ),
  
  qstats = target(
    combine_dynamic_diskframe(c(qstats_raw, qstats_demux, qstats_derep2)) %>%
      as.data.frame() %>%
      tibble::as_tibble(),
    transform = combine(qstats_derep2, qstats_demux)
  ),
  
  qstats_n = qstats %>%
    dplyr::group_by(file, step) %>%
    dplyr::summarize(nreads = dplyr::n()) %>%
    dplyr::ungroup(),
  
  qstats_length = qstats %>%
    dplyr::group_by(file, step, length) %>%
    dplyr::summarize(nreads = dplyr::n()) %>%
    dplyr::ungroup(),
  
  qstats_minq = qstats %>%
    dplyr::group_by(file, step, minq) %>%
    dplyr::summarize(nreads = dplyr::n()) %>%
    dplyr::ungroup(),
  
  # for floating points, we need to make bins so that there aren't
  # too many unique values.  Do transforms to preserve relevant distinctions
  # first
  qstats_eexp = qstats %>%
    dplyr::mutate(eexp = round(log(eexp), 2)) %>%
    dplyr::group_by(file, step, eexp) %>%
    dplyr::summarize(nreads = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(eexp = exp(eexp)),
  
  qstats_erate = qstats %>%
    dplyr::mutate(erate = round(log(erate) - log1p(-erate), 2)) %>%
    dplyr::group_by(file, step, erate) %>%
    dplyr::summarize(nreads = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(erate = exp(erate)/(1 + exp(erate))),
  
  qstats_pnoerr = qstats %>%
    dplyr::mutate(p.noerr = round(log(p.noerr) - log1p(-p.noerr), 2)) %>%
    dplyr::group_by(file, step, p.noerr) %>%
    dplyr::summarize(nreads = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p.noerr = exp(p.noerr)/(1 + exp(p.noerr))),
  
  trace = TRUE
)
tictoc::toc()

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
  ncpus <- local_cpus()
  make(plan, target = "raxml_epa_full", max_expand = config$smallsplit)
}

remove(snakemake)
remove(dconfig)

save(list = ls(),
     file = config$drakedata)
