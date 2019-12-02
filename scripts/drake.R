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
  # filter(tech != "Ion") %>%
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
  
  position_trace = target(
    get_trace(paste0("derep_by_", seq_run), positions),
    transform = map(positions, seq_run, .id = seq_run)
  ),
  
  position_map = target(
    dplyr::select(positions, "seq_run", "plate", "well") %>%
      unique(),
    transform = map(positions, .id = seq_run),
    dynamic = map(positions)
  ),
  
  # derep2----
  # extract regions, do quality filtering, and dereplicate
  # this is all one step because the extraction and filtering are fast,
  # but saving the outputs would be large.
  derep2 = target(
    extract_and_derep(
      positions = positions,
      trim_file = unique(positions$trim_file),
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

  # err----
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

  # dada----
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

  # seq_table----
  # Make a sample x ASV abundance matrix
  # (for each size range, if multiple)
  seq_table = target(
    dada2::makeSequenceTable(purrr::compact(dada)) %>%
      magrittr::extract(,nchar(colnames(.)) >= min_length_post &
                          nchar(colnames(.)) <= max_length_post,
                        drop = FALSE),
    transform = map(dada, .tag_in = step, .id = c(seq_run, region))
  ),

  # nochim----
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
    transform = map(seq_table, .id = c(seq_run, region))),

  # big_seq_table ----
  # Join all the sequence tables for each region
  big_seq_table = target(
    dada2::mergeSequenceTables(tables = list(seq_table)),
    transform = combine(seq_table, .tag_in = step, .by = region)),

  # big_fasta ----
  # write the big_seq_table as a fasta files so that they can be clustered by
  # VSEARCH.
  big_fasta = target(
    write_big_fasta(big_seq_table,
                    file_out(big_fasta_file)),
    transform = map(.data = !!region_meta, .id = region)
  ),

  # raw ----
  # Join the raw read info
  preconseq = target({
    dadalist <- drake_combine(dada)
    names(dadalist) <- gsub("dada_", "", names(dadalist))
    dereplist <- readd_list(derep2)
    names(dereplist) <- gsub("derep2_", "", names(dereplist))
    dereplist <- dereplist[names(dadalist)]
    dadalist <- do.call(c, c(dadalist, use.names = FALSE))
    dereplist <- do.call(c, c(dereplist, use.names = FALSE))
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

  # allseqs ----
  # make a data frame of all consensus sequences and ASVs; each row is an ITS2 ASV
  allseqs = target(
    make_allseq_table(list(conseq),
                      drake_combine(big_seq_table)) %>%
      dplyr::filter(!is.na(ITS2)) %>%
      dplyr::mutate(full = dplyr::coalesce(long, short)),
    transform = combine(conseq, big_seq_table),
    format = "fst"
  ),

  # dbprep ----
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

  # taxon ----
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
    transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID),
    format = "fst"),

  taxon_table = target(
    drake_combine(taxon) %>%
     combine_taxon_tables(allseqs),
    transform = combine(taxon)
  ),

  taxon_labels = target(
    make_taxon_labels(taxon_table)
  ),

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
    read_platemap(file_in(!!config$platemap), !!config$platemap_sheet),
    format = "fst"),

  # long_consensus----
  # get the long amplicon consensus and convert to RNAStringSet.
  # Use the sequence hashes as the name, so that names will be robust
  # make sure that we only use sequences that actually match the long sequence.
  cons_32S = allseqs %>%
    dplyr::select("32S", ITS1, long, hash) %>%
    dplyr::filter(
      complete.cases(.),
      endsWith(long, `32S`)
    ) %>%
    dplyr::select("32S", hash) %>%
    unique() %>%
    assertr::assert(assertr::is_uniq, hash) %$%
    rlang::set_names(`32S`, hash) %>%
    chartr(old = "T", new = "U") %>%
    Biostrings::RNAStringSet(),
  cons_LSU = allseqs %>%
    dplyr::select(LSU, long, ITS1, hash) %>%
    dplyr::filter(
      complete.cases(.),
      endsWith(long, LSU)
    ) %>%
    dplyr::select(LSU, hash) %>%
    unique() %$%
    rlang::set_names(LSU, hash) %>%
    chartr(old = "T", new = "U") %>%
    Biostrings::RNAStringSet(),
  cons_long = allseqs %>%
    dplyr::select(long, ITS1, hash) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::select(long, hash) %>%
    unique() %$%
    rlang::set_names(long, hash) %>%
    chartr(old = "T", new = "U") %>%
    Biostrings::RNAStringSet(),

  # cmaln_32S
  # Align conserved positions and annotate conserved secondary structure of the
  # 32S consensus using Infernal.
  cmaln_32S =
    cmalign(cmfile = file_in(!!config$cm_32S),
            seq = cons_32S,
            glocal = TRUE,
            cpu = ignore(ncpus)),

  # cmaln_long
  # Concatenate ITS1 (padded with gaps) to the beginning of the Infernal 32S
  # alignment to make a preliminary long amplicon alignment with conserved
  # 5.8S/LSU secondary structure annotations.
  cmaln_long =
    dplyr::inner_join(
      allseqs %>%
        dplyr::select(hash, long, ITS1) %>%
        dplyr::filter(complete.cases(.), startsWith(long, ITS1)) %>%
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
            max(nchar(.$ITS1)),
            "right",
            "-"),
          cmaln_32S$SS_cons
        ),
        # Don't take the alignment in the variable regions as
        # conserved.
        ref = stringr::str_pad(chartr("v", ".", cmaln_32S$RF),
                               max(nchar(.$aln)), "left", "-"),
        seq_names = .$hash,
        file = file_out(!!config$cmaln_long)
      )
    },

  # guidetree_32S ----
  # Take only the conserved, gap-free positions in the 32S alignment
  aln_32S_conserv =
    remove_nonconsensus_nongaps(cmaln_32S),

  # Generate a distance matrix from the conserved 32S alignment.
  dist_32S_conserv =
    DECIPHER::DistanceMatrix(
      aln_32S_conserv,
      type = "dist",
      processors = ignore(ncpus)
    ),

  # Make a UPGMA guide tree for mlocarna based on the aligned positions in
  # the Infernal alignment for 32S.
  guidetree_32S =
    DECIPHER::IdClusters(
      myDistMatrix = dist_32S_conserv,
      collapse = -1, #break polytomies
      type = "dendrogram",
      processors = ignore(ncpus)
    ) %T>%
    DECIPHER::WriteDendrogram(
      file = file_out(!!config$guide_tree)
    ),

  # create the pair probability files for mlocarna
  # this only needs to be redone if the sequences change, and it can be
  # efficiently done in parallel
  # It is made in a mirror directory, because mlocarna will write a file to the
  # same directory.
  mlocarna_pp_long = {
    mlocarna_realign(
      alignment = file_in(!!config$cmaln_long),
      target_dir = file_out(!!config$mlocarna_pp_dir),
      cpus = ignore(ncpus),
      only_dps = TRUE
    )
    mirror_dir(!!config$mlocarna_pp_dir, !!config$mlocarna_dir)
  },

  # realign the consensus sequences using mlocarna
  # this is a progressive alignment, where each locarna alignment is single
  # threaded, and the ones closer to the root of the tree take a lot longer.
  # it is efficient to do this in parallel near the beginning, but only one or
  # a few cores can be utilized at the end; mlocarna itself does not execute
  # them in parallel at all.
  # snakemakelocarna.sh is a wrapper for locarna which writes empty files for
  # every locarna call (thus spoofing mlocarna into thinking that it is working)
  # while collecting the dependency structure of the calls.  On the last call, it
  # submits all the calls as separate single-threaded jobs on SLURM using
  # Snakemake.
  # This requires a lot of configuration to work properly on a given cluster
  # system; see files in makelocarna_profile.
  # (default "config/snakemakelocarnaUPPMAX")
  # This will always fail on the first run, because it does not wait for all of
  # the locarna calls to complete.  Check the cluster queue and make sure
  # all of the jobs have completed, then this rule should pass on the next run
  # of drake.
  realign_long = {
    mlocarna_pp_long
    file_out(!!config$mlocarna_aln)
    mlocarna_realign(
      alignment = file_in(!!config$cmaln_long),
      guide_tree = file_in(!!config$guide_tree),
      target_dir = file_out(!!config$mlocarna_dir),
      stockholm = TRUE,
      consensus_structure = "alifold",
      cpus = 1,
      skip_pp = TRUE,
      pw_aligner = normalizePath(file_in(!!config$makelocarna)),
      pw_aligner_options = paste(
        "--profile", normalizePath(!!config$makelocarna_profile),
        "--conda", normalizePath(!!config$makelocarna_conda)
      )
    )
  },

  # read the mlocarna output
  aln_locarna_long =
    read_stockholm_msa(file_in(!!config$mlocarna_aln)),

  # make a tree based on the realigned consensus using RAxML
  raxml_locarna_long = {
      raxml_RNA(
        RNAaln = aln_locarna_long$alignment,
        S = aln_locarna_long$SS_cons,
        m = "GTRGAMMA",
        A = "S16",
        f = "a",
        N = "autoMRE_IGN",
        p = 12345,
        x = 827,
        k = TRUE,
        file = "locarna_long",
        dir = !!config$raxml_locarna_out_dir,
        exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
        threads = ignore(ncpus)
      )
    setwd(wd)
    result
  },

  aln_decipher_long =
    allseqs %>%
    dplyr::select(hash, long, ITS1, ITS2) %>%
    dplyr::filter(
      complete.cases(.),
      startsWith(long, ITS1),
      stringi::stri_detect_fixed(long, ITS2)
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
  
  aln_decipher_long_trim = trim_LSU_intron(aln_decipher_long),

  aln_decipher_LSU =
    allseqs %>%
    dplyr::select(hash, long, LSU, ITS1, ITS2) %>%
    dplyr::filter(
      complete.cases(.),
      startsWith(long, ITS1),
      endsWith(long, LSU),
      stringi::stri_detect_fixed(long, chartr("T", "U", ITS2))
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

  aln_decipher_LSU_trim = trim_LSU_intron(aln_decipher_LSU),

  raxml_decipher_LSU = {
    if (!dir.exists(!!config$raxml_decipher_dir))
      dir.create(!!config$raxml_decipher_dir, recursive = TRUE)
    wd <- setwd(!!config$raxml_decipher_dir)
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

  raxml_decipher_long = {
    if (!dir.exists(!!config$raxml_decipher_dir))
      dir.create(!!config$raxml_decipher_dir, recursive = TRUE)
    wd <- setwd(!!config$raxml_decipher_dir)
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
  
  aln_decipher_full =
    allseqs %>%
    dplyr::select(hash, full) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::arrange(dplyr::desc(nchar(full))) %>%
    dplyr::filter(!duplicated(full)) %>%
    dplyr::filter(!duplicated(hash)) %$%
    set_names(full, hash) %>%
    chartr("T", "U", .) %>%
    Biostrings::RNAStringSet() %>%
    DECIPHER::AlignSeqs(iterations = 10,
                        refinements = 10,
                        processors = ignore(ncpus)),
  
  raxml_decipher_full = {
    if (!dir.exists(!!config$raxml_decipher_dir))
      dir.create(!!config$raxml_decipher_dir, recursive = TRUE)
    wd <- setwd(!!config$raxml_decipher_dir)
    result <-
      aln_decipher_full %>%
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
        file = "decipher_full",
        exec = Sys.which("raxmlHPC-PTHREADS-SSE3"),
        threads = ignore(ncpus))
    setwd(wd)
    result
  },
    
    

  labeled_tree_decipher_long =
    relabel_tree(
      tree = raxml_decipher_long$bipartitions,
      old = taxon_labels$label,
      new = taxon_labels$tip_label
    ) %T>%
    castor::write_tree(file_out("data/trees/decipher_long.tree")),

  labeled_tree_mlocarna_long =
    relabel_tree(
      tree = raxml_mlocarna_long$bipartitions,
      old = taxon_labels$label,
      new = taxon_labels$tip_label
    ) %T>%
    castor::write_tree(file_out("data/trees/mlocarna_long.tree")),


  # Tulasnella is on a long branch because of a high divergence rate.
  # it ends up with Metazoa in the tree.  Exclude it.
  tulasnella =
    taxon_table %>%
    dplyr::filter(
      rank == "family",
      taxon == "Tulasnellaceae"
    ) %$%
    label %>%
    unique(),

  # find the sequences that were confidently identified as fungi
  # (and placed in a phylum).  We extract the common ancestor of these
  surefungi =
    taxon_table %>%
    dplyr::group_by(label) %>%
    dplyr::filter(
      rank == "phylum",
      n_tot >= 6,
      n_diff == 1,
      grepl("mycota", taxon)
    ) %$%
    label %>%
    unique() %>%
    setdiff(tulasnella),

  # find the most abundant sequence which we are confident is a plant
  bestplant =
    taxon_table %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank == "phylum",
                  taxon == "Chlorophyta",
                  n_tot == 6,
                  n_diff == 1) %>%
    dplyr::filter(n_reads == max(n_reads)) %$%
    label[1],

  # root the tree within plants.  This isn't accurate, but it should guarantee
  # that it is not rooted within fungi.
  tree_decipher_long =
    ape::root(raxml_decipher_long$bipartitions,
              bestplant),

  # Extract the minimally inclusive clade including all confidently identified
  # fungi (except Tulasnella).
  fungi_tree_decipher_long =
    ape::getMRCA(
      phy = tree_decipher_long,
      tip = intersect(surefungi, tree_decipher_long$tip.label)
    ) %>%
    ape::extract.clade(phy = tree_decipher_long),

  taxon_phy = phylotax(
    tree = fungi_tree_decipher_long,
    taxa = taxon_table
  ),

  physeq = assemble_physeq(
    platemap,
    datasets,
    relabel_seqtable(big_seq_table_ITS2),
    fungi_tree_decipher_long),

  longphyseq = physeq %>%
    phyloseq::prune_samples(
      samples = phyloseq::sample_data(.)$dataset == "long-pacbio"
      & phyloseq::sample_data(.)$sample_type == "Sample"
      & (phyloseq::sample_data(.)$qual == "X" | phyloseq::sample_data(.)$year == "2015")
      & rowSums(phyloseq::otu_table(.)) > 0
      ),

  shortphyseq = physeq %>%
    phyloseq::prune_samples(
      samples = phyloseq::sample_data(.)$dataset == "short-pacbio"
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
trace = TRUE)
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
#     if (!dir.exists(!!config$outdir)) dir.create(!!config$outdir, recursive = TRUE)
#     rmarkdown::render(
#       knitr_in(!!file.path(config$rmddir, "qual-check.Rmd")),
#       output_file = file_out(!!file.path(config$outdir, "qual-check.pdf")),
#       output_dir = !!config$outdir)},
#   # max_expand = if (interactive()) 9 else NULL,
#   trace = TRUE
# )
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
  make(plan, target = "raxml_decipher_long", max_expand = config$smallsplit)
}

remove(snakemake)

save(list = ls(),
     file = config$drakedata)
