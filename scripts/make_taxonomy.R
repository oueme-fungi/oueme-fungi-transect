library(futile.logger)
if (file.exists("ENTREZ_KEY")) {
  Sys.setenv(ENTREZ_KEY = ENTREZ_KEY)
}

# define or input file names and parameters
if (interactive()) {
  flog.info("Creating drake plan in interactive session...")
  library(here)
  r_dir <- "scripts"
  config_dir <- "config"
  ref_dir <- "reference"
  tedersoo_file <- file.path(ref_dir, "Tedersoo_Eukarya_classification.xlsx")
  regions_file <- file.path(config_dir, "regions.csv")
} else if (exists("snakemake")) {
  flog.info("Creating drake plan in snakemake session...")
  snakemake@source(".Rprofile", echo = FALSE)
  r_dir <- snakemake@config$rdir
  ref_dir <- snakemake@config$ref_root
  tedersoo_file <- snakemake@input$tedersoo
  regions_file <- snakemake@input$regions
  outputs <- snakemake@output
  logfile <- file(snakemake@log[[1]], open = "at")
  sink(logfile)
  sink(logfile, type = "message")
} else {
  flog.error("Can't find Snakemake object in non-interactive session!")
  stop()
}

library(magrittr)
library(drake)

source(file.path(r_dir, "taxonomy.R"))
source(file.path(r_dir, "parallel_helpers.R"))

setup_log("translate_references")

tedersoo_ranks <- c("kingdom", "subkingdom", "phylum", "subphylum", "class",
                    "order", "family", "genus")

dbs <- tibble::tribble(        ~db, ~header_format,
                       "rdp_train",          "rdp",
                          "warcup",          "rdp",
                           "unite",        "unite") %>%
  dplyr::mutate(header_file = file.path("reference", paste0(db, ".fasta.gz")),
                patch_file = file.path("reference", paste0(db, ".pre.sed")))

regions = readr::read_csv(regions_file, col_types = "cccciiiiic") %>%
  tidyr::separate_rows(reference, sep = ", ?")
  

db_out <- dplyr::select(regions, reference) %>%
  unique() %>%
  tidyr::separate(col = reference, into = c("db", "region"), sep = "\\.") %>%
  tidyr::crossing(method = c("dada2", "sintax")) %>%
  dplyr::mutate(new_header = paste0("new_header_", db) %>%
                  rlang::syms(),
                fasta_in = file.path(ref_dir, paste( db, region, "fasta.gz",
                                                     sep = ".")),
                fasta_out = file.path(ref_dir, paste(db, region, method,
                                                     "fasta.gz", sep = "."))) %>%
  dplyr::left_join(dbs, by = "db")

plan <- drake_plan(
  # read headers from the reference database fasta files
  raw_header = target(
    read_header(header_file = file_in(header_file),
                patch_file = file_in(patch_file),
                format = header_format),
    transform = map(.data = !!dbs,
                    .id = db)
    ),
  
  # the rdp training set does not have fully annotated taxonomy for nonfungi,
  # so these need to be looked up
  # find the accession numbers representing nonfungi in the RDP database
  rdp_nf_accno =
    dplyr::filter(raw_header_rdp_train, !grepl("Fungi", classifications)) %$%
    accno %>%
    unique(),
  
  rdp_nf_taxdata = 
    target(
      taxa::lookup_tax_data(rdp_nf_accno, type = "seq_id"),
      transform = split(rdp_nf_accno, slices = 50),
      retries = 1
    ),
  
  rdp_nf_ncbiheader =
    target(
      purrr::map_dfr(list(rdp_nf_taxdata),
                     accno_c12n_table),
      transform = combine(rdp_nf_taxdata)
    ),
  
  rdp_nf_taxa =
    target(
      list(rdp_nf_taxdata) %>%
        purrr::map_dfr(~ .$data$tax_data) %>%
        unique(),
      transform = combine(rdp_nf_taxdata)
    ),
  
  tedersoo_class =
    read_classification_tedersoo(file_in(tedersoo_file)),
  
  rdp_nf_reduced =
    target(
      reduce_ncbi_taxonomy(
        rdp_nf_ncbiheader,
        rdp_nf_taxa,
        ranks = c("kingdom", "phylum", "class", "order",
                  "family", "genus"),
        keytaxa = unique(tedersoo_class$taxon_names)
      )
    ),
  
  reduced_header_rdp_train = target(
    raw_header_rdp_train %>%
      dplyr::filter(!duplicated(accno)) %>%
    dplyr::left_join(dplyr::select(rdp_nf_reduced,
                                   accno,
                                   c_new = classifications),
                     by = "accno") %>%
      dplyr::mutate(
        classifications =
          dplyr::coalesce(c_new, classifications) %>%
          reduce_taxonomy()
      ) %>%
      dplyr::select(-c_new),
    transform = map(db = "rdp_train",
                    .tag_in = !!rlang::sym("reduced_header"),
                    .id = FALSE)
  ),
  
  reduced_header = target(
    raw_header %>%
      dplyr::filter(!duplicated(accno)) %>%
      dplyr::mutate_at("classifications", reduce_taxonomy),
    transform = map(raw_header = !!rlang::syms(c("raw_header_warcup", "raw_header_unite")), db = c("warcup", "unite"), .id = db)
  ),
  
  class = target(
    reduced_header$classifications %>%
      unique() %>%
      taxa::parse_tax_data() %>%
      taxa::get_data_frame(c("taxon_names", "classifications")),
    transform = map(reduced_header, .id = db)
  ),
  
  new_header = target(
    translate_taxonomy(reduced_header,
                       class,
                       tedersoo_class) %>%
    dplyr::mutate_at("classifications", regularize_taxonomy,
                     rank_in = tedersoo_ranks),
    transform = map(class, .id = db)
  ),
  
  write = target(
    replace_header(file_in(fasta_in), file_out(fasta_out), new_header,
                   in_format = header_format, out_format = method,
                   patch_file = file_in(patch_file)),
    transform = map(.data = !!db_out,
                      .id = c(db, region, method))
  ),
  trace = TRUE
)

outfiles <- as.character(plan$command) %>%
  stringr::str_match_all('file_out\\("([^"]+)"\\)') %>%
  purrr::keep(~length(.) > 0) %>%
  purrr::map(2) %>%
  unlist() %>%
  unique()

if (exists("outputs")) {
  assertthat::assert_that(all(outputs %in% outfiles))
}

make(plan, cache = drake_cache(".drake_taxonomy"))

# make sure to touch al the outputs
if (exists("outputs")) {
  for (o in outputs) Sys.setFileTime(o, Sys.time())
}