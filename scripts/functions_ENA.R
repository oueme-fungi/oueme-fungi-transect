# functions to generate European Nucleotide Archive  submission files
# author Brendan Furneaux

compile_soil_samples <- function(sample_data) {
    as(sample_data, "data.frame") %>%
    dplyr::filter(sample_type == "Sample", tech == "PacBio") %>%
    dplyr::left_join(
      readRDS("config/tags/all.rds") %>%
        tidyr::pivot_longer(
          cols = matches("(primer|barcode|adapter)_(seq|tag)_(fwd|rev)"),
          names_to = c(".value", "direction"),
          names_pattern = "(.+)_(fwd|rev)"
        ) %>%
        tidyr::unite("seq", ends_with("seq"), sep = "") %>%
        mutate_at(vars(ends_with("tag")), na_if, "unnamed") %>%
        mutate(tag = dplyr::coalesce(barcode_tag, primer_tag)) %>%
        dplyr::select(amplicon, well, direction, tag, seq) %>%
        filter(nchar(seq) > 0) %>%
        mutate_at("tag", str_replace, "its", "ITS") %>%
        mutate_at("tag", str_replace, "lr", "LR") %>%
        tidyr::unite("tag", tag, seq, sep = ": ") %>%
        group_by(amplicon, well, direction) %>%
        summarize(tag = paste(unique(tag), collapse = ", ")) %>%
        group_by(amplicon, well) %>%
        summarize(tag = paste(unique(tag), collapse = "; ")),
      by = c("amplicon", "well"))%>%
    dplyr::mutate(
      sample_alias = glue::glue("OT-{substr(year, 3, 4)}-{site}-{x}-{substr(buffer, 1, 1)}-{substr(amplicon, 1, 1)}"),
      tax_id = "410658",
      scientific_name = "soil metagenome",
      common_name = "soil metagenome",
      sample_title = sample_alias,
      sample_description = glue::glue("Ouémé Supérieur soil transect {year}, {ifelse(site == 'Ang', 'Angaradebou', 'Gando')} sample #{x}, preserved with {buffer}, {amplicon} amplicon"),
      "project name" = "ena-STUDY-UPPSALA UNIVERISTY-19-03-2020-10:28:50:771-6",
      "experimental factor" = paste(x, "m from transect origin"),
      "sample volume or weight for dna extraction" = 0.05,
      "nucleic acid extraction" = "Zymo Research Xpedition Soil Microbe DNA MiniPrep D6202",
      "target gene" = ifelse(amplicon == "Long", "ITS1, 5.8S, ITS2, LSU", "ITS2"),
      "pcr primers" = tag,
      "pcr conditions" = ifelse(
        amplicon == "Long",
        "initial denaturation:95degC_10min; denaturation:95degC_45s; annealing:59degC_45s; extension:72degC_90s; cycles:30; final extension:72degC_10min",
        "initial denaturation:95degC_10min; denaturation:95degC_60s; annealing:56degC_45s; extension:72degC_50s; cycles:35; final extension:72degC_3min"
      ),
      "sequencing method" = ifelse(
        amplicon == "Long",
        "PACBIO_SMRT",
        "ION_TORRENT, PACBIO_SMRT, ILLUMINA"
      ),
      "investigation type" = "metagenome",
      "collection date" = case_when(
        year == "2015" ~ "2015-05",
        year == "2016" ~ "2016-06"
      ),
      "geographic location (country and/or sea)" = "Benin",
      "geographic location (latitude)" = case_when(
        site == "Ang" ~ "N 9.75456",
        site == "Gan" ~ "N 9.75678"
      ),
      "geographic location (longitude)" = case_when(
        site == "Ang" ~ "W 2.14064",
        site == "Gan" ~ "W 2.31058"
      ),
      "geographic location (region and locality)" = case_when(
        site == "Ang" ~ "Borgou department, N'Dali Commune, Forêt Classée de l'Ouémé Supérieur near Angaradebou village",
        site == "Gan" ~ "Borgou department, N'Dali Commune, Forêt Classée de l'Ouémé Supérieur near Gando village"
      ),
      "soil environmental package" = "soil",
      "geographic location (depth)" = "0-0.05",
      "environment (biome)" = "tropical savannah ENVO:01000188",
      "environment (feature)" = "ectomycorrhizal woodland",
      "environment (material)" = "tropical soil ENVO:00005778",
      "geographic location (elevation)" = "320",
      "amount or size of sample collected" = "0.125 L",
      "sample weight for dna extraction" = "0.25",
      "storage conditions (fresh/frozen/other)" = case_when(
        buffer == "Xpedition" ~ "field lysis by bead beating in Xpedition lysis solution (Zymo Research)",
        buffer == "LifeGuard" ~ "stored for approx. two months in LifeGuard Soil Preservation Solution (Quiagen)"
      )
    )
}

compile_control_samples <- function(sample_data) {
    as(sample_data, "data.frame") %>%
    dplyr::filter(sample_type == "Pos", tech == "PacBio") %>%
    dplyr::left_join(
      readRDS("config/tags/all.rds") %>%
        tidyr::pivot_longer(
          cols = matches("(primer|barcode|adapter)_(seq|tag)_(fwd|rev)"),
          names_to = c(".value", "direction"),
          names_pattern = "(.+)_(fwd|rev)"
        ) %>%
        tidyr::unite("seq", ends_with("seq"), sep = "") %>%
        mutate_at(vars(ends_with("tag")), na_if, "unnamed") %>%
        mutate(tag = dplyr::coalesce(barcode_tag, primer_tag)) %>%
        dplyr::select(amplicon, well, direction, tag, seq) %>%
        filter(nchar(seq) > 0) %>%
        mutate_at("tag", str_replace, "its", "ITS") %>%
        mutate_at("tag", str_replace, "lr", "LR") %>%
        tidyr::unite("tag", tag, seq, sep = ": ") %>%
        group_by(amplicon, well, direction) %>%
        summarize(tag = paste(unique(tag), collapse = ", ")) %>%
        group_by(amplicon, well) %>%
        summarize(tag = paste(unique(tag), collapse = "; ")),
      by = c("amplicon", "well")) %>%
    dplyr::mutate(
      sample_alias = glue::glue("OT-control-{plate}-{substr(amplicon, 1, 1)}"),
      tax_id = "5341",
      scientific_name = "Agaricus bisporus",
      common_name = "common button mushroom",
      sample_title = sample_alias,
      sample_description = glue::glue("commercially purchased button mushroom as positive control, {amplicon} amplicon"),
      isolation_source = "purchased at local grocery store",
      "geographic location (country and/or sea)" = "Sweden",
      "geographic location (region and locality)" = "Uppsala",
      "project name" = "ena-STUDY-UPPSALA UNIVERISTY-19-03-2020-10:28:50:771-6",
      "target gene" = ifelse(amplicon == "Long", "ITS1, 5.8S, ITS2, LSU", "ITS2"),
      "pcr primers" = tag,
      "pcr conditions" = ifelse(
        amplicon == "Long",
        "initial denaturation:95degC_10min; denaturation:95degC_45s; annealing:59degC_45s; extension:72degC_90s; cycles:30; final extension:72degC_10min",
        "initial denaturation:95degC_10min; denaturation:95degC_60s; annealing:56degC_45s; extension:72degC_50s; cycles:35; final extension:72degC_3min"
      ),
      "sequencing method" = ifelse(
        amplicon == "Long",
        "PACBIO_SMRT",
        "ION_TORRENT, PACBIO_SMRT, ILLUMINA"
      )
    )
}

compile_blank_samples <- function(sample_data) {
  as(sample_data, "data.frame") %>%
  dplyr::filter(sample_type == "Blank", tech == "PacBio") %>%
  dplyr::left_join(
    readRDS("config/tags/all.rds") %>%
      tidyr::pivot_longer(
        cols = matches("(primer|barcode|adapter)_(seq|tag)_(fwd|rev)"),
        names_to = c(".value", "direction"),
        names_pattern = "(.+)_(fwd|rev)"
      ) %>%
      tidyr::unite("seq", ends_with("seq"), sep = "") %>%
      mutate_at(vars(ends_with("tag")), na_if, "unnamed") %>%
      mutate(tag = dplyr::coalesce(barcode_tag, primer_tag)) %>%
      dplyr::select(amplicon, well, direction, tag, seq) %>%
      filter(nchar(seq) > 0) %>%
      mutate_at("tag", str_replace, "its", "ITS") %>%
      mutate_at("tag", str_replace, "lr", "LR") %>%
      tidyr::unite("tag", tag, seq, sep = ": ") %>%
      group_by(amplicon, well, direction) %>%
      summarize(tag = paste(unique(tag), collapse = ", ")) %>%
      group_by(amplicon, well) %>%
      summarize(tag = paste(unique(tag), collapse = "; ")),
    by = c("amplicon", "well")) %>%
  dplyr::mutate(
    sample_alias = glue::glue("OT-blank-{plate}-{substr(amplicon, 1, 1)}"),
    tax_id = "2582415",
    scientific_name = "blank sample",
    common_name = "blank sample",
    sample_title = sample_alias,
    sample_description = "PCR blank",
    isolation_source = "Nuclease free water",
    "geographic location (country and/or sea)" = "Sweden",
    "geographic location (region and locality)" = "Uppsala",
    "project name" = "ena-STUDY-UPPSALA UNIVERISTY-19-03-2020-10:28:50:771-6",
    "target gene" = ifelse(amplicon == "Long", "ITS1, 5.8S, ITS2, LSU", "ITS2"),
    "pcr primers" = tag,
    "pcr conditions" = ifelse(
      amplicon == "Long",
      "initial denaturation:95degC_10min; denaturation:95degC_45s; annealing:59degC_45s; extension:72degC_90s; cycles:30; final extension:72degC_10min",
      "initial denaturation:95degC_10min; denaturation:95degC_60s; annealing:56degC_45s; extension:72degC_50s; cycles:35; final extension:72degC_3min"
    ),
    "sequencing method" = ifelse(
      amplicon == "Long",
      "PACBIO_SMRT",
      "ION_TORRENT, PACBIO_SMRT, ILLUMINA"
    )
  )
}

write_ena_samples <- function(samples, template_file, output_file) {
  template_names <- names(read_tsv(template_file, skip = 2))
  outdir <- dirname(output_file)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  file.copy(template_file, output_file, overwrite = TRUE)
  checkmate::assert_names(names(samples), must.include = template_names)
  select_at(samples, template_names) %>%
    write_tsv(output_file, append = TRUE, col_names = FALSE)
}

generate_pacbio_ion_manifest <- function(samples, short_files, long_files) {
  dplyr::inner_join(
    select(samples,
           sample_alias, "project name", plate, well, amplicon, sample_type),
    dplyr::bind_rows(
      short_files %>%
        mutate(plate = ifelse(
          tech == "Ion Torrent",
          list(c("001", "002")),
          as.list(plate)
        )) %>%
        unnest(plate),
      long_files
    ),
    by = c("plate", "well", "amplicon")
  ) %>%
    dplyr::transmute(
      STUDY = `project name`,
      SAMPLE = sample_alias,
      NAME = paste(sample_alias, tech, direction, sep = "-") %>% str_replace("-$", "") %>% str_replace(" ", ""),
      INSTRUMENT = paste(tech, machine) %>% str_replace("Ion S5", "S5"),
      LIBRARY_SOURCE = dplyr::case_when(
        sample_type == "Pos" ~ "GENOMIC",
        sample_type == "Sample" ~ "METAGENOMIC",
        sample_type == "Blank" ~ "OTHER"
      ),
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
          d$DESCRIPTION <- paste("also contains samples from", d$SAMPLE[2],
                                 "due to incomplete multiplexing")
          d$SAMPLE <- d$SAMPLE[1]
          d$NAME <- d$NAME[1]
          d$LIBRARY_SOURCE <- d$LIBRARY_SOURCE[1]
        }
        if (length(d$NAME) > 1) print(d)
        sink(file.path("output", "ENA", paste0(d$NAME, ".manifest")),
             append = FALSE)
        iwalk(d, ~ for (x in .x) {cat(.y, x); cat("\n")})
        sink()
      }
    )
}

generate_illumina_manifest <- function(samples, groups) {
  dplyr::inner_join(
    select(samples,
           sample_alias, "project name", plate, well, amplicon, sample_type),
    groups,
    by = c("plate", "well", "amplicon")
  ) %>%
    dplyr::transmute(
      STUDY = `project name`,
      SAMPLE = sample_alias,
      NAME = paste(sample_alias, tech, direction, sep = "-") %>%
        str_replace("-$", "") %>%
        str_replace(" ", ""),
      INSTRUMENT = paste(tech, machine),
      LIBRARY_SOURCE = dplyr::case_when(
        sample_type == "Pos" ~ "GENOMIC",
        sample_type == "Sample" ~ "METAGENOMIC",
        sample_type == "Blank" ~ "OTHER"
      ),
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
}

# get the ENA accession number from the submission receipt
get_reads_accno <- function(sample, reads_dir = "output/ENA/reads") {
  file.path(reads_dir, sample, "submit", "receipt.xml") %>%
    xml2::read_xml() %>%
    xml2::xml_find_first("RUN") %>%
    xml2::xml_attr("accession")
}

lookup_submitted_accnos <- function(path) {
  list.files(path = path, pattern = "receipt.xml", recursive = TRUE) %>%
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
    tidyr::pivot_wider(names_from = "tech", values_from = "accno")
}

group_by_taxon <- function(phylotax, all_labels = character(0)) {
  widen_taxonomy(phylotax$assigned) %>%
    dplyr::select(kingdom:genus, Taxonomy, label) %>%
    dplyr::full_join(tibble::tibble(label = all_labels), by = "label") %>%
    dplyr::mutate_at("Taxonomy", stringi::stri_replace_all_fixed, ";NA", "") %>%
    dplyr::mutate_at("Taxonomy", stringi::stri_replace_first_regex, ";[A-Za-z]+\\d+$", "") %>%
    unique() %>%
    tidyr::pivot_longer(cols = kingdom:genus, names_to = "rank", values_to = "taxon") %>%
    dplyr::filter(complete.cases(.), !stringr::str_detect(taxon, "\\d")) %>%
    dplyr::group_by(Taxonomy) %>%
    dplyr::summarize(
      taxon = dplyr::last(taxon),
      rank = dplyr::last(rank),
      label = paste(label, collapse = ",")
    )
}

# conversions for certain higher level taxa into the term used in the ENA
# taxonomy for environmental sequences
ena_taxon_replacements <- tibble::tribble(
  ~"original",          ~"ENA",
  "Fungi",              "fungus",
  "Alveolata",          "alveolate",
  "Apicomplexa",        "apicomplexan",
  "Gregarinasina",      "gregarine",
  "Neogregarinorida",   "gregarine",
  "Ciliophora",         "ciliate",
  "Haptorida",          "haptorid ciliate",
  "NA",                 "soil eukaryote",
  "Angiospermae",       "Magnoliophyta",
  "Stramenopila",       "Stramenopile",
  "Kickxellomycota",    "Kickxellomycotina",
  "Zoopagomycota",      "Zoopagomycotina",
  "Monoblepharomycota", "Monoblepharidomycetes",
  "Mortierellomycota",  "Mortierellomycotina",
  "Glomeromycota",      "Glomeromycotina",
  "Rozellomycota",      "Cryptomycota",
  "Bryopsida",          "Bryophyta",
  "Chromerida",         "Colpodellida",
  "Orbiliomycetes",     "Orbiliales",
  "Gigasporales",       "Diversisporales",
  "Paraglomeromycetes", "Paraglomerales",
  "Ichthyosporia",      "Ichthyosporea",
  "Symphyla",           "Arthropoda"
)

group_by_taxon_env <- function(phylotax, all_labels = character(0)) {
  widen_taxonomy(phylotax$assigned) %>%
    dplyr::select(kingdom:order, label) %>%
    dplyr::full_join(tibble::tibble(label = all_labels), by = "label") %>%
    dplyr::group_by_at(vars(kingdom:order)) %>%
    dplyr::summarize(label = paste(label, collapse = ",")) %>%
    tidyr::pivot_longer(cols = kingdom:order, names_to = "rank", values_to = "taxon") %>%
    dplyr::filter(complete.cases(.) | rank == "kingdom", !stringr::str_detect(taxon, "\\d") | is.na(taxon)) %>%
    dplyr::group_by(label) %>%
    dplyr::summarize(
      Taxonomy = paste(taxon, collapse = ";"),
      taxon = dplyr::last(paste(taxon)),
      rank = dplyr::last(rank)
    ) %>%
    dplyr::group_by(Taxonomy, taxon, rank) %>%
    dplyr::summarize(label = paste(label, collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_at(
      "taxon",
      stringi::stri_replace_all_regex,
      ena_taxon_replacements$original,
      ena_taxon_replacements$ENA,
      vectorize_all = FALSE
    ) %>%
    dplyr::mutate(taxon = ifelse(
      grepl("(Streptophyta|Metazoa)", Taxonomy) | taxon == "Viridiplantae",
      paste(taxon, "environmental sample"),
      paste("uncultured", taxon)
    ))
}

lookup_ncbi_taxon <- function(taxon, Taxonomy, rank, label, ...) {
  uid <- taxize::get_uid_(
    taxon,
    key = readLines("ENTREZ_KEY"),
    messages = FALSE
  )[[1]]
  if (is.null(uid)) {
    return(tibble::tibble(label = label, taxon = taxon, targetTaxonomy = Taxonomy, targetRank = rank))
  }
  uid$label <- label
  uid$taxon <- taxon
  uid$targetTaxonomy <- Taxonomy
  uid$targetRank <- rank
  uid$Taxonomy <- taxize::classification(
    uid$uid,
    db = "ncbi"
  ) %>%
    purrr::map_chr(~paste(.$name, collapse = ";"))
  uid$dist <- stringdist::stringdist(uid$Taxonomy, Taxonomy, method = "jaccard", q = 5)
  if (nrow(uid) <= 1) {
    return(uid)
  }
  # if (rank %in% uid$rank) {
  #   uid <- uid[uid$rank == rank, , drop = FALSE]
  #   if (nrow(uid) <= 1) {
  #     return(uid)
  #   }
  # }
  uid
}

ena_taxon_url <- "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"

lookup_ena_taxon <- function(taxon, Taxonomy, rank, label, ...) {
  url <- paste0(ena_taxon_url, curl::curl_escape(taxon))
  result <- httr::GET(url)
  result <- tryCatch(httr::content(result), error = function(e) NULL)
  if (is.null(result)) {
    return(
      tibble::tibble(
        label = label,
        taxon = taxon,
        targetTaxonomy = Taxonomy,
        targetRank = rank
      )
    )
  }
  result <- dplyr::bind_rows(result)
  result$label <- label
  result$taxon <- taxon
  result$targetTaxonomy <- Taxonomy
  result$targetRank <- rank
  result <- dplyr::rename(result, Taxonomy = lineage)
  result$Taxonomy <- stringr::str_replace_all(result$Taxonomy, "; ", ";")
  result$Taxonomy <- stringr::str_c(result$Taxonomy, result$scientificName)
  result$dist <- stringdist::stringdist(
    result$Taxonomy,
    Taxonomy,
    method = "jaccard",
    q = 5
  )
  result <- dplyr::arrange(result, dist)
  head(result, 1)
}

write_ENA_ASVs <- function(seqs, which_asvs, taxids, file, is_short) {
  # copy over the template to match the header
  writeLines("#template_accession ERT000009", file)
  # get all the labels and corresponding IDs.
  dplyr::select(taxids, label, taxId) %>%
    tidyr::separate_rows(label) %>%
    # choose only ASV labels belonging to this set
    dplyr::filter(label %in% which_asvs) %>%
    # add sequences
    dplyr::left_join(seqs, by = c("label" = "hash")) %>%
    dplyr::arrange(label) %>%
    # generate the submission data
    dplyr::transmute(
      ENTRYNUMBER = seq.int(nrow(.)),
      ORGANISM = taxId,
      CLONE = label,
      ENVSAM = "yes",
      ISOSOURCE = "soil",
      COUNTRY = "Benin",
      AREA = "Borgou",
      LOCALITY = "Foret Classee de l'Oueme Superieur",
      LATLON = "9.75 N 2.2 E",
      COLDATE = "2016",
      COLBY = "Brendan Furneaux",
      # sequences from the long amplicon will have ITS1 or LSU (and usually both)
      PFNAME1 = if (is_short) "gITS7" else "ITS1",
      PFSEQ1 = if (is_short) "GTGARTCATCGARTCTTTG" else "TCCGTAGGTGAACCTGC",
      PRNAME1 = if (is_short) "ITS4" else "LR5",
      PRSEQ1 = if (is_short) "TCCTCCGCTTATTGATATGC" else "TCCTGAGGGAAACTTCG",
      PRNAME2 = if (is_short) "ITS4A" else NULL,
      PRSEQ2 = if (is_short) "TCCTCGCCTTATTGATATGC" else NULL,
      # include the longest possible sequence
      SEQUENCE = chartr("U", "T", best),
      # which regions are present?
      `18S` = ifelse(is.na(ITS1) | is_short, "no", "yes"),
      ITS1 = ifelse(is.na(ITS1) | is_short, "no", "yes"),
      `5.8S` = ifelse(is.na(`5_8S`) & is.na(short), "no", "yes"),
      ITS2 = "yes",
      `28S` = ifelse(is.na(LSU1) & is.na(short), "no", "yes")
    ) %>%
    dplyr::filter(nchar(SEQUENCE) > 100) %>%
    # put SEQUENCE at the end where it belongs
    # (it needed to be calculated before ITS2 was changed from the ITS2 sequence to "yes")
    dplyr::select(!SEQUENCE, SEQUENCE) %T>%
    readr::write_delim(file, delim = "\t", append = TRUE, col_names = TRUE)
}
