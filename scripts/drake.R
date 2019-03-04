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
  ref.dir <- here("reference")
  rmd.dir <- here("writing")
  out.dir <- here("output")
  in.files <- list.files(path = trim.dir,
                         pattern = "\\.trim\\.fastq\\.gz$",
                         full.names = TRUE)
  dataset.file <- file.path(lab.dir, "datasets.csv")
  regions.file <- file.path(lab.dir, "regions.csv")
  platemap.file <- file.path(lab.dir, "Brendan_soil2.xlsx")
  platemap.sheet <- "Concentration samples"
  plan.file <- "plan.rds"
  meta1.file <- "meta1.rds"
  meta2.file <- "meta2.rds"
  meta3.file <- "meta3.rds"
  meta4.file <- "meta4.rds"
  meta4.csv <- "meta4.csv"
  pid.file <- "pids.txt"
  tid.file <- "tids.txt"
  drakedata.file <- "drake.Rdata"
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
  rmd.dir <- snakemake@config$rmddir
  out.dir <- snakemake@config$outdir
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
  cat("Creating drake plan in command line session...\n")
  r.dir <- Sys.getenv("RDIR")
  dataset.file <- Sys.getenv("DATASET")
  regions.file <- Sys.getenv("REGIONS")
  platemap.file <- Sys.getenv("PLATEMAPFILE")
  platemap.sheet <- Sys.getenv("PLATEMAPSHEET")
  target <- Sys.getenv("TARGETLIST")
  trim.dir <- Sys.getenv("TRIMDIR")
  region.dir <- Sys.getenv("REGIONDIR")
  filter.dir <- Sys.getenv("FILTERDIR")
  ref.dir <- Sys.getenv("REF_ROOT")
  rmd.dir <- Sys.getenv("RMDDIR")
  out.dir <- Sys.getenv("OUTDIR")
  # get the prereq files on stdin
  # there are too many to pass as an environmental variable!
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  prereqs <- stringr::str_split(prereqs, " ") %>% unlist
  
  in.files <- stringr::str_subset(prereqs,
                         pattern = "\\.trim\\.fastq\\.gz$")
}

library(magrittr)
library(tidyverse)
library(rlang)
library(glue)
library(drake)
library(assertr)
library(parallel)

#### parallel setup ####

# are we running slurm?
is_slurm <- function() nchar(Sys.which("sbatch")) > 0
is_local <- function() !is_slurm()

# how many cpus do we have on the local machine?
# if we're not running on the cluster, leave one cpu free.
local_cpus <- function() {
  if (is_slurm()) {
    out <- as.integer(Sys.getenv("CORES_PER_TASK"))
    assertthat::assert_that(assertthat::is.count(out))
    out
  } else {
    max(parallel::detectCores() - 1, 1)
  }
}

# how many cpus can we get on one node?
max_cpus <- function() {
  if (is_slurm()) {
    slurminfo <- system("sinfo -N --long", intern = TRUE)[-1]
    slurmcolnames <- strsplit(slurminfo[1], "\\s+")[[1]]
    slurminfo <- readr::read_fwf(slurminfo,
                          readr::fwf_empty(slurminfo, col_names = slurmcolnames),
                          skip = 1)
    return(max(slurminfo[["CPUS"]]))
  } else {
    return(local_cpus())
  }
}

#### If we were passed a target, what is it?
get_target <- function(default) {
  if (interactive()) {
    return(default)
  } else if (exists("snakemake")) {
    return(sub("^.", "", snakemake@output))
  } else {
    return(sub("^.", "", Sys.getenv("TARGETLIST")))
  }
}

#### what is our logfile called?
get_logfile <- function(default) {
  if (interactive()) {
    return(file.path("logs", paste0(default, ".log")))
  } else if (exists("snakemake")) {
    return(file(snakemake@log[[1]], open = "at"))
  } else {
    return(file(file.path(Sys.getenv("LOGDIR"), paste0(default, ".log"))))
  }
}

#### Start logging ####
setup_log <- function(default) {
  logfile <- get_logfile(default)
  sink(logfile, split = TRUE)
  sink(logfile, type = "message")
}

#### read scripts and configs ####
# load various pipeline functions
source(file.path(r.dir, "combine_derep.R"))
source(file.path(r.dir, "extract_regions.R"))
source(file.path(r.dir, "dada.R"))
source(file.path(r.dir, "plate_check.R"))

# load the dataset and region definitions
datasets <- read_csv(dataset.file)
regions <- read_csv(regions.file)

# "meta" dataframes containing definitions for each instance of mapped targets

#### meta1 ####
# meta1 has one row per demultiplexed, primer-trimmed fastq.gz file
# The "index" is FID (File ID)
meta1 <- datasets %>%
  mutate(Primer.pair = paste(str_replace(Forward, "_tag.*$", ""),
                            str_replace(Reverse, "_tag.*$", ""),
                            sep = "-"),
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
         FID = glue("{Seq.Run}{Plate}{Well}{Direction}") %>%
           str_replace_all("[:punct:]", ""),
         positions = make.names(glue("positions_{FID}")),
         join_derep = make.names(glue("join_derep_.{Primer.pair}.")),
         itsxtrim = make.names(glue("itsxtrim_.{Primer.pair}."))) %>%
  filter(file.exists(file.path(trim.dir, Trim.File))) %>%
  filter(file.size(file.path(trim.dir, Trim.File)) > 40) %>%
  verify(file.path(trim.dir, Trim.File) %in% in.files) %>%
  mutate_at(c("join_derep", "itsxtrim", "FID"), syms)

#### meta2 ####
# meta2 has one row per region per input file
# The index is RID (Region ID)
# PID is defined because it will be the index for meta3
meta2 <- meta1 %>%
  mutate_at("Regions", str_split, ",") %>%
  unnest(Region = Regions) %>%
  left_join(regions, by = "Region") %>%
  mutate(RID = glue("{Seq.Run}{Plate}{Well}{Direction}{Region}"),
         PID = glue("{Seq.Run}_{Plate}_{Region}"), 
         Region.File = glue("{Seq.Run}_{Plate}-{Well}{Direction}-{Region}.trim.fastq.gz"),
         Filter.File = glue("{Seq.Run}_{Plate}-{Well}{Direction}-{Region}.qfilt.fastq.gz")) %>%
  mutate_at(c("positions", "RID", "PID"), syms)

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
  saveRDS(meta1, "meta1.rds")
  saveRDS(meta2, "meta2.rds")
}

#### meta3 ####
# meta3 has one row per region per plate
# The index is PID (Plate ID)
meta3 <- meta2 %>%
  select(Seq.Run, Plate, Region, PID) %>%
  unique() %>%
  mutate_at("PID", as.character) %>%
  arrange(PID) %>%
  mutate(region_derep = glue("region_derep_{PID}")) %>%
  left_join(datasets %>% select(-Regions), by = "Seq.Run") %>%
  left_join(regions, by = "Region") %>%
  mutate_at(c("region_derep", "Region", "PID"), syms)
if (!interactive()) {
  saveRDS(meta3, "meta3.rds")
  readr::write_lines(meta3$PID, "pids.txt")
}

#### meta4 ####
# meta4 has one row per region
# The index is TID (Taxonomy ID)
meta4 <- meta3 %>%
  mutate_if(is.list, as.character) %>%
  group_by_at(names(regions)) %>%
  summarize(PID = paste(PID, collapse = ",")) %>%
  ungroup() %>%
  mutate(Reference = strsplit(Reference, " *, *") %>%
                                map(unlist)) %>%
  unnest(Reference, .preserve = Region) %>%
  mutate(TID = glue("{Region}_{Reference}"),
         big_seq_table = glue("big_seq_table_{Region}"),
         Reference.File = Reference) %>%
  arrange(TID) %>%
  mutate_at(c("TID", "big_seq_table", "Region"), syms)
if (!interactive()) {
  saveRDS(meta4, "meta4.rds")
  mutate_if(meta4, is.list, as.character) %>%
  readr::write_csv("meta4.csv")
  readr::write_lines(meta4$TID, "tids.txt")
}

#### drake plan ####
cat("\nbuilding plan...\n")
tictoc::tic()
plan <- drake_plan(
  
  datasets = readr::read_csv(file_in(!!dataset.file)),
  
  # dereplicate input sequences before running ITSx
  # this is done independently in parallel
  derep1 = target(
    tibble::tibble(
      file = Trim.File,
      derep = list(dada2::derepFastq(file_in(!!file.path(trim.dir, Trim.File)),
                                     qualityType = "FastqQuality",
                                     n = 1e4) %>% inset2("quals", NULL))),
    transform = map(.data = !!select(meta1, Trim.File, FID, Primer.pair),
                    .tag_in = step,
                    .id = FID)),
  
  # join the results from dereplicating into a list of unique sequences
  # and a map from the original sequences to the uniques
  join_derep = target(
    combine_derep(derep1),
    transform = combine(derep1,
                        .tag_in = step,
                        .by = Primer.pair)),
  
  # break the list of unique files into equal sized chunks for ITSx
  split_fasta = target(
    split(join_derep$fasta, seq_along(join_derep$fasta) %% !!bigsplit + 1),
    transform = map(join_derep,
                    .tag_in = step, .id = Primer.pair)),
  
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
                      .id = Primer.pair)),
  
  # combine the chunks into a master positions list.
  itsxtrim = target(
    dplyr::bind_rows(itsx_shard),
    transform = combine(itsx_shard,
                        .tag_in = step,
                        .by = Primer.pair)),
  
  # find the entries in the positions list which correspond to the sequences
  # in each file
  positions = target(
    dplyr::filter(join_derep$map, file == Trim.File) %>%
      dplyr::left_join(itsxtrim, by = c("newmap" = "seq")),
    transform = map(.data = !!select(meta1, join_derep, itsxtrim, Trim.File,
                                     FID),
                    .tag_in = step,
                    .id = FID)),
  
  # cut the required regions out of each sequence
  regions = target(
    extract_region(infile = file_in(!!file.path(trim.dir, Trim.File)),
                   outfile = file_out(!!file.path(region.dir, Region.File)),
                   region = Region,
                   positions = positions),
    transform = map(.data = !!select(meta2, Trim.File, Region.File, Region,
                                     positions, RID),
                    .tag_in = step,
                    .id = RID)),
  
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
                                     MaxLength, MinLength, MaxEE, RID),
                    .tag_in = step,
                    .id = RID)),
  
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
                  .id = RID)),
  
  region_derep = target(
    purrr::compact(c(derep2)),
    transform = combine(derep2, .by = PID, .tag_in = step),
  ),
  
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
                    .tag_in = step, .id = PID)),
  
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
                    .id = PID)),
  
  # Make maps from individual sequences in source files to dada ASVs
  dada_map = target(
    dadamap(region_derep, dada),
    transform = map(dada, .data = !!meta3,
                    .tag_in = step, .id = PID)),
  
  # Make a sample x ASV abundance matrix
  seq_table = target(
    dada2::makeSequenceTable(dada),
    transform = map(dada, .data = !!meta3, .tag_in = step, .id = PID)),
  
  # Remove likely bimeras
  nochim = target(
    dada2::removeBimeraDenovo(seq_table, method = "consensus",
                              multithread = ignore(dada_cpus)),
    transform = map(seq_table, .data = !!meta3, .tag_in = step, .id = PID)),
  
  # Join all the sequence tables for each region
  big_seq_table = target(
    dada2::mergeSequenceTables(tables = list(nochim)),
    transform = combine(nochim, .tag_in = step, .by = Region)),
  
  # Assign taxonomy to each ASV
  taxon = target(
    taxonomy(big_seq_table,
             reference = file_in(!!file.path(ref.dir, paste0(Reference, ".fasta.gz"))),
             multithread = ignore(dada_cpus)),
    transform = map(.data = !!meta4, .tag_in = step, .id = TID)),
  
  # Download the FUNGuild database
  funguild_db = FUNGuildR::get_funguild_db(),
  
  # Assign ecological guilds to the ASVs based on taxonomy.
  guilds_table = target(
    FUNGuildR::funguild_assign(taxon, db = funguild_db),
    transform = map(taxon, TID, .tag_in = step, .id = TID)),
  
  platemap = read_platemap(file_in(!!platemap.file), platemap.sheet),
  
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
  
  # merge all the sequence counts into one data.frame
  seq_counts = target(
    dplyr::bind_rows(seq_count),
    transform = combine(seq_count)),

  # calculate quality stats for the different regions
  qstats = target(
    q_stats(file_in(!!file.path(region.dir, Region.File))),
    transform = map(.data = !!meta2,
                    .tag_in = step,
                    .id = RID)),
  # join all the quality stats into one data.frame
  qstats_join = target(
    dplyr::bind_rows(qstats) %>%
      tidyr::extract(
        col = "file",
        into = c("Seq.Run", "Plate", "Well", "Direction", "Region"),
        regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.trim\\.fastq\\.gz") %>%
      dplyr::left_join(datasets),
    transform = combine(qstats)),
  # knit a report about the quality stats.
  qstats_knit = {
    if (!dir.exists(!!out.dir)) dir.create(!!out.dir, recursive = TRUE)
    rmarkdown::render(
      knitr_in(!!file.path(rmd.dir, "qual-check.Rmd")),
      output_file = file_out(!!file.path(out.dir, "qual-check.pdf")),
      output_dir = !!out.dir)},
  trace = TRUE
)
tictoc::toc()

if (!interactive()) saveRDS(plan, "plan.rds")

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
     file = "drake.Rdata")

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
