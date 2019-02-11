library(magrittr)
library(dplyr)
library(tidyr)
library(dada2)
library(here)
library(drake)
library(assertthat)
library(assertr)
library(rlang)
library(readr)
library(stringr)
library(purrr)
library(glue)
library(FUNGuildR)
library(rITSx)

if (interactive()) {
  base.dir <- here()
  r.dir <- here("scripts")
  data.dir <- here("data")
  lab.dir <- here("config")
  seq.dir <- here("sequences")
  trim.dir <- file.path(seq.dir, "trim")
  region.dir <- file.path(seq.dir, "regions")
  filter.dir <- file.path(seq.dir, "filter")
  ref.dir <- here("reference")
  out.dir <- here("output")
  in.files <- list.files(path = trim.dir,
                         pattern = "\\.trim\\.fastq\\.gz$",
                         full.names = TRUE)
  dataset.file <- file.path(lab.dir, "datasets.csv")
  regions.file <- file.path(lab.dir, "regions.csv")
  ncpu <- max(1, parallel::detectCores() - 1)
} else {
  r.dir <- Sys.getenv("RDIR")
  dataset.file <- Sys.getenv("DATASET")
  regions.file <- Sys.getenv("REGIONS")
  target <- Sys.getenv("TARGETLIST")
  ncpu <- as.integer(Sys.getenv("NCORES"))
  trim.dir <- Sys.getenv("TRIMDIR")
  region.dir <- Sys.getenv("REGIONDIR")
  filter.dir <- Sys.getenv("FILTERDIR")
  ref.dir <- Sys.getenv("REF_ROOT")
  out.dir <- Sys.getenv("OUTDIR")
  # get the prereq files on stdin
  # there are too many to pass as an environmental variable!
  con <- file("stdin")
  open(con, blocking = TRUE)
  prereqs <- readLines(con)
  close(con)
  prereqs <- str_split(prereqs, " ") %>% unlist
  cat("prereqs: ", prereqs, "\n")
  
  in.files <- str_subset(prereqs,
                         pattern = "\\.trim\\.fastq\\.gz$")
  cat("in.files: ", in.files, "\n")
}

source(file.path(r.dir, "combine_derep.R"))
source(file.path(r.dir, "extract_regions.R"))
source(file.path(r.dir, "dada.R"))


datasets <- read_csv(dataset.file)
regions <- read_csv(regions.file)

# Generate the list of trimmed, demultiplexed files we expect
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
  unnest %>%
  mutate(Trim.File = glue("{Seq.Run}_{Plate}-{Well}{Direction}.trim.fastq.gz"),
         FID = glue("{Seq.Run}{Plate}{Well}{Direction}") %>%
           str_replace_all("[:punct:]", ""),
         positions = make.names(glue("positions_{FID}")),
         join_derep = make.names(glue("join_derep_.{Primer.pair}.")),
         itsxtrim = make.names(glue("itsxtrim_.{Primer.pair}."))) %>%
  filter(file.exists(file.path(trim.dir, Trim.File))) %>%
  verify(Trim.File %in% basename(in.files)) %>%
  mutate_at(c("join_derep", "itsxtrim", "FID"), syms)
saveRDS(meta1, "meta1.rds")

meta2 <- meta1 %>%
  mutate_at("Regions", str_split, ",") %>%
  unnest(Region = Regions) %>%
  left_join(regions) %>%
  mutate(RID = glue("{Seq.Run}{Plate}{Well}{Direction}{Region}"),
         PID = glue("{Seq.Run}_{Plate}_{Region}"), 
         Region.File = glue("{Seq.Run}_{Plate}-{Well}{Direction}-{Region}.trim.fastq.gz"),
         Filter.File = glue("{Seq.Run}_{Plate}-{Well}{Direction}-{Region}.qfilt.fastq.gz")) %>%
  mutate_at(c("positions", "RID", "PID"), syms)
saveRDS(meta2, "meta2.rds")

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
}

meta3 <- meta2 %>%
  group_by(Seq.Run, Plate, Region) %>%
  summarize(Filter.File = paste(as.character(Filter.File), collapse = " ")) %>%
  ungroup %>%
  left_join(datasets %>% select(-Regions)) %>%
  left_join(regions) %>%
  mutate(PID = glue("{Seq.Run}_{Plate}_{Region}")) %>%
  mutate_at(c("PID", "Region"), syms)
saveRDS(meta3, "meta3.rds")

meta4 <- regions %>%
  mutate(Reference = strsplit(Reference, " *, *") %>%
                                map(unlist)) %>%
  unnest(Reference, .preserve = Region) %>%
  mutate(TID = glue("{Region}_{Reference}"),
         big_seq_table = glue("big_seq_table_{Region}"),
         Reference.File = Reference) %>%
  mutate_at(c("TID", "big_seq_table", "Region"), syms)
saveRDS(meta4, "meta4.rds")

plan <- drake_plan(
  
  datasets = read_csv(file_in(!!dataset.file)),
  
  derep1 = target(
    tibble(file = Trim.File,
           derep = list(derepFastq(file_in(!!file.path(trim.dir, Trim.File))))),
    transform = map(.data = !!select(meta1, Trim.File, FID, Primer.pair),
                    .id = FID)),
  
  join_derep = target(
    combine_derep(derep1),
    transform = combine(derep1, .by = Primer.pair)),
  
  itsxtrim = target(
    itsx(in.file = join_derep$fasta,
         read_function = Biostrings::readDNAStringSet,
         fasta = FALSE, summary = FALSE, graphical = FALSE,
         save_regions = "none", not_found = FALSE,
         complement = FALSE, cpu = ignore(!!ncpu)),
    transform = map(join_derep, Primer.pair, .id = Primer.pair)),
  
  positions = target(
    filter(join_derep$map, file == Trim.File) %>%
      left_join(itsxtrim$positions, by = c("newmap" = "seq")),
    transform = map(.data = !!select(meta1, join_derep, itsxtrim, Trim.File, FID),
                    .id = FID)),
  
  regions = target(
    extract_region(infile = file_in(!!file.path(trim.dir, Trim.File)),
                   outfile = file_out(!!file.path(region.dir, Region.File)),
                   region = Region,
                   positions = positions),
    transform = map(.data = !!select(meta2, Trim.File, Region.File, Region,
                                     positions, RID),
                    .id = RID)),
  
  filter = target({
    outfile <- file_out(!!file.path(filter.dir, Filter.File))
    unlink(outfile)
    dada2::filterAndTrim(fwd = file_in(!!file.path(region.dir, Region.File)),
                         filt = file_out(!!file.path(filter.dir, Filter.File)),
                         maxLen = MaxLength,
                         minLen = MinLength,
                         maxEE = MaxEE)
    if (!file.exists(outfile)) {
      ShortRead::writeFastq(ShortRead::ShortReadQ(),
                            outfile)
    }},
    transform = map(.data = !!select(meta2, Region.File, Filter.File,
                                    MaxLength, MinLength, MaxEE, RID),
                    .id = RID)),
  
  derep2 = target({
    infiles <- file_in(!!file.path(filter.dir, unlist(str_split(Filter.File, " "))))
    infiles <- purrr::discard(infiles, ~file.size(.) < 50)
    dada2::derepFastq(fls = infiles)},
    transform = map(.data = !!meta3,
                    .id = PID)),
  
  err = target({
    err.fun <- ifelse(Tech == "PacBio",
                      PacBioErrfun,
                      loessErrfun)
    learnErrors(fls = derep2, nbases = 1e8,
                multithread=ignore(!!ncpu), randomize=TRUE,
                errorEstimationFunction = err.fun,
                verbose = TRUE)},
    transform = map(derep2, Tech, PID, .id = PID)),
  
  dada = target(
    dada(derep = derep2, err = err,
         multithread = ignore(!!ncpu), 
         HOMOPOLYMER_GAP_PENALTY = eval(parse_expr(Homopolymer.Gap.Penalty)),
         BAND_SIZE = Band.Size,
         pool = eval(parse_expr(Pool))),
    transform = map(derep2, err, PID, Homopolymer.Gap.Penalty, Band.Size, Pool,
                    .id = PID)),
  
  dada_map = target(
    dadamap(derep2, dada),
    transform = map(derep2, dada)),
  
  seq_table = target(
    makeSequenceTable(dada),
    transform = map(dada, PID, .id = PID)),
  
  nochim = target(
    dada2::removeBimeraDenovo(seq_table, method = "consensus",
                              multithread = ignore(!!ncpu)),
    transform = map(seq_table, PID, .id = PID)),
  
  big_seq_table = target(
    join_seqs(nochim),
    transform = combine(nochim, .by = Region)),
  
  taxon = target(
    taxonomy(big_seq_table,
             reference = file_in(!!file.path(ref.dir, paste0(Reference, ".fasta.gz"))),
             multithread = ignore(!!ncpu)),
    transform = map(.data = !!meta4, .id = TID)),
  
  funguild_db = get_funguild_db(),
  
  guilds_table = target(
    funguild_assign(taxon, db = funguild_db),
    transform = map(taxon, TID, .id = TID)),
  seq_count = target(
    tibble(
      file = file_in(fq.file),
      reads = system(glue("zcat {file} | grep -c '^@'", file = fq.file),
                     intern = TRUE) %>%
        as.integer),
    transform = map(fq.file = !!c(file.path(trim.dir, meta1$Trim.File),
                                 file.path(region.dir, meta2$Region.File),
                                 file.path(filter.dir, meta2$Filter.File)))),
  seq_counts = target(
    bind_rows(seq_count),
    transform = combine(seq_count)),

  qstats = target(
    q_stats(file_in(!!file.path(region.dir, Region.File))),
    transform = map(.data = !!meta2,
                    .id = RID)),
  qstats_join = target(
    bind_rows(qstats),
    transform = combine(qstats)),
  qstats_knit = {
    if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
    rmarkdown::render(
      knitr_in(!!file.path(r.dir, "qual-check.Rmd")),
      output_file = file_out(!!file.path(out.dir, "qual-check.pdf")))},
  trace = TRUE
)

saveRDS(plan, "plan.rds")

drake_plan_source(plan)
#predict_runtime(dconfig, jobs = ncpu)
if (interactive()) {
  dconfig <- drake_config(plan)
  vis_drake_graph(dconfig)
}
future::plan("multiprocess")
cat("\n First drake::make... (future)\n")
# make embarassing targets at the beginning
make(plan,
     parallelism = "future",
     jobs = ncpu, jobs_preprocess = ncpu,
     caching = "worker",
     targets = str_subset(plan$target, "^join_derep_")
)
cat("\n Second drake::make... (loop)\n")
# itsx is internally parallel
make(plan,
     parallelism = "loop",
     jobs = 1, jobs_preprocess = ncpu,
     caching = "worker",
     targets = str_subset(plan$target, "^itsxtrim_")
)
# embarrasing targets after itsx
cat("\n Third drake::make... (future)\n")
make(plan,
     parallelism = "future",
     jobs = ncpu, jobs_preprocess = ncpu,
     caching = "worker",
     targets = c(str_subset(plan$target, "^derep2_"), "qstats_knit", "seq_counts")
)

# dada is internally parallel
cat("\n Fourth drake::make... (loop)\n")
make(plan,
     parallelism = "loop",
     jobs = 1, jobs_preprocess = ncpu,
     caching = "worker",
     targets = str_subset(plan$target, "^taxon_")
)
# finish up parallel
cat("\n Fifth drake::make... (future)\n")
make(plan,
     parallelism = "future",
     jobs = ncpu, jobs_preprocess = ncpu,
     caching = "worker"
)
#   times[[ncpu]] <- tictoc::toc()
# }
if (interactive()) vis_drake_graph(dconfig)

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
