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
  r.dir <- here("R")
  data.dir <- here("data")
  lab.dir <- file.path(data.dir, "lab_setup")
  seq.dir <- here("sequences")
  trim.dir <- file.path(seq.dir, "trim")
  region.dir <- file.path(seq.dir, "regions")
  filter.dir <- file.path(seq.dir, "filter")
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

datasets <- read_csv(dataset.file)
regions <- read_csv(regions.file)

base.dadaOpt <- list(pool = TRUE)
optreplace <- function(x) list_modify(base.dadaOpt, !!!x)

# Generate the list of trimmed, demultiplexed files we expect
meta1 <- datasets %>%
  mutate(PrimerPair = paste(str_replace(Forward, "_tag.*$", ""),
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
  mutate(TrimFile = glue("{Seq.Run}_{Plate}-{Well}{Direction}.trim.fastq.gz"),
         fid_ = glue("{Seq.Run}{Plate}{Well}{Direction}") %>%
           str_replace_all("[:punct:]", ""),
         positions_ = make.names(glue("positions_{fid_}")),
         join_derep_ = make.names(glue("join_derep_.{PrimerPair}.")),
         itsxtrim_ = make.names(glue("itsxtrim_.{PrimerPair}."))) %>%
  filter(file.exists(file.path(trim.dir, TrimFile))) %>%
  verify(TrimFile %in% basename(in.files))
saveRDS(meta1, "meta1.rds")

meta2 <- meta1 %>%
  mutate_at("Regions", str_split, ",") %>%
  unnest(Region = Regions) %>%
  left_join(regions) %>%
  mutate(rid_ = glue("{fid_}{Region}"),
         pid_ = glue("{Seq.Run}_{Plate}_{Region}"), 
         regionfile_ = glue("{Seq.Run}_{Plate}-{Well}{Direction}-{Region}.trim.fastq.gz"),
         filterfile_ = glue("{Seq.Run}_{Plate}-{Well}{Direction}-{Region}.qfilt.fastq.gz"))
saveRDS(meta2, "meta2.rds")

if (interactive()) {
  #use a subset
  meta1 %<>%
    group_by(Seq.Run, Plate) %>%
    mutate(filesize = file.size(file.path(trim.dir, TrimFile))) %>%
    arrange(desc(filesize), .by_group = TRUE) %>%
    dplyr::slice(1:4) %>%
    ungroup()
  
  meta2 %<>%
    semi_join(meta1, by = "TrimFile")
}

meta3 <- meta2 %>%
  group_by(Seq.Run, Plate, Region) %>%
  summarize(filterfile_ = paste(as.character(filterfile_), collapse = " ")) %>%
  left_join(datasets %>% select(-Regions)) %>%
  left_join(regions) %>%
  mutate(pid_ = glue("{Seq.Run}_{Plate}_{Region}"))
saveRDS(meta3, "meta3.rds")

plan <- drake_plan(
  datasets = read_csv(file_in(!!dataset.file)),
  derep1 = target(
    tibble(file = f,
           derep = list(derepFastq(file_in(!!file.path(trim.dir, f))))),
    transform = map(f = !!meta1$TrimFile,
                    fid = !!syms(meta1$fid_),
                    primer.pair = !!meta1$PrimerPair,
                    .id = fid)),
  
  join_derep = target(
    combine_derep(derep1),
    transform = combine(derep1, .by = primer.pair)),
  
  itsxtrim = target(
    itsx(in.file = join_derep$fasta,
         read_function = Biostrings::readDNAStringSet,
         fasta = FALSE, summary = FALSE, graphical = FALSE,
         save_regions = "none", not_found = FALSE,
         complement = FALSE, cpu = ignore(ncpu)),
    transform = map(join_derep, primer.pair, .id = primer.pair)),
  
  positions = target(
    filter(join_derep$map, file == f) %>%
      left_join(itsxtrim$positions, by = c("newmap" = "seq")),
    transform = map(join_derep = !!syms(meta1$join_derep_),
                    itsxtrim = !!syms(meta1$itsxtrim_),
                    primer.pair = !!meta1$PrimerPair,
                    f = !!meta1$TrimFile,
                    fid = !!syms(meta1$fid_),
                    .id = fid)),
  
  regions = target(
    extract_region(infile = file_in(!!file.path(trim.dir, trimf)),
                   outfile = file_out(!!file.path(region.dir, regf)),
                   region = r,
                   positions = positions),
    transform = map(trimf = !!meta2$TrimFile,
                    regf = !!meta2$regionfile_,
                    qfilf = !!meta2$filterfile_,
                    r = !!meta2$Region,
                    positions = !!syms(meta2$positions_),
                    max_length = !!meta2$MaxLength,
                    min_length = !!meta2$MinLength,
                    maxEE = !!meta2$MaxEE,
                    rid = !!syms(meta2$rid_),
                    .id = rid)),
  
  filter = target({
    outfile <- file_out(!!file.path(filter.dir, qfilf))
    unlink(outfile)
    dada2::filterAndTrim(fwd = file_in(!!file.path(region.dir, regf)),
                         filt = file_out(!!file.path(filter.dir, qfilf)),
                         maxLen = max_length,
                         minLen = min_length,
                         maxEE = max_EE)
    if (!file.exists(outfile)) {
      ShortRead::writeFastq(ShortRead::ShortReadQ(),
                            outfile)
    }},
    transform = map(trimf = !!meta2$TrimFile,
                    regf = !!meta2$regionfile_,
                    qfilf = !!meta2$filterfile_,
                    r = !!meta2$Region,
                    positions = !!syms(meta2$positions_),
                    max_length = !!meta2$MaxLength,
                    min_length = !!meta2$MinLength,
                    max_EE = !!meta2$MaxEE,
                    rid = !!syms(meta2$rid_),
                    .id = rid)),
  
  derep2 = target({
    infiles <- file_in(!!file.path(filter.dir, unlist(str_split(filtf, " "))))
    infiles <- purrr::discard(infiles, ~file.size(.) < 50)
    dada2::derepFastq(fls = infiles)},
    transform = map(filtf = !!meta3$filterfile_,
                    Tech = !!meta3$Tech,
                    HGP = !!meta3$Homopolymer.Gap.Penalty,
                    BandSize = !!meta3$Band.Size,
                    Pool = !!meta3$Pool,
                    pid = !!syms(meta3$pid_),
                    .id = pid)),
  
  err = target({
    err.fun <- ifelse(Tech == "PacBio",
                      PacBioErrfun,
                      loessErrfun)
    learnErrors(fls = derep2, nbases = 1e8,
                multithread=ignore(ncpu), randomize=TRUE,
                errorEstimationFunction = err.fun,
                verbose = TRUE)},
    transform = map(derep2,
                    Tech,
                    pid,
                    .id = pid)),
  
  dada = target(
    dada(derep = derep2, err = err,
         multithread = ignore(ncpu), 
         HOMOPOLYMER_GAP_PENALTY = eval(parse_expr(HGP)),
         BAND_SIZE = BandSize,
         pool = eval(parse_expr(Pool))),
    transform = map(derep2, err, HGP, BandSize, Pool, pid, .id = pid)),
  
  qstats = target(
    q_stats(file_in(!!file.path(region.dir, f))),
    transform = map(f = !!meta2$regionfile_,
                    rid = !!syms(meta2$rid_),
                    .id = rid)),
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
     targets = str_subset(plan$target, "^derep2_")
)
# itsx is internally parallel
cat("\n Fourth drake::make... (loop)\n")
make(plan,
     parallelism = "loop",
     jobs = 1, jobs_preprocess = ncpu,
     caching = "worker",
     targets = str_subset(plan$target, "^dada_")
)
# finish up parallel
cat("\n Fifth drake::make... (future)\n")
make(plan,
     parallelism = "future",
     jobs = ncpu, jobs_preprocess = ncpu,
     caching = "worker",
     targets = str_subset(plan$target, "^derep2_")
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
