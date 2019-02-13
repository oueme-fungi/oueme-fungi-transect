library(magrittr)
library(dplyr)
library(tidyr)
library(dada2)
library(here)
library(drake)
library(future)
library(future.batchtools)
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
  ncpu <- as.integer(Sys.getenv("CORES_PER_TASK"))
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
#  cat("prereqs: ", prereqs, "\n")
  
  in.files <- str_subset(prereqs,
                         pattern = "\\.trim\\.fastq\\.gz$")
#  cat("in.files: ", in.files, "\n")
}

cat("ncpu = ", ncpu, "\n")
is_slurm = nchar(Sys.which("sbatch")) > 0

source(file.path(r.dir, "combine_derep.R"))
source(file.path(r.dir, "extract_regions.R"))
source(file.path(r.dir, "dada.R"))


datasets <- read_csv(dataset.file)
regions <- read_csv(regions.file)

# "meta" dataframes containing definitions for each instance of mapped targets

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
  unnest %>%
  mutate(Trim.File = glue("{Seq.Run}_{Plate}-{Well}{Direction}.trim.fastq.gz"),
         FID = glue("{Seq.Run}{Plate}{Well}{Direction}") %>%
           str_replace_all("[:punct:]", ""),
         positions = make.names(glue("positions_{FID}")),
         join_derep = make.names(glue("join_derep_.{Primer.pair}.")),
         itsxtrim = make.names(glue("itsxtrim_.{Primer.pair}."))) %>%
  filter(file.exists(file.path(trim.dir, Trim.File))) %>%
  filter(file.size(file.path(trim.dir, Trim.File)) > 40) %>%
  verify(Trim.File %in% basename(in.files)) %>%
  mutate_at(c("join_derep", "itsxtrim", "FID"), syms)

# meta2 has one row per region per input file
# The index is RID (Region ID)
# PID is defined because it will be the index for meta3
meta2 <- meta1 %>%
  mutate_at("Regions", str_split, ",") %>%
  unnest(Region = Regions) %>%
  left_join(regions) %>%
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

# meta3 has one row per region per plate
# The index is PID (Plate ID)
meta3 <- meta2 %>%
  group_by(Seq.Run, Plate, Region) %>%
  summarize(Filter.File = paste(as.character(Filter.File), collapse = " ")) %>%
  ungroup %>%
  left_join(datasets %>% select(-Regions)) %>%
  left_join(regions) %>%
  mutate(PID = glue("{Seq.Run}_{Plate}_{Region}")) %>%
  mutate_at(c("PID", "Region"), syms)
if (!interactive()) saveRDS(meta3, "meta3.rds")

# meta4 has one row per region
# The index is TID (Taxonomy ID)
meta4 <- regions %>%
  mutate(Reference = strsplit(Reference, " *, *") %>%
                                map(unlist)) %>%
  unnest(Reference, .preserve = Region) %>%
  mutate(TID = glue("{Region}_{Reference}"),
         big_seq_table = glue("big_seq_table_{Region}"),
         Reference.File = Reference) %>%
  mutate_at(c("TID", "big_seq_table", "Region"), syms)
if (!interactive()) saveRDS(meta4, "meta4.rds")

# the hmmersearch in itsx is very slow, but super parallel.
# we can break it into a lot of pieces and run it on several nodes.
# snowy has 16 cores per node, rackham has 20.  80 is the LCM.
# (5 nodes on snowy, 4 nodes on rackham).
bigsplit = 80

# dada2 commands can use internal parallelism, but with diminishing returns.
# instead of devoting a large number of cores to each, we will run them
# with only a fraction of the node, in several parallel jobs.
dadacores = min(4, ncpu)
dadajobs = ncpu %/% dadacores

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
  
  split_fasta = target(
    split(join_derep$fasta, seq_along(join_derep$fasta) %% !!bigsplit + 1),
    transform = map(join_derep, Primer.pair, .id = Primer.pair)),
  
  itsx_shard = target(
    itsx(in.file = split_fasta[[shard]],
         read_function = Biostrings::readDNAStringSet,
         fasta = FALSE, summary = FALSE, graphical = FALSE,
         save_regions = "none", not_found = FALSE,
         complement = FALSE, cpu = 1)[["positions"]],
    transform = cross(split_fasta, shard = !!(1:bigsplit), .id = FALSE)),
  
  itsxtrim = target(
    bind_rows(itsx_shard),
    transform = combine(itsx_shard, .by = Primer.pair)),
  
  positions = target(
    filter(join_derep$map, file == Trim.File) %>%
      left_join(itsxtrim, by = c("newmap" = "seq")),
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
    dada2::derepFastq(fls = infiles, n = 1e4)},
    transform = map(.data = !!meta3,
                    .id = PID)),
  
  err = target({
    err.fun <- ifelse(Tech == "PacBio",
                      PacBioErrfun,
                      loessErrfun)
    learnErrors(fls = derep2, nbases = 1e8,
                multithread=ignore(!!dadacores), randomize=TRUE,
                errorEstimationFunction = err.fun,
                verbose = TRUE)},
    transform = map(derep2, Tech, PID, .id = PID)),
  
  dada = target(
    dada(derep = derep2, err = err,
         multithread = ignore(!!dadacores), 
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
                              multithread = ignore(!!dadacores)),
    transform = map(seq_table, PID, .id = PID)),
  
  big_seq_table = target(
    join_seqs(nochim),
    transform = combine(nochim, .by = Region)),
  
  taxon = target(
    taxonomy(big_seq_table,
             reference = file_in(!!file.path(ref.dir, paste0(Reference, ".fasta.gz"))),
             multithread = ignore(!!dadacores)),
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

if (!interactive()) saveRDS(plan, "plan.rds")

#drake_plan_source(plan)
#predict_runtime(dconfig, jobs = ncpu)
if (interactive()) {
  dconfig <- drake_config(plan)
  vis_drake_graph(dconfig)
}
future::plan("multiprocess")
cat("\n Making pre-itsx targets (multiprocess)...\n")
# make embarassing targets at the beginning
make(plan,
     parallelism = "future",
     jobs = ncpu, jobs_preprocess = ncpu,
     retries = 2,
     keep_going = TRUE,
     caching = "worker",
     targets = str_subset(plan$target, "^split_fasta_")
)

# itsx (actually hmmer) does have an internal parallel option, but it isn't very efficient at
# using all the cores.  Instead, we divide the work into a large number of 
# shards and submit them all as seperate jobs on SLURM.
# failing that, do all the shards locally on the cores we have.
if (is_slurm) {
  cat("\n Making itsx_shard (SLURM)...\n")
  future::plan(batchtools_slurm, template = "slurm_itsx.tmpl",
               workers = sum(startsWith(plan$target, "itsx_shard")))
  make(plan,
       parallelism = "future",
       jobs = sum(startsWith(plan$target, "itsx_shard")),
       jobs_preprocess = ncpu,
       caching = "worker",
       targets = str_subset(plan$target, "^itsx_shard")
       )
  future::plan("multiprocess")
}

# embarrasing targets after itsx
cat("\n Making pre-dada targets (multiprocess)...\n")
make(plan,
     parallelism = "future",
     jobs = ncpu, jobs_preprocess = ncpu,
     retries = 2,
     keep_going = TRUE,
     caching = "worker",
     targets = c(str_subset(plan$target, "^derep2_"), "qstats_knit", "seq_counts")
)

# dada is internally parallel
cat("\n Making dada and taxonomy targets (loop)...\n")
make(plan,
     parallelism = "future",
     jobs = dadajobs, jobs_preprocess = ncpu,
     retries = 1,
     keep_going = TRUE,
     caching = "worker",
     targets = str_subset(plan$target, "^taxon_")
)
# finish up parallel
cat("\n Making all remaining targets (multiprocess)...\n")
make(plan,
     parallelism = "future",
     jobs = ncpu, jobs_preprocess = ncpu,
     retries = 2,
     keep_going = TRUE,
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
