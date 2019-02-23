library(magrittr)
library(tidyverse)
library(rlang)
library(glue)
library(readxl)

library(here)
library(drake)
library(future)

library(assertthat)
library(assertr)

library(FUNGuildR)
library(rITSx)
library(ShortRead)
library(dada2)

# This tells us the total number of cores on our current machine.

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
  rmd.dir <- here("writing")
  out.dir <- here("output")
  in.files <- list.files(path = trim.dir,
                         pattern = "\\.trim\\.fastq\\.gz$",
                         full.names = TRUE)
  dataset.file <- file.path(lab.dir, "datasets.csv")
  regions.file <- file.path(lab.dir, "regions.csv")
  platemap.file <- file.path(lab.dir, "Brendan_soil2.xlsx")
  platemap.sheet <- "Concentration samples"
} else {
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
  prereqs <- str_split(prereqs, " ") %>% unlist
#  cat("prereqs: ", prereqs, "\n")
  
  in.files <- str_subset(prereqs,
                         pattern = "\\.trim\\.fastq\\.gz$")
#  cat("in.files: ", in.files, "\n")
}

#### parallel setup ####

# how many cpus do we have on the local machine?
local_cpu <- detectCores()

# are we running slurm?
is_slurm <- nchar(Sys.which("sbatch")) > 0
is_local <- !is_slurm

# how many cpus can we get on one node?
if (is_slurm) {
  options(
    clustermq.scheduler = "slurm",
    clustermq.template = "slurm_clustermq.tmpl"
  )
  slurminfo <- system("sinfo -N --long", intern = TRUE)[-1]
  slurmcolnames <- strsplit(slurminfo[1], "\\s+")[[1]]
  slurminfo <- read_fwf(slurminfo,
                        fwf_empty(slurminfo, col_names = slurmcolnames),
                        skip = 1)
  max_cpu <- max(slurminfo[["CPUS"]])
  local_cpu <- as.integer(Sys.getenv("CORES_PER_TASK"))
} else {
  max_cpu <- max(local_cpu - 1, 1)
  local_cpu <- max_cpu
}
cat("max_cpu =", max_cpu, "\n")
cat("local_cpu =", local_cpu, "\n")

# the hmmersearch in itsx is very slow, but super parallel.
# we can break it into a lot of pieces and run it on several nodes.
bigsplit <- 80
# changing this value will invalidate most of the drake cache

# the computationally intensive parts of the dada pipeline have an internal
# parallel implementation.  This speeds them up, but with diminishing
# returns.  How many cores do we want?
dada_cpu <- 8
# reduce that if it's not possible
dada_cpu <- min(dada_cpu, max_cpu)
# what is the maximum number of dada jobs we can run?
# (we will later limit it to the number of regions*sequencing runs)
dadajobs <- Inf
# are we going to run the dada pipeline on this machine?
local_dada <- is_local # for now
# if we're runing on one machine, then we need to limit the number of jobs.
if (local_dada) {
  dadajobs <- max_cpu %/% dada_cpu
  # if this will leave us with extra cores, use them too.
  dada_cpu <- max_cpu %/% dadajobs
}

cat("dada_cpu =", dada_cpu, "\n")
cat("dadajobs =", dadajobs, "\n")

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

#### meta2 ####
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

#### meta3 ####
# meta3 has one row per region per plate
# The index is PID (Plate ID)
meta3 <- meta2 %>%
  select(Seq.Run, Plate, Region, PID) %>%
  unique() %>%
  mutate(region_derep = glue("region_derep_{PID}")) %>%
  left_join(datasets %>% select(-Regions)) %>%
  left_join(regions) %>%
  mutate_at(c("region_derep", "Region"), syms)
if (!interactive()) saveRDS(meta3, "meta3.rds")

#### meta4 ####
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

#### drake plan ####
cat("\nbuilding plan...\n")
tictoc::tic()
plan <- drake_plan(
  
  datasets = read_csv(file_in(!!dataset.file)),
  
  # dereplicate input sequences before running ITSx
  # this is done independently in parallel
  derep1 = target(
    tibble(file = Trim.File,
           derep = list(derepFastq(file_in(!!file.path(trim.dir, Trim.File)),
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
    transform = map(join_derep, Primer.pair,
                    .tag_in = step, .id = Primer.pair)),
  
  # find the location of LSU, 5.8S, and SSU (and this ITS1 and ITS2)
  # in the sequences from each chunk
  itsx_shard = target(
    itsx(in.file = split_fasta[[shard]],
         read_function = Biostrings::readDNAStringSet,
         fasta = FALSE, summary = FALSE, graphical = FALSE,
         save_regions = "none", not_found = FALSE,
         complement = FALSE, cpu = 1)[["positions"]],
    transform = cross(split_fasta,
                      shard = !!(1:bigsplit),
                      .tag_in = step,
                      .id = FALSE)),
  
  # combine the chunks into a master positions list.
  itsxtrim = target(
    bind_rows(itsx_shard),
    transform = combine(itsx_shard,
                        .tag_in = step,
                        .by = Primer.pair)),
  
  # find the entries in the positions list which correspond to the sequences
  # in each file
  positions = target(
    filter(join_derep$map, file == Trim.File) %>%
      left_join(itsxtrim, by = c("newmap" = "seq")),
    transform = map(.data = !!select(meta1, join_derep, itsxtrim, Trim.File, FID),
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
      out <- dada2::derepFastq(fls = infiles, n = 1e4,
                               qualityType = "FastqQuality")
      fqs <- FastqStreamer(infile, n = 1e4, qualityType = "FastqQuality")
      seqids <- character(0)
      while (length(fq <- yield(fqs))) {
        seqids <- c(seqids, as.character(fq@id))
      }
      close(fqs)
      names(out[["map"]]) <- seqids
    } else {
      out = NULL
    }
    list(Filter.File = out)
  },
  transform = map(.data = !!meta2,
                  .tag_in = step,
                  .id = RID)),
  
  region_derep = target(
    compact(c(derep2)),
    transform = combine(derep2, .by = PID, .tag_in = step),
  ),
  
  # Calibrate dada error models
  err = target({
    err.fun <- ifelse(Tech == "PacBio",
                      PacBioErrfun,
                      loessErrfun)
    learnErrors(fls = region_derep, nbases = 1e8,
                multithread = ignore(!!dada_cpu), randomize = TRUE,
                errorEstimationFunction = err.fun,
                verbose = TRUE,
                qualityType = "FastqQuality")},
    transform = map(.data = !!meta3,
                    .tag_in = step, .id = PID)),
  
  # Run dada denoising algorithm
  dada = target(
    dada(derep = region_derep, err = err,
         multithread = ignore(!!dada_cpu), 
         HOMOPOLYMER_GAP_PENALTY = eval(parse_expr(Homopolymer.Gap.Penalty)),
         BAND_SIZE = Band.Size,
         pool = eval(parse_expr(Pool))),
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
    makeSequenceTable(dada),
    transform = map(dada, .data = !!meta3, .tag_in = step, .id = PID)),
  
  # Remove likely bimeras
  nochim = target(
    dada2::removeBimeraDenovo(seq_table, method = "consensus",
                              multithread = ignore(!!dada_cpu)),
    transform = map(seq_table, .data = !!meta3, .tag_in = step, .id = PID)),
  
  # Join all the sequence tables for each region
  big_seq_table = target(
    dada2::mergeSequenceTables(tables = nochim),
    transform = combine(nochim, .tag_in = step, .by = Region)),
  
  # Assign taxonomy to each ASV
  taxon = target(
    taxonomy(big_seq_table,
             reference = file_in(!!file.path(ref.dir, paste0(Reference, ".fasta.gz"))),
             multithread = ignore(!!dada_cpu)),
    transform = map(.data = !!meta4, .tag_in = step, .id = TID)),
  
  # Download the FUNGuild database
  funguild_db = get_funguild_db(),
  
  # Assign ecological guilds to the ASVs based on taxonomy.
  guilds_table = target(
    funguild_assign(taxon, db = funguild_db),
    transform = map(taxon, TID, .tag_in = step, .id = TID)),
  
  platemap = read_platemap(file_in(!!platemap.file), platemap.sheet),
  
  # Count the sequences in each fastq file
  seq_count = target(
    tibble(
      file = file_in(fq.file),
      reads = system(glue("zcat {file} | grep -c '^@'", file = fq.file),
                     intern = TRUE) %>%
        as.integer),
    transform = map(fq.file = !!c(file.path(trim.dir, meta1$Trim.File),
                                 file.path(region.dir, meta2$Region.File),
                                 file.path(filter.dir, meta2$Filter.File)),
                    .tag_in = step)),
  
  # merge all the sequence counts into one data.frame
  seq_counts = target(
    bind_rows(seq_count),
    transform = combine(seq_count)),

  # calculate quality stats for the different regions
  qstats = target(
    q_stats(file_in(!!file.path(region.dir, Region.File))),
    transform = map(.data = !!meta2,
                    .tag_in = step,
                    .id = RID)),
  # join all the quality stats into one data.frame
  qstats_join = target(
    bind_rows(qstats) %>%
      tidyr::extract(
        col = "file",
        into = c("Seq.Run", "Plate", "Well", "Direction", "Region"),
        regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.trim\\.fastq\\.gz") %>%
      left_join(datasets),
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

plansteps <- unique(na.omit(plan[["step"]]))

if (!interactive()) saveRDS(plan, "plan.rds")

cat("\nCalculating outdated targets...\n")
tictoc::tic()
dconfig <- drake_config(plan, jobs_preprocess = local_cpu)
od <- outdated(dconfig)
tictoc::toc()

if (interactive()) {
  vis_drake_graph(dconfig,
                  # group = "step", clusters = plansteps,
                  targets_only = TRUE)
}
remove(dconfig)

if (is_slurm) {
  options(clustermq.scheduler = "slurm",
          clustermq.template = "slurm_clustermq.tmpl")
} else {
  options(clustermq.scheduler = "multicore")
}

timestamp <- strftime(Sys.time(), '%Y%m%d%H%M%S')

#### pre-ITSx ####
# make embarassingly parallel targets at the beginning
# if running on SLURM, use the clustermq backend.
# if running on the local machine, future has less overhead.
derep_targets <- str_subset(od, "^derep1_")
if (length(derep_targets)) {
  if (is_slurm) {
    derep_parallelism = "clustermq"
    # the widest part of the workflow is dereplication of each file,
    # so we can determine how many workers to use based on this.
    # These jobs are relatively quick though, so we don't need a single
    # worker each.
    derep_jobs = min(bigsplit, length(derep_targets) %/% 4)
    derep_template = list(log_file = glue("logs/derep-{timestamp}%a.log"))
    cat("\n Dereplicating input files (SLURM)...\n")
  } else {
    derep_parallelism = "future"
    future::plan(strategy = "multiprocess")
    derep_jobs = local_cpu
    derep_template <- list()
    cat("\n Dereplicating input files (multiprocess)...\n")
  }
  tictoc::tic()
  make(plan,
       parallelism = derep_parallelism,
       template = derep_template,
       jobs = derep_jobs,
       jobs_preprocess = local_cpu,
       retries = 2,
       elapsed = 3600, # 1 hour per target is more than enough
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = derep_targets
  )
  tictoc::toc()
  if (length(failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n All initial dereplication targets are up-to-date.\n")

#### Pre-ITSx targets ####
# These are computationally easy, but take a lot of memory, and would be
# inefficient to send to clustermq workers.  It's better to just do them locally.
preitsx_targets <- str_subset(od, "^split_fasta_")
if (length(preitsx_targets)) {
  cat("\nMaking targets to prepare for ITSx...\n")
  tictoc::tic()
  future::plan(strategy = "multiprocess")
  make(plan,
       parallelism = "future",
       jobs = max(local_cpu %/% dada_cpu, 1),
       jobs_preprocess = local_cpu,
       retries = 2,
       elapsed = 3600, # 1 hour
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = preitsx_targets
  )
  tictoc::toc()
  if (length(failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n All pre-itsx targets are up-to-date.\n")


#### ITSx ####
# itsx (actually hmmer) does have an internal parallel option, but it isn't
# very efficient at using all the cores.
# Instead, we divide the work into a large number of 
# shards and submit them all as seperate jobs on SLURM.
# failing that, do all the shards locally on the cores we have.
itsx_targets <- str_subset(od, "^itsx_shard")
if (is_slurm && length(itsx_targets)) {
  itsx_jobs <- min(max(ceiling(length(itsx_targets) / 2), bigsplit), length(itsx_targets))
  cat("\n Making itsx_shard (SLURM with ", itsx_jobs, "workers)...\n")
  tictoc::tic()
  make(plan,
       parallelism = "clustermq",
       template = list(log_file = glue("logs/itsx-{timestamp}%a.log"),
                       memory = 7*1024),
       jobs = itsx_jobs,
       elapsed = 3600*6, #6 hours
       jobs_preprocess = local_cpu,
       caching = "worker",
       cache_log_file = TRUE,
       targets = itsx_targets
  )
  tictoc::toc()
  if (length(failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
}

#### pre-DADA2 ####
# single-threaded targets after itsx
# for local runs, ITSx targets will also run here.
derep2_targets <- str_subset(od, "^filter_")
predada_targets <- c(str_subset(od, "^filter_"),
                     str_subset(od, "^(qstats_knit|seq_counts)$"))
if (length(predada_targets)) {
  
  predada_parallelism <- "future"
  future::plan(strategy = "multiprocess")
  predada_jobs <- local_cpu
  predada_template <- list()
  
  if (is_slurm 
      && nrow(meta3) > 2 * local_cpu) {
    predada_jobs <- nrow(meta3)
    predada_parallelism <- "clustermq"
    predada_template <- list(
      log_file = glue("logs/predada-{timestamp}%a.log"),
      memory = 7*1024) # 7 gb is 1 processor on a fat node or 2 on a regular
  }
  cat("\n Making pre-dada targets (multiprocess)...\n")
  tictoc::tic()
  make(plan,
       parallelism = predada_parallelism,
       template = predada_template,
       jobs = predada_jobs,
       jobs_preprocess = local_cpu,
       retries = 2,
       elapsed = 3600, #1 hour
       keep_going = FALSE,
       caching = "worker",
       cache_log_file = TRUE,
       targets = predada_targets
  )
  tictoc::toc()
  if (length(failed())) {
    if (interactive()) stop() else quit(status = 1)
  }
} else cat("\n Pre-DADA2 targets are up-to-date. \n")


#### DADA2 pipeline ####
# dada is internally parallel, so these need to be sent to nodes with multiple
# cores (and incidentally a lot of memory) or, if local, use all the cores of
# the local machine

dada_targets <- str_subset(od, "^taxon_")
if (length(dada_targets)) {
  if (is_slurm) {
    dada_parallelism <- "clustermq"
    dada_template <- list(ncpus = dada_cpu, 
                          log_file = glue("logs/dada-{timestamp)}%a.log"))
  } else {
    dada_parallelism <- if (dadajobs > 1) "future" else "loop"
    future::plan(strategy = "multiprocess")
    dada_template <- list()
  }
cat("\n Making DADA and taxonomy targets (multiprocess with", dada_cpu, "cores per target)...\n")
tictoc::tic()
make(plan,
     parallelism = dada_parallelism,
     jobs = dadajobs,
     jobs_preprocess = local_cpu,
     template = dada_template,
     retries = 1,
     elapsed = 3600*6, #6 hours
     keep_going = FALSE,
     caching = "worker",
     cache_log_file = TRUE,
     targets = dada_targets
)
tictoc::toc()
if (length(failed())) {
  if (interactive()) stop() else quit(status = 1)
}
} else cat("\n DADA2 pipeline targets are up-to-date.\n")

#### Finish ####
# For now the later steps are not very intensive, so they can be done
# using the resources of the master computer.
if (length(od)) {
cat("\n Making all remaining targets (multiprocess)...\n")
tictoc::tic()
future::plan(strategy = "multiprocess")
make(plan,
     parallelism = if (local_cpu > 1) "future" else "loop",
     jobs = local_cpu,
     jobs_preprocess = local_cpu,
     retries = 2,
     elapsed = 600, # 10 minutes
     keep_going = FALSE,
     caching = "worker",
     cache_log_file = TRUE
)
tictoc::toc()
if (length(failed())) {
  if (interactive()) stop() else quit(status = 1)
}
}
cat("\nAll targets are up-to-date.\n")


if (interactive()) {
  dconfig <- drake_config(plan)
  vis_drake_graph(dconfig)
}

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
