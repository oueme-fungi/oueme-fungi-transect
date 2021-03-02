##### drake plan for lightweight targets #####
# Targets in this plan don't require the resources of a computing cluster
# The whole plan should run in less than an hour on a modern desktop/laptop.
# This plan is run after drake_plan_heavy.R, which does the heavyweight targets,
# including region extraction, denoising, recombination, alignment,
# ML phylogeny, and primary taxonomic assignments.
# It performs statistical analyses on the results of the heavyweight targets,
# and also generates tables, figures, and the PDF of the manuscript.
#
# author Brendan Furneaux

#### Load config files ####
# the main config file will be sent from Snakemake, but we can otherwise load it directly
if (exists("snakemake")) {
  load(snakemake@input[["drakedata"]])
} else {
  load(here::here("data/plan/drake.Rdata"))
}

# ENTREZ key
# I'm not sure this is used anymore, since relevant code was exported to
# the reannotate project.
# It's not included in the git repository, because it is a personal key.
# Supply your own in the project root directory.
if (file.exists("ENTREZ_KEY")) Sys.setenv(ENTREZ_KEY = readLines("ENTREZ_KEY"))

# load packages quietly
suppressPackageStartupMessages({
  library(drake)
  library(magrittr)
  library(tidyverse)
  library(rlang)
  library(glue)
  library(drake)
  library(assertr)
  library(disk.frame)
  library(ape)
  library(ggplot2)
  library(futile.logger)
})

# these files define our datasets, rDNA regions, and plate layout
datasets <- read_csv(config$dataset, col_types = "cccccccicccccccccccicc")
regions <- read_csv(config$regions, col_types = "cccciiiiic")
platemap <- read_platemap(config$platemap, config$platemap_sheet)

# names of the LaTeX files and supporting directories required for manuscript
# submission
latexfiles <- c("transect_paper.tex", "transect_supplement.tex")
latexdirs <- file.path(c("transect_paper_files", "transect_supplement_files"), "figure-latex")

#### Load functions ####
# I've tried to split code that takes more than a few lines out of the main plan
# into one "functions_*.R" file per section of the plan.
for (f in list.files(config$rdir, "functions_.+.R", full.names = TRUE)) source(f)

# start logging
setup_log("drake_light")

#### meta plans ####
# Sometimes drake's plan transformations are not enough to correctly link up
# the dependencies between targets.
# In those cases, they are defined using these metaplans, via
# transform = map(.data = !![section]_meta, [...])

# Data frame for mapping between targets for building phyloseq experiment objects
physeq_meta <-
  tidyr::crossing(
    dplyr::select(datasets, "seq_run", "tech", "amplicon"),
    guild = c("fungi", "ecm"), #, "ecm2", "ecm3"),
    tibble::tibble(
      algorithm = c("PHYLOTAX", "Consensus"),
      taxset = c("hybrid", "short")
    )
  ) %>%
  dplyr::filter(
    seq_run != "is_057", # couldn't be demultiplexed
    amplicon == "Short" | taxset == "hybrid" # short consensus only on Short
  ) %>%
  mutate(
    fungi = glue::glue("fungi_{algorithm}"),
    ecm = glue::glue("ecm_{algorithm}"),
    ecm2 = glue::glue("ecm2_{algorithm}"),
    ecm3 = glue::glue("ecm3_{algorithm}"),
    taxa = glue::glue("{tolower(algorithm)}_hybrid")
  ) %>%
  mutate_at(vars(starts_with("ecm"), starts_with("fungi"), "taxa"), compose(syms, make.names))

#### Drake plan ####
plan2 <- drake_plan(
  #### Import "heavy" targets from the first plan ####
  file_meta = target(
    transform = map(primer_ID = !!unique(itsx_meta$primer_ID), .tag_in = step),
    trigger = trigger(mode = "blacklist")
  ),
  illumina_group = target(
    transform = map(
      seq_run = !!unique(illumina_meta$seq_run),
      .tag_in = step,
      .id = seq_run
    ),
    trigger = trigger(mode = "blacklist")
  ),
  allseqs = target(trigger = trigger(mode = "blacklist")),
  taxon = target(
    transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID),
    trigger = trigger(mode = "blacklist")
  ),
  raxml_decipher_unconst_long = target(trigger = trigger(mode = "blacklist")),
  big_seq_table = target(
    transform = map(region = !!(region_meta$region), .tag_in = step, .id = region),
    trigger = trigger(mode = "blacklist")
  ),
  chimeras = target(
    transform = map(.data = !!(dada_meta[,c("seq_run", "region")]), .tag_in = step, .id = c(seq_run, region)),
    trigger = trigger(mode = "blacklist")
  ),
  allchimeras = target(
    transform = map(region = !!(region_meta$region), .tag_in = step, .id = region),
    trigger = trigger(mode = "blacklist")
  ),
  qstats = target(trigger = trigger(mode = "blacklist")),

  #### Write additional data files ####
  # big_fasta
  # write the big_seq_table as a fasta files so that they can be clustered by
  # VSEARCH.
  big_fasta = target(
    write_big_fasta(big_seq_table,
                    file_out(big_fasta_file)),
    transform = map(.data = !!region_meta, .tag_in = step, .id = region)
  ),

  # Output Supplementary tables
  # write out tables
  # [tibble], creates file for supplement 2
  tagmap_long = platemap %>%
    filter(primer_pair == "ITS1_LR5") %>%
    select(year, site, buffer, x, sample_type, primer_pair, plate, well) %>%
    left_join(
      read_rds(file_in(!!file.path("config", "tags", "long.rds"))),
      by = "well"
    ) %>%
    select(year, site, buffer, x, sample_type, amplicon, everything()) %>%
    dplyr::mutate(
      sample_alias = dplyr::case_when(
        sample_type == "Sample" ~
          paste("OT", substr(as.character(year), 3,4), site, x, substr(buffer, 1, 1),
                substr(amplicon, 1, 1), sep = "-"),
        sample_type == "Pos" ~
          paste("OT-control", formatC(plate, width = 3, flag = "0"), substr(amplicon, 1, 1),
                sep = "-"),
        sample_type == "Blank" ~
          paste("OT-blank", formatC(plate, width = 3, flag = "0"), substr(amplicon, 1, 1),
                sep = "-")
      )
    )  %>%
    dplyr::left_join(read_accnos) %>%
    dplyr::select(-"Illumina-f", -"Illumina-r", -"IonTorrent") %T>%
    write_tsv(file_out(!!file.path("output", "supp_file_2.tsv")), na = ""),

  # [tibble], creates file for supplement 3
  tagmap_short = platemap %>%
    filter(primer_pair == "gITS7_ITS4") %>%
    select(year, site, buffer, x, sample_type, primer_pair, plate, well) %>%
    left_join(
      read_rds(file_in(!!file.path("config", "tags", "short.rds"))),
      by = "well"
    ) %>%
    dplyr::mutate(
      sample_alias = dplyr::case_when(
        sample_type == "Sample" ~
          paste("OT", substr(as.character(year), 3,4), site, x, substr(buffer, 1, 1),
                substr(amplicon, 1, 1), sep = "-"),
        sample_type == "Pos" ~
          paste("OT-control", formatC(plate, width = 3, flag = "0"), substr(amplicon, 1, 1),
                sep = "-"),
        sample_type == "Blank" ~
          paste("OT-blank", formatC(plate, width = 3, flag = "0"), substr(amplicon, 1, 1),
                sep = "-")
      )
    )  %>%
    dplyr::left_join(read_accnos) %>% {
      dplyr::left_join(
        dplyr::select(., -IonTorrent),
        dplyr::select(., well, IonTorrent) %>% dplyr::filter(!is.na(IonTorrent)),
        by = "well"
      )
    } %>%
    select(sample_alias, year, site, buffer, x, sample_type, amplicon, everything()) %T>%
    write_tsv(file_out(!!file.path("output", "supp_file_1.tsv")), na = ""),

  #### Phylogeny ####
  # Some of these targets depend on argets from the next (taxonomy) section

  # [phylo object]
  unrooted_tree = raxml_decipher_unconst_long$bipartitions,

  # Tulasnella is on a long branch because of a high divergence rate,
  # and its placement is highly unstable.
  # [character vector] tip labels (ASV hashes)
  tulasnella =
    dplyr::filter(taxon_table_primary, rank == "family", taxon == "Tulasnellaceae") %$%
    label %>%
    unique(),

  # find ASVs for which we can be confident of the kingdom assignment
  # at least 5 assignments with a confidence > 0.5
  # [character vector] tip labels (ASV hashes)
  sure_kingdoms =
    identify_taxa(taxon_table_hybrid, unrooted_tree,
                  rank = "kingdom", confidence = 0.5, n = 5,
                  ignore = tulasnella),

  # grab out the bikonts
  # in our dataset the confidently identified ones are all in Alveolata and
  # Viridiplantae, there are also some unconfident assignments for Rhizaria
  # and Stramenopila.
  # [a character vector] tip labels (ASV hashes)
  sure_bikonts = extract_bikonta(sure_kingdoms, unrooted_tree),

  # root the tree with the bikonts
  # [phylo object]
  rooted_tree = target(
    ape::root.phylo(
      unrooted_tree,
      sure_bikonts,
      resolve.root = TRUE,
      edgelabel = FALSE
    ),
    transform = map(amplicon = "Long", .tag_out = tree, .id = FALSE)
  ),

  # find the sequences that were confidently identified as fungi
  # [character vector]
  sure_fungi = dplyr::filter(sure_kingdoms, taxon == "Fungi")$label,

  # label the tree
  # [phylo object], also writes file for inspection (should go in Dryad)
  labeled_tree =
    phylotax::relabel_tree(
      tree = rooted_tree,
      old = taxon_labels$old,
      new = taxon_labels$new,
      drop = allchimeras_ITS2
    ) %T>%
    castor::write_tree(file_out(
      !!glue::glue("data/trees/labeled.tree", group = group)
    )),

  # [phylo object], also writes file for inspection (should go in Dryad)
  phylo_labeled_tree =
    phylotax::relabel_tree(
      tree = rooted_tree,
      old = phylotax_labels$old,
      new = phylotax_labels$new,
      drop = allchimeras_ITS2,
      quote = TRUE
    ) %T>%
    castor::write_tree(file_out(
      !!glue::glue("data/trees/labeled_phylotax.tree", group = group)
    )),

  #### Taxonomy ####

  # put together all the primary taxonomy assignments
  # [tibble] as returned by phylotax::taxtable
  taxon_table_primary = target(
    drake_combine(taxon) %>%
      combine_taxon_tables(allseqs) %>%
      add_warcup_kingdom(),
    transform = combine(taxon)
  ),

  # only assignments which came from the long amplicons
  # taxon_table_long = target(
  #   dplyr::filter(taxon_table_primary, region %in% c("LSU", "ITS")),
  #   transform = map(taxset = "long", .tag_out = taxon_table, .id = FALSE)
  # ),

  # only assignments which came from the short amplicons
  # [tibble] as returned by phylotax::taxtable
  taxon_table_short = target(
    dplyr::filter(taxon_table_primary, .data[["region"]] == "short"),
    transform = map(taxset = "short", region = "short", ref_region = "ITS",
                    .tag_out = taxon_table, .id = FALSE)
  ),

  # ASVs which were found in both long and short amplicon datasets were
  # identified twice by each of the ITS algorithm/reference combinations:
  # once as the full short amplicon, and once as the ITS region extracted
  # from the long amplicon.
  # Take the full-length ITS as likely to be the more accurate one
  # This version of the taxon table is used to make the best possible
  # assignments for the short reads.
  # [tibble] as returned by phylotax::taxtable
  taxon_table_hybrid = target(
    dplyr::group_by(taxon_table_primary, label) %>%
      dplyr::filter(!any(.data[["region"]] == "ITS") | .data[["region"]] != "short") %>%
      dplyr::ungroup(),
    transform = map(taxset = "hybrid", region = "All", ref_region = "All",
                    .tag_out = taxon_table, .id = FALSE)
  ),

  # Only ITS-based assignments. These are used for the metacoder trees
  # [tibble] as returned by phylotax::taxtable
  taxon_table_ITS = target(
    dplyr::filter(taxon_table_hybrid, .data[["region"]] %in% c("ITS", "short")),
    transform = map(taxset = "ITS", region = "All", ref_region = "ITS",
                    .tag_out = taxon_table, .id = FALSE)
  ),

  # Create labels for the tree(s) which show the assigned taxonomy before
  # refinement
  # [tibble] columns "old" and "new", giving old and new tip labels
  taxon_labels =
    phylotax::make_taxon_labels(taxon_table_hybrid, "n_reads", TRUE),

  # calculate (phylogenetic) consensus taxonomy
  # [phylotax object]
  phylotax_hybrid = target(
    phylotax::phylotax(
      tree = rooted_tree,
      taxa = taxon_table_hybrid,
      method = c(region = "All", reference = "All",
                 ref_region = "All", method = algorithm)
    ),
    transform = map(taxset = "hybrid", algorithm = "PHYLOTAX",
                    .tag_in = step, .tag_out = taxa, .id = FALSE)
  ),

  # [phylotax object]
  consensus = target(
    phylotax::phylotax(
      taxa = taxon_table,
      method = c(region = "All", reference = "All",
                 ref_region = "All", method = algorithm)
    ),
    transform = map(taxon_table, region, ref_region, algorithm = "Consensus",
                    .tag_in = step, .tag_out = taxa, .id = taxset)
  ),

  # make labels for phylogenetic consensus taxa
  # [tibble] columns "old" and "new", giving old and new tip labels
  phylotax_labels =
    phylotax_hybrid %$%
    bind_rows(assigned, retained) %>%
    phylotax::make_taxon_labels(),

  # Guild assignment ----

  # Load the FUNGuild database
  # Downloaded Feb 10, 2020
  # [tibble] as required by FUNGuildR::funguild_assign
  funguild_db = readRDS(file_in(!!file.path("reference", "funguild.rds"))),

  # [tibble] columns: name, region, reference, ref_region, method, label,
  # confidence, n_reads, kingdom:genus, Taxonomy, taxon, taxonomicLevel,
  # trophicMode, guild, growthForm, trait, confidenceRanking, notes,
  # citationSource, ECM
  guilds_PHYLOTAX = target(
    phylotax::extract_taxon(phylotax_hybrid, "Fungi", true_members = sure_fungi) %$%
      FUNGuildR::funguild_assign(widen_taxonomy(assigned), funguild_db) %>%
      dplyr::mutate(ECM = grepl("Ectomycorrhizal", guild)),
    transform = map(algorithm = "PHYLOTAX", .tag_out = guilds, .id = FALSE)
  ),

  # [tibble] columns: name, region, reference, ref_region, method, label,
  # confidence, n_reads, kingdom:genus, Taxonomy, taxon, taxonomicLevel,
  # trophicMode, guild, growthForm, trait, confidenceRanking, notes,
  # citationSource, ECM
  guilds_Consensus = target(
    phylotax::extract_taxon(consensus_short, "Fungi") %$%
      FUNGuildR::funguild_assign(widen_taxonomy(assigned), funguild_db) %>%
      dplyr::mutate(ECM = grepl("Ectomycorrhizal", guild)),
    transform = map(algorithm = "Consensus", .tag_out = guilds, .id = FALSE)
  ),

  #should be a no-op
  # [tibble] columns: name, region, reference, ref_region, method, label,
  # confidence, n_reads, kingdom:genus, Taxonomy, taxon, taxonomicLevel,
  # trophicMode, guild, growthForm, trait, confidenceRanking, notes,
  # citationSource
  fungi = target(
    dplyr::filter(guilds, kingdom == "Fungi"),
    transform = map(guilds, .tag_in = step, .id = algorithm)
  ),

  # [tibble] columns: name, region, reference, ref_region, method, label,
  # confidence, n_reads, kingdom:genus, Taxonomy, taxon, taxonomicLevel,
  # trophicMode, guild, growthForm, trait, confidenceRanking, notes,
  # citationSource
  # fungi_combined = target(
  #   unique(dplyr::bind_rows(fungi)),
  #   transform = combine(fungi)
  # ),

  # [tibble] columns: name, region, reference, ref_region, method, label,
  # confidence, n_reads, kingdom:genus, Taxonomy, taxon, taxonomicLevel,
  # trophicMode, guild, growthForm, trait, confidenceRanking, notes,
  # citationSource
  ecm = target(
    dplyr::filter(guilds, ECM),
    transform = map(guilds, .tag_in = step, .id = algorithm)
  ),

  # [tibble] columns: name, region, reference, ref_region, method, label,
  # confidence, n_reads, kingdom:genus, Taxonomy, taxon, taxonomicLevel,
  # trophicMode, guild, growthForm, trait, confidenceRanking, notes,
  # citationSource
  # ecm_combined = target(
  #   unique(dplyr::bind_rows(ecm)),
  #   transform = combine(ecm)
  # ),

  # [tibble] columns: name, region, reference, ref_region, method, label,
  # confidence, n_reads, kingdom:genus, Taxonomy, taxon, taxonomicLevel,
  # trophicMode, guild, growthForm, trait, confidenceRanking, notes,
  # citationSource
  ecm2 = target(
    dplyr::filter(ecm, !taxon %in% c("Peziza", "Pezizaceae", "Pyronemataceae")),
    transform = map(ecm, .tag_in = step, .id = algorithm)
  ),

  # [tibble] columns: name, region, reference, ref_region, method, label,
  # confidence, n_reads, kingdom:genus, Taxonomy, taxon, taxonomicLevel,
  # trophicMode, guild, growthForm, trait, confidenceRanking, notes,
  # citationSource
  ecm3 = target(
    dplyr::filter(ecm, confidenceRanking != "Possible"),
    transform = map(ecm, .tag_in = step, .id = algorithm)
  ),

  #### Phyloseq objects ####
  # see functions_phyloseq.R

  # [phyloseq::sample_data object]
  sample_data = assemble_sample_data(platemap, datasets),

  # [phyloseq::otu_table object]
  physeq_otu_table =
    assemble_otu_table(sample_data, relabel_seqtable(big_seq_table_ITS2), allchimeras_ITS2),

  # this phyloseq object contains all the data
  # [phyloseq::phyloseq object]
  proto_physeq = phyloseq::phyloseq(sample_data, physeq_otu_table),

  # trim down the phyloseq data to include only data from a single
  # sequencing technology, amplicon length, and guild (ECM vs all fungi) as
  # determined by taxonomic annotations from a particular algorithm.
  # [phyloseq::phyloseq object]
  physeq = target(
    build_physeq(taxa, physeq_otu_table, sample_data,
                 amplicon, tech, guild,
                 fungi, ecm, ecm2, ecm3),
    transform = map(.data = !!physeq_meta, .tag_in = step,
                    .id = c(guild, tech, amplicon, algorithm))
  ),

  #### Correlogram ####
  # see functions_mantel.R

  # [vegan::mantel.correlog object]
  correlog = target(
    correlog(physeq, metric, timelag),
    transform = cross(physeq, metric = c("bray", "wunifrac"), timelag = c(0, 1),
                      .tag_in = step,
                      .id = c(guild, metric, tech, amplicon, algorithm, timelag))
  ),

  #### Spatiotemportal variogram ####
  # see functions_variogram.R

  # empirical spatial variogram
  # [tibble] columns: dist, gamma, bin
  # each observation gets its own row.
  variog = target(
    variog(physeq, metric, breaks = c(1:25 - 0.5, 31000)),
    transform = cross(physeq, metric = c("bray", "wunifrac"),
                      .tag_in = step,
                      .id = c(guild, metric, tech, amplicon, algorithm))
  ),

  # [list]
  #   pars: [list]
  #     nugget [numeric scalar]
  #     sill   [numeric scalar]
  #     range  [numeric scalar]
  #   summary: [tibble] columns: term, estimate, std.error, statistic, p.value,
  #                              conf.low, conf.high
  #   predict: [tibble] columns: dist, timelag, gamma
  variofit2 = target(
    fit_variogram(variog),
    transform = map(variog, .tag_in = step,
                    .id = c(guild, metric, tech, amplicon, algorithm))
  ),

  # empirical soatiotemporal variogram
  # [tibble] columns: dist, gamma, timelag, bin
  # each observation gets its own row.
  variogST = target(
    variogST(physeq, metric, breaks = c(1:25 - 0.5, 30000)),
    transform = cross(physeq, metric = c("bray", "wunifrac"),
                      .tag_in = step,
                      .id = c(guild, metric, tech, amplicon, algorithm))
  ),

  # [list]
  #   pars: [list]
  #     timerange [numeric scalar]
  #     nugget    [numeric scalar]
  #     sill      [numeric scalar]
  #     range     [numeric scalar]
  #   summary: [tibble] columns: term, estimate, std.error, statistic, p.value,
  #                              conf.low, conf.high
  #   predict: [tibble] columns: dist, timelag, gamma
  variofitST2 = target(
    fit_variogramST(variogST, variofit2),
    transform = map(variogST, variofit2, .tag_in = step, .id = c(guild, metric, tech, amplicon, algorithm))
  ),

  ##### Tables and Figures #####################################################
  # Targets below intended to compile data in useful formats
  # for making the tables, figures, and in-line results in the paper.
  # see functions_figures.R

  # draw the rDNA, primers, and amplicons (Fig1a)
  amplicon_plot = draw_amplicon_plot(file_in("writing/rDNA_amplicons.csv")),

  # list of all sequences by region and sequencing run, along with number of
  # reads and sequence length
  # [tibble] columns: seq_run, region, seq, reads, length
  reads_table = target(
    compile_reads_table(big_seq_table),
    transform = combine(big_seq_table)
  ),

  # Summary information for each region
  # [tibble] columns: seq_run, region, ASVs, length_min, length_q1, length_mean,
  #                   length_med, length_q3, length_max, reads
  region_table = reads_table %>%
    dplyr::group_by(seq_run, region) %>%
    dplyr::summarize(
      ASVs = dplyr::n_distinct(seq),
      length_min = min(length),
      length_q1 = reldist::wtd.quantile(length, 0.25, weight = reads),
      length_med = reldist::wtd.quantile(length, 0.5, weight = reads),
      length_mean = weighted.mean(length, w = reads),
      length_q3 = reldist::wtd.quantile(length, 0.75, weight = reads),
      length_max = max(length),
      reads = sum(reads)
    ),

  # plot of ASV richness vs. denoising success vs. length
  length_richness_plot = draw_length_richness_plot(region_table, readcounts),

  # Table of reads per ASV for each sequencing run
  # [tibble] columns: seq, is_057, pb_483, pb_500, SH-2257
  asv_table = compile_asv_table(big_seq_table_ITS2),

  # Table of reads per OTU for each sequencing run
  # [tibble] columns: seq, is_057, pb_483, pb_500, SH-2257
  otu_table =
   compile_otu_table(file_in("data/clusters/ITS2.table"), big_fasta_ITS2),

  # Find the central sequence for each 97% OTU
  # [tibble] identity, ASVseq, OTUseq
  # ASVseq and OTUseq are ASV sequence hashes.
  otu_map = generate_otu_map(file_in(here::here("data/clusters/ITS2.uc"))),

  #### Alpha diversity ####
  # functions used in this section are in richness.R

  # for ASVs and OTUs, make inputs for iNEXT
  # the format is a list, where names are the sample names,
  # and values are integer vectors giving the number of reads
  # for all ASVs/OTUs with nonzero reads in that sample
  # [list]
  pre_iNEXT_ASV = target(
    nonzero_asv_counts(big_seq_table_ITS2, platemap),
    transform = map(cluster = "ASV", .tag_out = pre_iNEXT, .id = FALSE)
  ),

  # [list]; see pre_iNEXT_ASV
  pre_iNEXT_OTU = target(
    nonzero_otu_counts(file_in(here::here("data/clusters/ITS2.table")), platemap),
    transform = map(cluster = "OTU", .tag_out = pre_iNEXT, .id = FALSE)
  ),

  # [list of data.frame] one data.frame per sample.
  #    columns:
  #     m: number of samples
  #     method: interpolated, observed, or extrapolated
  #     order; order of Hill number; 0 for richness
  #     qD: estimated Hill number
  #     qD.LCL: lower confidence limit for qD
  #     qD.UCL: upper confidence limit for qD
  #     SC: estimated sample coverage
  #     SC.LCL: lower confidence interval for SC
  #     SC.UCL: upper confidence interval for SC
  sample_iNEXT = target(
    lapply(pre_iNEXT, function(x) iNEXT:::iNEXT.Ind(
      x,
      m = unique(c(round(10^seq(0, log10(sum(x)), 0.1)), sum(x))),
      q = 0)),
    transform = map(pre_iNEXT, .id = cluster),
    dynamic = map(pre_iNEXT)
  ),

  # rarefied species (ASV/OTU) richness
  # iNEXT::estimateD causes problems when there is only one species
  # it also does more calculation than we actually need (i.e. Shannon and Simpson)
  # iNEXT:::iNEXT.Ind causes problem when there is only one value of m
  # [list] rarefied richness for each sample
  sample_rarefy = target(
    lapply(pre_iNEXT, iNEXT:::iNEXT.Ind, m = c(1,100), q = 0) %>%
      vapply(magrittr::extract, 0, 2, "qD"),
    transform = map(pre_iNEXT, .id = cluster),
    dynamic = map(pre_iNEXT)
  ),

  # matched richness values for different strategies in the same sample well
  # [tibble] columns: well, strategy_x, value_x, strategy_y, value_y
  rarefy_data = target(
    reshape_rarefy_data(readd(sample_rarefy)),
    transform = map(sample_rarefy, .id = cluster)
  ),

  # [tibble] strategy_x, strategy_y, slope, lower 0.95, upper 0.95, r_squared,
  #          slope_label, r2_label
  demingfits = target(
    fit_deming_slope(rarefy_data),
    transform = map(rarefy_data, .id = cluster)
  ),

  # plot of corresponding richness values for sifferent strategies.
  # For ASVs as a main figure, for OTUS as supplementary
  alpha_plot = target(
    make_alpha_plot(rarefy_data, demingfits, cluster),
    transform = map(demingfits, rarefy_data, cluster, .id = cluster)
  ),

  # Species (OTU/ASV) accumulation plot
  # ASVs as main figure, OTUs in supplement
  accum_plot = target(
    make_accum_plot(readd(sample_iNEXT), cluster),
    transform = map(sample_iNEXT, cluster, .id = cluster)
  ),

  #### Quality stats ####
  # functions used in this section in qstats.R

  # parse the sample names in the qstats
  # [tibble] columns: seq_run, plate, well, read, step, stat, value, nreads, region
  # (also a column for each stat, but these are there by accident)
  parsed_qstat = parse_qstat(qstats),

  demuxlength =  parsed_qstat %>%
    filter(is.na(read) | read != "_R2", stat == "length") %>%
    filter(is.na(step) | step == "demux", !is.na(well), ) %>%
    group_by(seq_run) %>%
    summarize(value = reldist::wtd.quantile(value, weight = nreads)) %>%
    deframe(),

  demuxqual = parsed_qstat %>%
    filter(is.na(read) | read != "_R2", stat == "erate") %>%
    filter(is.na(step) | step == "demux", !is.na(well)) %>%
    group_by(seq_run) %>%
    summarize(value = round(weighted.mean(value, w = nreads, na.rm = TRUE), 4)) %>%
    deframe(),

  readcounts = filter(parsed_qstat, is.na(read) | read != "_R2", stat == "n"),

  rawcounts = readcounts %>%
    filter(is.na(step) | step == "raw", is.na(well)) %>%
    group_by(seq_run) %>%
    summarize(
      # nreads =  prettyNum(sum(nreads), big.mark = ",")
      nreads = sum(nreads)
    ) %>%
    deframe(),

  demuxcounts = readcounts %>%
    filter(is.na(step) | step == "demux", !is.na(well)) %>%
    group_by(seq_run) %>%
    summarize(
      # nreads =  prettyNum(sum(nreads), big.mark = ",")
      nreads = sum(nreads)
    ) %>%
    deframe(),

  filtercounts_full = readcounts %>%
    filter(step == "filter", is.na(region) | region %in% c("long", "short")) %>%
    group_by(seq_run, region) %>%
    summarize(
      # nreads =  prettyNum(sum(nreads), big.mark = ",")
      nreads = sum(nreads)
    ) %>%
    select(seq_run, nreads) %>%
    deframe(),

  filtercounts_ITS2 = readcounts %>%
    filter(step == "filter", region == "ITS2") %>%
    group_by(seq_run, region) %>%
    summarize(
      # nreads =  prettyNum(sum(nreads), big.mark = ",")
      nreads = sum(nreads)
    ) %>%
    select(seq_run, nreads) %>%
    deframe(),

  regioncounts = readcounts %>%
    filter(step == "lsux", region == "5_8S") %>%
    group_by(seq_run) %>%
    summarize(
      # nreads =  prettyNum(sum(nreads), big.mark = ",")
      nreads = sum(nreads)
    ) %>%
    deframe(),

  asvcounts = count_asvs(asv_table, allseqs),

  bioinf_table = compile_bioinf_table(
    rawcounts, demuxcounts, filtercounts_full, regioncounts, filtercounts_ITS2,
    region_table, asv_table, allseqs, datasets
  ),

  bioinf_header =
    tibble(name = names(bioinf_table)) %>%
    separate(name, into = c("tech", "amplicon", "type"), sep = "_") %>%
    mutate(
      tech = paste(tech, plyr::mapvalues(tech, datasets$tech, datasets$machine, FALSE))
    ),

  ##### Venn diagrams ####
  # functions used in this section are in venn.R

  venn_ASV = venndata(asv_table, ASVs, cols = c("pb_500", "pb_483", "SH-2257", "is_057")),

  venn_OTU = venndata(otu_table, OTUs, cols = c("pb_500", "pb_483", "SH-2257", "is_057")),

  vennplot_tech_ASV = vennplot_tech(asv_table, ASVs, letter = "a"),

  vennplot_tech_OTU = vennplot_tech(otu_table, OTUs, letter = "b"),

  vennplot_amplicon_ASV = vennplot_amplicon(asv_table, ASVs, letter = "c"),

  vennplot_amplicon_OTU = vennplot_amplicon(otu_table, OTUs, letter = "d"),

  ##### Read abundance comparisons ####
  big_table = {
    out <- bind_rows(
      mutate(asv_table, type = "ASV"),
      mutate(otu_table, type = "OTU")
    ) %>%
      select(-seq) %>%
      select(type, everything())

    names(out)[-1] <-
      tibble(seq_run = names(out)[-1]) %>%
      left_join(datasets, by = "seq_run") %>%
      glue::glue_data("{tech} {machine} - {amplicon}")
    out
  },

  multi_table =
    datasets %>%
    mutate_at("tech", factor, levels = c("PacBio", "Illumina", "Ion Torrent"), ordered = TRUE) %>%
    mutate_at("amplicon", factor, levels = c("Long", "Short"), ordered = TRUE) %>%
    arrange(tech, amplicon) %>%
    mutate(var = paste(tech, machine) %>%
             factor(ordered = TRUE, levels = unique(.))) %>% {
               crossing(
                 select(., x_var = var, x_amplicon = amplicon),
                 select(., y_var = var, y_amplicon = amplicon),
                 type = c("ASV", "OTU")
               )
             } %>%
    filter(
      x_var < y_var | (x_var == y_var & x_amplicon < y_amplicon),
      paste(x_var, x_amplicon, sep = " - ") %in% names(big_table),
      paste(y_var, y_amplicon, sep = " - ") %in% names(big_table)
    ) %>%
    group_by(x_var, x_amplicon, y_var, y_amplicon, type) %>%
    group_map(choosevars, .data = big_table) %>%
    bind_rows() %>%
    group_by(x_var, x_amplicon, y_var, y_amplicon, type),

  comparisons = multi_table %>%
    group_modify(
      function(d, g) {
        filter_all(d, ~. > 0) %>%
          mutate_all(log10) %>%
          {bind_cols(
            lm(y ~ x, .) %>%
              broom::glance(.),
            lm(y ~ offset(x), .) %>%
              broom::tidy(.) %>%
              select(term, estimate) %>%
              pivot_wider(names_from = term, values_from = estimate)
          )}
      }
    ) %>%
    mutate(
      label = paste("R^2 ==", formatC(r.squared, format = "f", digits = 2))
    ),

  ##### Taxonomy chart ####
  # see functions_taxonomy.R

  # fraction of reads for each ASV in each sequencing run
  taxon_reads_proto = make_taxon_reads_proto(proto_physeq, datasets),

  taxon_reads_primary =
    compile_primary_taxon_reads(taxon_reads_proto, taxon_table_primary) %>%
    format_taxon_reads(funguild_db),

  # PHYLOTAX for long amplicon PacBio
  # region = NA, reference = NA, ref_region = NA; method = "PHYLOTAX"
  taxon_reads_long_phylotax =
    select_taxon_reads(phylotax_hybrid, taxon_reads_proto,
                       seq_run == "pb_500") %>%
    format_taxon_reads(funguild_db),

  # LCA consensus;
  # region = "All", reference = "All", ref_region = "All"; method = "Consensus"
  # Short and long amplicons;
  taxon_reads_hybrid_consensus =
    select_taxon_reads(consensus_hybrid, taxon_reads_proto,
                       region = "hybrid") %>%
    format_taxon_reads(funguild_db),

  # LCA consensus calculated only on short amplicon region;
  # Apply only to Short amplicons
  # region = "short", reference = "All", ref_region = "ITS"; method = "Consensus"
  taxon_reads_short_consensus =
    select_taxon_reads(consensus_short, taxon_reads_proto, amplicon == "Short",
                       region = "short") %>%
    format_taxon_reads(funguild_db),

  # Phylotax if available; otherwise LCA consensus
  # Apply only to Short amplicons
  # region = "short", reference = "All", ref_region = "All", method = "combined"
  taxon_reads_hybrid_phylotax =
    select_taxon_reads(phylotax_hybrid, taxon_reads_proto, amplicon == "Short",
                       region = "hybrid") %>%
    format_taxon_reads(funguild_db),

  # LCA consensus using only ITS references
  # region = "All", reference = "All", ref_region = "ITS", method = "Consensus"
  taxon_reads_ITS_consensus =
    select_taxon_reads(consensus_ITS, taxon_reads_proto, region = "ITS") %>%
    format_taxon_reads(funguild_db),

  taxon_reads = bind_rows(
      taxon_reads_primary,
      taxon_reads_long_phylotax,
      taxon_reads_hybrid_consensus,
      taxon_reads_short_consensus,
      taxon_reads_hybrid_phylotax
    ),

  # Fraction of reads and ASVs in each dataset which were assigned using the
  # tree
  phylotax_long_reads = group_by(taxon_reads_long_phylotax, tech, amplicon) %>%
    mutate(ASVs = 1, ASV_frac = 1/n()) %>%
    filter(label %in% phylotax_hybrid$tree$tip.label) %>%
    summarize_at(c("reads", "ASVs", "ASV_frac"), sum),
  phylotax_hybrid_reads = group_by(taxon_reads_hybrid_phylotax, tech, amplicon) %>%
    mutate(ASVs = 1, ASV_frac = 1/n()) %>%
    filter(label %in% phylotax_hybrid$tree$tip.label) %>%
    summarize_at(c("reads", "ASVs", "ASV_frac"), sum),


  taxon_chart = compile_taxon_chart(taxon_reads),

  # show success of taxonomic assignment for different datasets
  taxon_chart_plot = draw_taxon_chart(taxon_chart),

  # show class composition of fungal community in different datasets
  taxon_plot_fungi_class = draw_fungal_classes(taxon_reads, datasets),

  # show kingdoms for reads with an extra-short ITS2
  taxon_plot_short_kingdom =
    draw_short_kingdoms(big_seq_table_ITS2, taxon_reads, datasets),

  # show ECM status assignments for different datasets
  taxon_plot_ecm = draw_ecm_assignment(taxon_reads, datasets),

  #### Stats for text ####
  agaricus =
    taxon_table_primary %>%
    group_by(label, rank) %>%
    filter(
      rank == "genus",
      any(taxon == "Agaricus"),
      n_distinct(taxon, na.rm = TRUE) == 1
    ) %$%
    unique(label),

  agaricus_reads =
    pos_control_physeq %>%
    phyloseq::prune_taxa(taxa = agaricus) %>%
    phyloseq::prune_taxa(taxa = colSums(phyloseq::otu_table(.)) > 0) %>%
    phyloseq::otu_table(),

  pos_control_physeq =
    phyloseq::subset_samples(
      proto_physeq,
      sample_type == "Pos" |
        (seq_run == "is_057" & well == "E7")
    ),

  neg_control_physeq = phyloseq::subset_samples(
    proto_physeq,
    sample_type == "Blank"
  ),

  nonpos_control_physeq =
    phyloseq::subset_samples(
      proto_physeq,
      sample_type != "Pos" &
        !(seq_run == "is_057" & well == "E7")
    ),

  pos_control =
    agaricus_reads %>%
    as("matrix") %>%
    as_tibble(rownames = "sample") %>%
    pivot_longer(-1, names_to = "seq", values_to = "reads") %>%
    filter(reads > 0) %>%
    extract(
      sample,
      into = c("seq_run", "plate", "well"),
      regex = "([a-zA-Z]{2}[-_]\\d{3,4})(00\\d)([A-H]1?\\d)"
    ) %>%
    filter(seq_run != "is_057") %$%
    unique(seq),

  pos_control_reads =
    nonpos_control_physeq %>%
    phyloseq::prune_taxa(taxa = pos_control) %>%
    phyloseq::otu_table() %>%
    as.matrix() %>%
    rowSums(),

  all_reads =
    nonpos_control_physeq %>%
    phyloseq::otu_table() %>%
    as.matrix() %>%
    rowSums(),

  pos_control_data =
    nonpos_control_physeq %>%
    phyloseq::sample_data() %>%
    as("data.frame") %>%
    mutate(
      pc_reads = pos_control_reads,
      all_reads = all_reads,
      pc_frac = pc_reads/all_reads
    ) %>%
    select(seq_run, well, plate, pc_frac, pc_reads, all_reads),

  pc_pc_reads =
    pos_control_physeq %>%
    phyloseq::prune_taxa(taxa = pos_control) %>%
    phyloseq::otu_table() %>%
    as.matrix() %>%
    rowSums(),

  pc_all_reads =
    pos_control_physeq %>%
    phyloseq::otu_table() %>%
    as.matrix() %>%
    rowSums(),

  pc_pos_control_data =
    pos_control_physeq %>%
    phyloseq::sample_data() %>%
    as("data.frame") %>%
    mutate(
      pc_reads = pc_pc_reads,
      all_reads = pc_all_reads,
      nonpc_reads = all_reads - pc_reads,
      nonpc_frac = nonpc_reads/all_reads
    ) %>%
    select(seq_run, well, plate, nonpc_frac, nonpc_reads, all_reads),

  neg_control_reads =
    neg_control_physeq %>%
    phyloseq::otu_table() %>%
    as.matrix() %>%
    rowSums(),

  neg_control_data =
    neg_control_physeq %>%
    phyloseq::sample_data() %>%
    as("data.frame") %>%
    mutate(
      all_reads = neg_control_reads
    ) %>%
    select(seq_run, well, plate, all_reads),

  agaricus_fasta =
    allseqs %>%
    filter(hash %in% phyloseq::taxa_names(agaricus_reads)) %$%
    set_names(short, hash) %>%
    discard(is.na) %>%
    Biostrings::RNAStringSet() %>%
    Biostrings::writeXStringSet(file_out("temp/agaricus.fasta")),

  agaricus_order =
    phyloseq::merge_samples(proto_physeq, "seq_run") %>%
    phyloseq::otu_table() %>%
    apply(MARGIN = 1, order, decreasing = TRUE) %>%
    apply(MARGIN = 2, match, x = 1),

  #### Variogram plots ####
  variog_data = target({
    rlang::quo(list(correlog, variog, variogST, variofit2, variofitST2))
    ignore(plan2) %>%
      filter(step == "correlog") %>%
      select("correlog", timelag, guild, metric, tech, amplicon, algorithm) %>%
      mutate_at(
        "correlog",
        lapply,
        readd,
        character_only = TRUE,
        path = ignore(cache_dir)
      ) %>%
      mutate_at("correlog", map, "mantel.res") %>%
      mutate_at("correlog", map_lgl, ~ any(.x[,"Pr(corrected)"] < 0.05 & .x[, "Mantel.cor"] > 0, na.rm = TRUE)) %>%
      pivot_wider(names_from = "timelag", values_from = "correlog", names_prefix = "mantel_") %>%
      left_join(
        ignore(plan2) %>%
          filter(step == "variofit2") %>%
          select("variofit2", guild, metric, tech, amplicon, algorithm),
        by = c("guild", "metric", "tech", "amplicon", "algorithm")
      ) %>%
      left_join(
        ignore(plan2) %>%
          filter(step == "variofitST2") %>%
          select("variofitST2", "variogST", guild, metric, tech, amplicon, algorithm),
        by = c("guild", "metric", "tech", "amplicon", "algorithm")
      ) %>%
      mutate(
        variofit = ifelse(mantel_1, .data[["variofitST2"]], .data[["variofit2"]])
      ) %>%
      mutate_at(
        c("variofit", "variogST", "variofitST2", "variofit2"),
        map,
        readd,
        character_only = TRUE,
        path = ignore(cache_dir)
      ) %>%
      select(mantel_0, mantel_1, "variofitST2", "variofit2", "variofit",
             "variogST", guild, metric, tech, amplicon, algorithm) %>%
      mutate_at(c("guild", "metric", "tech", "amplicon", "algorithm"),
                str_replace_all, "\"", "") %>%
      mutate(
        guild = plyr::mapvalues(guild, c("fungi", "ecm"), c("All Fungi", "ECM")),
        metric = plyr::mapvalues(metric, c("bray", "wunifrac"), c("Bray-Curtis", "W-UniFrac")),
        amplicon = factor(
          amplicon,
          levels = c("Short", "Long")
        ),
        algorithm = factor(
          algorithm,
          levels = c("Consensus", "PHYLOTAX", "combined"),
          labels = c("Cons", "PHYLO", "PHYLO")
        )
      )
  },
  transform = combine(correlog, variog, variogST, variofit2, variofitST2)
  ),

  variog_points = unnest(variog_data, "variogST") %>%
    mutate_at("timelag", replace_na, "0") %>%
    filter(dist <= 25) %>%
    mutate_at("timelag", factor),

  variog_empir = variog_points %>%
    # mutate(bin = cut(dist, c(1:13 - 0.5, 14.5, 16.5, 19.5, 24.5))) %>%
    group_by(guild, metric, tech, amplicon, algorithm, dist, timelag) %>%
    summarize(
      gamma = mean(gamma),
      # dist = mean(dist),
      n = n()
    ),

  variog_confint =
    variog_points %>%
    group_by(guild, metric, tech, amplicon, algorithm, dist, timelag) %>%
    summarize(
      lower = quantile(gamma, 0.025),
      upper = quantile(gamma, 0.975)
    ),

  variog_fit = variog_data %>%
    mutate_at("variofit", purrr::map, "predict") %>%
    unnest("variofit") %>%
    filter(mantel_1 | timelag == "0", mantel_0) %>%
    mutate_at("timelag", factor),

  variofit_table =
    variog_data %>%
    select(mantel_0, mantel_1, variofit_0 = variofit2, variofit_1 = variofitST2,
           metric, amplicon, tech, guild, algorithm) %>%
    pivot_longer(1:4, names_to = c(".value", "timelag"), names_sep = "_") %>%
    mutate(params = map(variofit, "summary")) %>%
    unnest(params) %>%
    filter(term %in% c("range", "timerange")),

  #### Heat trees ####
  # see functions_heat_trees.R

  taxdata = generate_taxdata(taxon_reads_ITS_consensus, asv_table),

  taxdata_ECM = generate_ecm_taxdata(taxon_reads_ITS_consensus, asv_table),

  heattree = draw_full_heattree(taxdata, file_out("temp/heattree.pdf")),

  heattree_length_compare =
    draw_amplicons_heattree(taxdata, file_out("temp/heattree_amplicons.pdf")),

  ecm_heattree =
    draw_ecm_heattree(taxdata_ECM, file_out("temp/ecm_heattree.pdf")),

  buffer_compare_taxmap = {
    taxmap <- metacoder::parse_phyloseq(buffer_compare_physeq)
    taxmap$data$otu_table <- metacoder::calc_obs_props(taxmap, "otu_table")
    taxmap$data$read_count <- metacoder::calc_taxon_abund(taxmap, "otu_table")
    taxmap
  },

  buffer_compare_long =  draw_buffer_compare_heattree(
    buffer_compare_taxmap,
    amplicon = "Long",
    file_out("temp/compare_long.pdf")
  ),

  buffer_compare_short = draw_buffer_compare_heattree(
    buffer_compare_taxmap,
    amplicon = "Short",
    file_out("temp/compare_short.pdf")
  ),

  #### Supplementary figures ####
  # functions for targets in this section are in functions_figures.R

  phylotax_fig = draw_phylotax_figure(),

  length_table = compile_length_table(asv_table, allseqs, datasets),

  full_dens_short = length_density_plot(length_table, tech ~ .,
                                        amplicon == "Short", region == "short"),

  full_dens_long = length_density_plot(length_table, tech ~ .,
                                       amplicon == "Long", region == "long"),

  full_ecdf_short = full_length_ecdf_plot(length_table, "Short"),

  full_ecdf_long = full_length_ecdf_plot(length_table, "Long", "Length (bp)"),

  its2_dens = length_density_plot(length_table, tech + amplicon ~ .,
                                  region == "ITS2"),

  its2_ecdf = its2_length_ecdf_plot(length_table),

  region_length_fig = draw_region_length_figure(reads_table, regions),

  buffer_compare_physeq = make_buffer_compare_physeq(proto_physeq, taxon_reads),

  # Generate ASV table for comparisons between buffer solutions
  buffer_asv_physeq = generate_buffer_asv_physeq(proto_physeq, fungi_PHYLOTAX$label),

  # Create distance matrix for the buffer solution comparison
  buffer_asv_bc_dist =
    # Normalize reads to 1
    phyloseq::transform_sample_counts(buffer_asv_physeq, function(x) x / sum(x)) %>%
    # Bray-Curtis distance
    phyloseq::distance("bray"),

  # Calculate PERMANOVA
  buffer_adonis = target(
    vegan::adonis2(
      buffer_asv_bc_dist ~ paste(site, x) + paste(tech, amplicon) + buffer + year,
      data = as(phyloseq::sample_data(buffer_asv_physeq), "data.frame"),
      by = "margin",
      permutations = 9999
    ) %>%
      broom::tidy() %>%
      column_to_rownames("term"),
    seed = 12345
  ),

  # Calculate PCoA for different buffer samples
  buffer_asv_pcoa =
    vegan::capscale(
      buffer_asv_bc_dist ~ 1,
      data = as(phyloseq::sample_data(buffer_asv_physeq), "data.frame")
    ),

  # Generate ASV table for comparisons between sequencing methods
  tech_asv_physeq = generate_tech_asv_physeq(proto_physeq, fungi_PHYLOTAX$label),

  # Add taxonomy to the ASV table and cluster at the class level
  tech_class_physeq =
    generate_tech_class_physeq(tech_asv_physeq, fungi_PHYLOTAX),

  # Sample data for easy reference
  tech_class_sample_data =
    as(phyloseq::sample_data(tech_class_physeq), "data.frame"),

  # Calculate Bray-Curtis distance for comparing technologies
  tech_class_bc_dist = phyloseq::distance(tech_class_physeq, "bray"),

  # PERMANOVA test for technologies
  tech_class_adonis = target(
    vegan::adonis2(
      tech_class_bc_dist ~ paste(site, x, year) + tech + amplicon,
      data = tech_class_sample_data,
      by = "margin",
      permutations = 9999
    ) %>%
      broom::tidy() %>%
      column_to_rownames("term"),
    seed = 12345
  ),

  # PCoA for technologies
  tech_class_pcoa =
    vegan::capscale(
      phyloseq::otu_table(tech_class_physeq) ~ Condition(paste(site, x, year)),
      data = tech_class_sample_data,
      distance = "bray"
    ),

  # Species (actually class) scores for PCoA
  tech_class_pcoa_scores =
    vegan::scores(tech_class_pcoa, display = "species") %>%
    as_tibble(rownames = "class") %>%
    # Take only scores which are at least 10% of the maximum
    mutate(L = sqrt(MDS1^2 + MDS2^2)) %>%
    filter(L >= max(L) / 10) %>%
    # Abbreviate using the first two letters
    mutate(abbrev = substr(class, start = 1, stop = 2)),

  # PCoA plot
  tech_class_pcoa_plot = draw_tech_class_pcoa_plot(
    pcoa = tech_class_pcoa,
    sample_data = tech_class_sample_data,
    scores = tech_class_pcoa_scores
  ),

  # Generate ASV table for comparisons between sequencing methods for ECM
  tech_ecm_asv_physeq =
    tech_asv_physeq %>%
    # Only ECM
    phyloseq::prune_taxa(taxa = ecm_PHYLOTAX$label),

  # Add taxonomy to the ECM ASV table and cluster at the family level
  tech_ecm_fam_physeq =
    generate_tech_ecm_fam_physeq(tech_asv_physeq, ecm_PHYLOTAX),

  # Sample data for easy reference
  tech_ecm_fam_sample_data =
    as(phyloseq::sample_data(tech_ecm_fam_physeq), "data.frame"),

  # Calculate Bray-Curtis distance for comparing technologies
  tech_ecm_fam_bc_dist = phyloseq::distance(tech_ecm_fam_physeq, "bray"),

  # PERMANOVA test for technologies
  tech_ecm_fam_adonis = target(
    vegan::adonis2(
      tech_ecm_fam_bc_dist ~ paste(site, x, year) + tech + amplicon,
      data = tech_ecm_fam_sample_data,
      by = "margin",
      permutations = 9999
    ) %>%
      broom::tidy() %>%
      column_to_rownames("term"),
    seed = 12345
  ),

  # PCoA for technologies
  tech_ecm_fam_pcoa =
    vegan::capscale(
      phyloseq::otu_table(tech_ecm_fam_physeq) ~ Condition(paste(site, x, year)),
      data = tech_ecm_fam_sample_data,
      distance = "bray"
    ),

  # Species (actually class) scores for PCoA
  tech_ecm_fam_pcoa_scores =
    vegan::scores(tech_ecm_fam_pcoa, display = "species") %>%
    as_tibble(rownames = "family") %>%
    # Take only scores which are at least 10% of the maximum
    mutate(L = sqrt(MDS1^2 + MDS2^2)) %>%
    filter(L >= max(L) / 10) %>%
    # Abbreviate using the first two letters
    mutate(abbrev = substr(family, start = 1, stop = 2)),

  # PCoA plot
  tech_ecm_fam_pcoa_plot =
    left_join(
      vegan::scores(tech_ecm_fam_pcoa, display = "sites") %>%
        as_tibble(rownames = "sample"),
      tech_ecm_fam_sample_data %>%
        as_tibble(rownames = "sample"),
      by = "sample"
    ) %>%
    mutate(
      group = paste(tech, "–", amplicon) %>%
        factor(
          levels = c("Illumina – Short", "PacBio – Short", "PacBio – Long"),
          ordered = TRUE
        )
    ) %>%
    ggplot(aes(x = MDS1, y = MDS2, color = group, pch = group)) +
    geom_vline(xintercept = 0, color = "gray80") +
    geom_hline(yintercept = 0, color = "gray80") +
    geom_point() +
    geom_text(
      aes(x = MDS1 * 3, y = MDS2 * 3, label = abbrev),
      data = tech_ecm_fam_pcoa_scores,
      inherit.aes = FALSE
    ) +
    scale_color_discrete(name = NULL) +
    scale_discrete_manual(aesthetics = "pch", name = NULL, values = 1:3) +
    coord_equal() +
    xlab("PCoA1") +
    ylab("PCoA2"),

  length_kw =
    reads_table %>%
    filter(region == "short") %>%
    group_by(length, seq_run) %>%
    summarize(reads = sum(reads)) %>%
    ungroup() %>%
    slice(rep(1:n(), reads)) %>%
    mutate_at("seq_run", factor) %$%
    kruskal.test(length, g = seq_run),

  fungi_read =
    taxdata %>%
    taxa::filter_taxa(
      taxon_names == "Fungi",
      reassign_obs = FALSE,
      reassign_taxa = FALSE
    ) %>%
    extract2("data") %>%
    extract2("tax_read") %>%
    set_names(names(.) %>%
                str_extract("(PacBio|Illumina|Ion.Torrent)_(Long|Short)")) %>%
    unlist(),

  fungi_asv =
    taxdata %>%
    taxa::filter_taxa(
      taxon_names == "Fungi",
      reassign_obs = FALSE,
      reassign_taxa = FALSE
    ) %>%
    extract2("data") %>%
    extract2("tax_asv") %>%
    set_names(names(.) %>%
                str_extract("(PacBio|Illumina|Ion.Torrent)_(Long|Short)")) %>%
    unlist(),

  ecm_counts =
    taxon_reads %>%
    filter(Algorithm == "PHYLOTAX") %>%
    group_by(seq_run) %>%
    mutate(ASVs = 1/n()) %>%
    filter(ECM %in% c("Possible ECM", "Probable ECM", "Highly Probable ECM")) %>%
    summarize_at(c("ASVs", "reads"), sum),

  ecm_read = select(ecm_counts, seq_run, reads) %>% deframe(),
  ecm_asv = select(ecm_counts, seq_run, ASVs) %>% deframe(),

  short_otu_table = {
    big_fasta_short
    read_tsv(
      file_in("data/clusters/short.table"),
      col_types = cols(
        .default = col_integer(),
        `#OTU ID` = col_character()
      )
    ) %>%
      column_to_rownames("#OTU ID") %>%
      as.matrix() %>%
      t() %>%
      as_tibble(rownames = "sample") %>%
      tidyr::extract(
        col = "sample",
        into = c("seq_run", "plate", "well"),
        regex = "([a-zA-Z]{2}[-_]\\d{3,4})_(\\d{3})([A-H]1?\\d)"
      ) %>%
      tidyr::gather(key = "seq", value = "reads", -(1:3)) %>%
      dplyr::group_by(seq_run, seq) %>%
      dplyr::summarize(reads = sum(reads)) %>%
      dplyr::left_join(
        tibble(
          seq = colnames(big_seq_table_short),
          length = nchar(seq)
        ) %>%
          mutate_at("seq", chartr, old = "T", new = "U") %>%
          mutate_at("seq", tzara::seqhash),
        by = "seq"
      )
  },

  # Plot of reads along transect
  reads_plot = parsed_qstat %>%
    dplyr::filter(stat == "n") %>%
    tidyr::extract(col = "well", into = c("row", "col"),
                   regex = "([A-H])(\\d+)") %>%
    dplyr::filter(
      is.na(step) | step == "demux",
      !is.na(seq_run),
      read != "_R2" | is.na(read)
    ) %>%
    dplyr::group_by(seq_run, plate, row, col) %>%
    dplyr::summarize(raw_reads = sum(nreads)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_at("col", as.integer) %>%
    dplyr::left_join(
      proto_physeq %>% {
        dplyr::mutate(
          phyloseq::sample_data(.),
          reads = phyloseq::otu_table(.) %>% rowSums()
        )
      } %>%
        dplyr::select(seq_run, plate, row, column, site, year, qual, sample_type,
                      x, buffer, tech, amplicon, reads) %>%
        dplyr::mutate(row = LETTERS[row]) %>%
        unique(),
      by = c("seq_run", "plate", "row", "col" = "column")
    ) %>%
    dplyr::filter(tech != "Ion Torrent", !is.na(year), !is.na(seq_run),
                  !is.na(x), !is.na(buffer)) %>%
    tidyr::complete(tidyr::nesting(seq_run, tech, amplicon), site, x,
                    tidyr::nesting(year, buffer)) %>%
    tidyr::pivot_longer(
      cols = c("reads", "raw_reads"),
      names_to = "type",
      values_to = "reads"
    ) %>%
    dplyr::mutate(
      type = factor(type, levels = c("raw_reads", "reads")),
      buffer = factor(buffer, levels = c("Xpedition", "LifeGuard")),
      tech = factor(tech, levels = c("PacBio", "Illumina")),
      "Sample group" = paste(year, buffer) %>%
        ifelse(type == "raw_reads", "Raw", .) %>%
        forcats::fct_relevel("Raw")
    ) %>%
    dplyr::arrange(`Sample group`) %>%
    ggplot(aes(x = x, y = reads, group = `Sample group`, fill = `Sample group`)) +
    geom_col(position = "identity") +
    ggnomics::facet_nested(
      cols = vars(tech, amplicon),
      rows = vars(year, buffer, site),
      nest_line = TRUE,
      resect = unit(1, "mm")
    ) +
    scale_fill_manual(
      values = c(
        RColorBrewer::brewer.pal(3, "Set2") %>%
          set_names(c("2015 Xpedition", "2016 LifeGuard", "2016 Xpedition")),
        Raw = "gray20"
      ),
      breaks = c("2015 Xpedition", "2016 LifeGuard", "2016 Xpedition", "Raw"),
      labels = c("2015 Xpedition", "2016 LifeGuard", "2016 Xpedition", "Filtered")) +
    scale_y_log10(breaks = c(10, 1e3, 1e5),
                  labels = c("10", "1k", "100k"), position = "right") +
    xlab("x (m)") +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      panel.spacing = unit(1, "mm"),
      plot.margin = margin(l = 0.5, unit = "mm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    ) +
    geom_hline(yintercept = 100, linetype = 2, alpha = 0.5),

  conc_plot =
    proto_physeq %>%
    phyloseq::sample_data(.) %>%
    as("data.frame") %>%
    select(plate, row, column, site, year, sample_type, x, buffer, amplicon, dna_conc, pcr_conc) %>%
    arrange(plate, row, column, amplicon) %>%
    unique() %>%
    mutate_at("dna_conc", as.double) %>%
    pivot_longer(cols = c("dna_conc", "pcr_conc"), names_to = "step", values_to = "conc", names_pattern = "([[:alpha:]]+)_conc") %>%
    filter(!is.na(year)) %>%
    mutate(
      amplicon = ifelse(step == "dna", "", amplicon),
      buffer = factor(buffer, levels = c("Xpedition", "LifeGuard")),
      step = factor(step, levels = c("dna", "pcr"), labels = c("DNA", "PCR product")),
      conc = 1000 * conc,
      "Sample group" = paste(year, buffer)
    ) %>%
    unique() %>%
    complete(nesting(amplicon, step), site, x, nesting(year, buffer, `Sample group`)) %>%
    arrange(`Sample group`) %>%
    ggplot(aes(x = x, y = conc, group = `Sample group`, fill = `Sample group`)) +
    geom_col() +
    ggnomics::facet_nested(
      cols = vars(step, amplicon),
      rows = vars(year, buffer, site),
      nest_line = TRUE,
      resect = unit(1, "mm"),
      switch = "y"
    ) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set2") %>% set_names(c("2015 Xpedition", "2016 LifeGuard", "2016 Xpedition")), guide = "none") +
    scale_y_log10(
      breaks = c(10, 1e3, 1e5),
      labels = c("0.1", "1", "100"),
      name = "DNA concentration (ng/µL)"
    ) +
    xlab("x (m)") +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "mm"),
      plot.margin = margin(r = 0.5, unit = "mm"),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    ),

  short_length_table =
    reads_table %>%
    filter(region == "short") %>%
    group_by(seq_run, length) %>%
    summarize(reads = sum(reads)) %>%
    tidyr::complete(seq_run, length, fill = list(reads = 0)),

  # summary stats for quality check pdf, and supplementary plot
  parsed_qstat_raw = parsed_qstat %>%
    dplyr::group_by(seq_run, region, step, stat, value, region, read) %>%
    dplyr::summarize_at("nreads", sum) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_at(
      "region",
      factor,
      levels = c(NA, "ITS2", "ITS1", "5_8S", "LSU1", "LSU2", "LSU3", "LSU4",
                 "D1", "D2", "D3", "short", "ITS", "LSU", "32S", "long")
    ) %>%
    dplyr::left_join(
      dplyr::select(datasets, dataset, seq_run, tech, machine, amplicon),
      by = "seq_run"
    ) %>%
    dplyr::mutate_at(
      "dataset", factor,
      levels = c("short-ion", "short-illumina", "short-pacbio", "long-pacbio")
    ) %>%
    tidyr::replace_na(list(read = "")) %>%
    dplyr::mutate(
      `Read type` = paste(tech, amplicon, substr(read, 2, 3)) %>%
        trimws()
    ),

  ## Error rate -- strategies

  erate_plot = dplyr::filter(parsed_qstat_raw, stat == "erate", step == "raw") %>%
    dplyr::arrange(value) %>%
    dplyr::group_by(`Read type`) %>%
    dplyr::mutate(
      ecdf = cumsum(nreads)/sum(nreads),
      label = ifelse(cumsum(value >= 1e-5) == 1, as.character(`Read type`), ""),
      label = gsub("(Illumina|Ion Torrent) Short", "\\1", label)
    ) %>%
    dplyr::filter(value >= 1e-5) %>%
    ggplot(aes(
      x = value,
      y = ecdf,
      group = `Read type`,
      color = `Read type`,
      label = label
    )) +
    geom_line() +
    ggrepel::geom_text_repel(
      # nudge_x = 0.2,
      # xlim = c(3, NA),
      ylim = c(0.05, NA),
      min.segment.length = 0,
      size = 3,
      segment.alpha = 0.3
      # hjust = 0
    ) +
    scale_x_continuous(
      name = "Expected error-rate",
      trans = reverselog_trans(10),
      limits = c(NA, 3e-6)
    ) +
    scale_y_continuous(name = "Fraction passing") +
    scale_color_read(guide = "none"),

  ## Expected number of errors --strategies

  eexp_plot = dplyr::filter(parsed_qstat_raw, stat == "eexp", step == "raw") %>%
    dplyr::arrange(value) %>%
    dplyr::group_by(`Read type`) %>%
    dplyr::mutate(
      ecdf = cumsum(nreads)/sum(nreads),
      label = ifelse(cumsum(value >= 0.1) == 1, as.character(`Read type`), ""),
      label = gsub("(Illumina|Ion Torrent) Short", "\\1", label)
    ) %>%
    dplyr::filter(value >= 0.1) %>%
    ggplot(aes(
      x = value,
      y = ecdf,
      group = `Read type`,
      color = `Read type`,
      label = label
    )) +
    geom_vline(xintercept = 3, linetype = "dashed", color = "gray50") +
    geom_line() +
    ggrepel::geom_text_repel(
      nudge_x = 0.2,
      xlim = c(1, NA),
      min.segment.length = 0,
      size = 3,
      segment.alpha = 0.3,
      hjust = 0) +
    scale_x_continuous(
      name = "Expected number of errors",
      trans = reverselog_trans(10),
      limits = c(NA, 0.03)
    ) +
    scale_y_continuous(name = "Fraction passing") +
    scale_color_read(guide = "none"),

  ## Expected number of errors - regions
  eexp_region_plot =
    dplyr::filter(parsed_qstat_raw, stat == "eexp", step == "lsux", amplicon == "Long") %>%
    dplyr::arrange(value) %>%
    dplyr::mutate_at("region", forcats::fct_relevel, "ITS", "ITS1", "32S", "5_8S",
                     "ITS2", "LSU", "LSU1", "D1",
                     "LSU2", "D2", "LSU3", "D3", "LSU4", "long") %>%
    dplyr::group_by(region) %>%
    dplyr::mutate(
      ecdf = cumsum(nreads)/sum(nreads),
      `Region type` = dplyr::recode(
        region,
        ITS1 = "extracted",
        ITS2 = "extracted",
        `5_8S` = "extracted",
        LSU1 = "extracted",
        LSU2 = "extracted",
        LSU3 = "extracted",
        LSU4 = "extracted",
        D1 = "extracted",
        D2 = "extracted",
        D3 = "extracted",
        ITS = "concatenated",
        LSU = "concatenated",
        `32S` = "concatenated",
        long = "full read"
      ),
      label = ifelse(cumsum(value >= 0.1) == 1, as.character(region), ""),
      label = gsub("_", ".", label)
    ) %>%
    dplyr::filter(value >= 0.1) %>%
    ggplot(aes(
      x = value,
      y = ecdf,
      linetype = `Region type`,
      group = region,
      color = region,
      label = label
    )) +
    geom_vline(xintercept = 3, linetype = "dashed", color = "gray50") +
    geom_line() +
    scale_x_continuous(
      name = "Expected number of errors",
      trans = reverselog_trans(10),
      limits = c(NA, 0.03)
    ) +
    ggrepel::geom_text_repel(
      nudge_x = 0.2,
      xlim = c(1, NA),
      size = 3,
      min.segment.length = 0,
      segment.alpha = 0.3,
      hjust = 0
    ) +
    # directlabels::geom_dl(method = "last.qp", size = 1) +
    scale_y_continuous(name = "Fraction passing") +
    scale_color_discrete(guide = "none") +
    scale_linetype_manual(
      values = c(extracted = "dashed", concatenated = "dotdash", `full read` = "solid")
    ) +
    theme(
      legend.position = c(0.01, 0.01),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.key = element_blank()
    ),

  length_model_data =
    reads_table %>%
    filter(region == "short") %>%
    mutate_at("seq", tzara::seqhash) %>%
    mutate_at("seq", factor) %>%
    tidyr::complete(seq_run, nesting(seq, length), fill = list(reads = 0)),

  ## Render paper ----
  # write *.aux and *.tex files.  Remove extra bibliographies from *.tex files
  supplement_tex =
    withr::with_dir(
      "writing",
      {
        file_in(!!file.path("writing", "all.bib"))
        knitr_in(!!file.path("writing", "transect_supplement.Rmd"))
        file_in(!!file.path("writing", "preamble_supplement.tex"))
        file_in(!!file.path("scripts", "functions_output.R"))
        rmarkdown::render(
          "transect_supplement.Rmd",
          clean = FALSE
        )
        file_out(!!file.path("writing", "transect_supplement.tex"))
        readLines("transect_supplement.tex") %>%
          magrittr::extract(!grepl("\\\\printbibliography", .)) %>%
          writeLines("transect_supplement.tex")
      }
    ),

  article_tex = withr::with_dir(
    "writing",
    {
      knitr_in(!!file.path("writing", "transect_paper.Rmd"))
      file_in(!!file.path("writing", "preamble.tex"))
      file_in(!!file.path("scripts", "functions_output.R"))
      rmarkdown::render(
        "transect_paper.Rmd",
        clean = FALSE
      )
      file_out(!!file.path("writing", "transect_paper.tex"))

      readLines("transect_paper.tex") %>%
        magrittr::extract(!grepl("\\\\printbibliography", .) | !duplicated(.)) %>%
        writeLines("transect_paper.tex")
    }
  ),

  # use difftex to mark changes from the original submission.
  supplement_diff_tex = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_supplement_original.tex"))
      file_in(!!file.path("writing", "transect_supplement_rev1.tex"))
      file_out(!!file.path("writing", "transect_supplement_diff.tex"))
      tinytex::tlmgr_install(setdiff(c("latexdiff", "listings"), tinytex::tl_pkgs()))
      system2(
        command = "latexdiff",
        args = list(
          "--type=CFONT",
          "--add-to-config",
          '"FLOATENV=sidewaysfigure"',
          "--graphics-markup=2",
          "transect_supplement_original.tex",
          "transect_supplement_rev1.tex"
        ),
        stdout = TRUE
      ) %>%
        stringr::str_replace_all(
          "transect_paper",
          "transect_paper_diff"
        ) %>%
        writeLines("transect_supplement_diff.tex")
    }
  ),

  supplement_diff2_tex = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_supplement_rev1.tex"))
      file_in(!!file.path("writing", "transect_supplement.tex"))
      file_out(!!file.path("writing", "transect_supplement_diff2.tex"))
      tinytex::tlmgr_install(setdiff(c("latexdiff", "listings"), tinytex::tl_pkgs()))
      system2(
        command = "latexdiff",
        args = list(
          "--type=CFONT",
          "--add-to-config",
          '"FLOATENV=sidewaysfigure"',
          "--graphics-markup=2",
          "transect_supplement_rev1.tex",
          "transect_supplement.tex"
        ),
        stdout = TRUE
      ) %>%
        stringr::str_replace_all(
          "transect_paper",
          "transect_paper_diff2"
        ) %>%
        writeLines("transect_supplement_diff2.tex")
    }
  ),

  article_diff_tex = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_paper_original.tex"))
      file_in(!!file.path("writing", "transect_paper_rev1.tex"))
      file_out(!!file.path("writing", "transect_paper_diff.tex"))
      tinytex::tlmgr_install(setdiff(c("latexdiff", "listings"), tinytex::tl_pkgs()))
      system2(
        command = "latexdiff",
        args = list(
          "--type=CFONT",
          "transect_paper_original.tex",
          "transect_paper_rev1.tex"
        ),
        stdout = TRUE
      ) %>%
        stringr::str_replace_all(
          "transect_supplement",
          "transect_supplement_diff"
        ) %>%
        writeLines("transect_paper_diff.tex")
    }
  ),

  article_diff2_tex = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_paper_rev1.tex"))
      file_in(!!file.path("writing", "transect_paper.tex"))
      file_out(!!file.path("writing", "transect_paper_diff2.tex"))
      tinytex::tlmgr_install(setdiff(c("latexdiff", "listings"), tinytex::tl_pkgs()))
      system2(
        command = "latexdiff",
        args = list(
          "--type=CFONT",
          "transect_paper_rev1.tex",
          "transect_paper.tex"
        ),
        stdout = TRUE
      ) %>%
        stringr::str_replace_all(
          "transect_supplement",
          "transect_supplement_diff2"
        ) %>%
        writeLines("transect_paper_diff2.tex")
    }
  ),

  # render the latex to PDFs
  supplement_pdf = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_supplement.tex"))
      file_in(!!file.path("writing", "transect_paper.tex"))
      file_out(!!file.path("writing", "transect_supplement.pdf"))
      tinytex::xelatex(
        "transect_supplement.tex",
        bib_engine = "biber",
        clean = FALSE
      )
    }
  ),

  # supplement_diff_pdf = withr::with_dir(
  #   "writing",
  #   {
  #     file_in(!!file.path("writing", "transect_supplement_diff.tex"))
  #     file_in(!!file.path("writing", "transect_paper_diff.tex"))
  #     file_out(!!file.path("writing", "transect_supplement_diff.pdf"))
  #     tinytex::xelatex(
  #       "transect_supplement_diff.tex",
  #       bib_engine = "biber",
  #       clean = FALSE
  #     )
  #   }
  # ),

  supplement_diff2_pdf = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_supplement_diff2.tex"))
      file_in(!!file.path("writing", "transect_paper_diff2.tex"))
      file_out(!!file.path("writing", "transect_supplement_diff2.pdf"))
      tinytex::xelatex(
        "transect_supplement_diff2.tex",
        bib_engine = "biber",
        clean = FALSE
      )
    }
  ),

  article_pdf = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_paper.tex"))
      file_in(!!file.path("writing", "transect_supplement.tex"))
      file_in(!!file.path("writing", "transect_supplement.pdf"))
      file_out(!!file.path("writing", "transect_paper.pdf"))
      tinytex::xelatex(
        "transect_paper.tex",
        bib_engine = "biber",
        clean = FALSE
      )
    }
  ),

  article_diff_pdf = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_paper_diff.tex"))
      file_in(!!file.path("writing", "transect_supplement_diff.tex"))
      file_in(!!file.path("writing", "transect_supplement_diff.pdf"))
      file_out(!!file.path("writing", "transect_paper_diff.pdf"))
      tinytex::xelatex(
        "transect_paper_diff.tex",
        bib_engine = "biber",
        clean = FALSE
      )
    }
  ),

  article_diff2_pdf = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_paper_diff2.tex"))
      file_in(!!file.path("writing", "transect_supplement_diff2.tex"))
      file_in(!!file.path("writing", "transect_supplement_diff2.pdf"))
      file_out(!!file.path("writing", "transect_paper_diff2.pdf"))
      tinytex::xelatex(
        "transect_paper_diff2.tex",
        bib_engine = "biber",
        clean = FALSE
      )
    }
  ),

  # zip up the latex source files for submission
  latex_dir =
    withr::with_dir(
      "writing",
      {
        file_in(!!file.path("writing", c(latexfiles, latexdirs)))
        file_out(!!file.path("output", "latex.tar.gz"))
        tar(
          tarfile = "latex.tar.gz",
          files = c(
            latexfiles,
            list.files(path = latexdirs, full.names = TRUE)
          ),
          compression = "gzip",
          tar = "tar"
        )
        file.rename("latex.tar.gz", file.path("..", "output", "latex.tar.gz"))
      }
    ),

  # word count
  wordcount = {
    if (!"texcount" %in% tinytex::tl_pkgs()) tinytex::tlmgr_install("texcount")
    system2(
      list.files(
        file.path(tinytex::tinytex_root(), "bin"),
        pattern = "texcount",
        recursive = TRUE,
        full.names = TRUE
      ),
      args = file_in(!!file.path("writing", "transect_paper.tex")),
      stdout = file_out(!!file.path("writing", "transect_paper.wc"))
    )
  },

  # take only the fungi from the PHYLOTAX and Consensus calls
  phylotax_fungi = phylotax::extract_taxon(
    phylotax = phylotax_hybrid,
    taxon = "Fungi",
    true_members = sure_fungi
  ),

  #### Output tree ####
  tree_fig = tree_figure(
    phylotax = phylotax_fungi,
    outfile = file_out(!!file.path("output", "supp_file_3.pdf")),
    charscale = 0.015,
    rankgap = 0.015
  ),

  #### prepare submission materials for ENA ####
  ENA_samples = target(
    compile_samples(sample_data),
    transform = map(
      compile_samples = list(
        compile_soil_samples,
        compile_control_samples,
        compile_blank_samples
      ),
      sample_type = c("soil", "control", "blank"),
      .id = sample_type
    )
  ),

  write_samples = target(
    write_ena_samples(
      ENA_samples,
      file_in(!!paste0("config/ENA_", sample_type, "_sample_template.tsv")),
      file_out(!!paste0("output/ENA_", sample_type, "_samples.tsv"))
    ),
    transform = map(ENA_samples, sample_type, .id = sample_type)
  ),

  pacbio_ion_manifest = target(
    generate_pacbio_ion_manifest(
      dplyr::bind_rows(ENA_samples),
      file_meta_gits7its4,
      file_meta_its1lr5
    ),
    transform = combine(ENA_samples)
  ),

  illumina_manifest = target(
    generate_illumina_manifest(
      dplyr::bind_rows(ENA_samples),
      illumina_group_SH.2257
    ),
    transform = combine(ENA_samples)
  ),

  read_accnos = lookup_submitted_accnos(path = file_in("output/ENA/reads")),
  pretaxid = target(
    .fun(phylotax_hybrid),
    transform = map(
      .fun = list(group_by_taxon, group_by_taxon_env),
      .lookupfun = list(lookup_ncbi_taxon, lookup_ena_taxon),
      taxon_type = c("ind", "env"),
      .id = taxon_type
      )
  ),

  taxid = target(
    do.call(.lookupfun, pretaxid),
    transform = map(pretaxid, .lookupfun, .id = taxon_type),
    dynamic = map(pretaxid),
    retries = 2
  ),

  taxid_all = target(
    {
      taxids <- bind_rows(readd(taxid))
      write_csv(taxids, file_out(!!sprintf("output/%s_taxa.csv", taxon_type)))
      dplyr::select(taxids, label, targetTaxonomy, targetRank) %>%
        write_csv(file_out(!!sprintf("output/%s_taxa_target.csv", taxon_type)))
      dplyr::select(taxids, label, Taxonomy, rank) %>%
        write_csv(file_out(!!sprintf("output/%s_taxa_ena.csv", taxon_type)))
      taxids
    },
    transform = map(taxid, taxon_type, .id = taxon_type)
  ),
  # transform = FALSE,
  trace = TRUE
) %>%
  # bind_plans(plan2, .) %>%
  # transform_plan() %>%
  dplyr::filter(ifelse(amplicon == '"Short"', metric != '"wunifrac"', TRUE) %|% TRUE)

saveRDS(plan2, "data/plan/drake_light.rds")
njobs <- max(local_cpus() %/% 2, 1)
parallelism <- if (njobs > 1) "clustermq" else "loop"
options(clustermq.scheduler = "multicore")
#if (interactive()) {
cache_dir <- ".drake"
cache <- drake_cache(cache_dir)
dconfig <- drake_config(plan2, cache = cache)
vis_drake_graph(dconfig, targets_only = TRUE,
                group = "step",
                clusters = c("correlog", "variog", "variogST", "variofit2",
                             "variofitST2", "physeq", "big_seq_table",
                             "big_fasta", "taxon", "chimeras", "allchimeras"))
make(
  plan2,
  cache = cache,
  parallelism = parallelism,
  jobs = njobs,
  # targets = c("article_pdf", "supplement_pdf", "article_diff_pdf", "supplement_diff_pdf", "wordcount"),
  # targets = c("phylotax_long_reads", "phylotax_hybrid_reads", "bioinf_table"),
  # keep_going = TRUE,
  lazy_load = FALSE,
  memory_strategy = "autoclean",
  garbage_collection = TRUE,
  cache_log_file = "drake_light_cache.csv",
  console_log_file = file.path(normalizePath("logs"), "drake_light.log")
)
