if (exists("snakemake")) {
  load(snakemake@input[["drakedata"]])
} else {
  load("data/plan/drake.Rdata")
}
if (file.exists("ENTREZ_KEY")) Sys.setenv(ENTREZ_KEY = readLines("ENTREZ_KEY"))
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

setup_log("drake_light")

source(file.path(config$rdir, "dada.R"))
source(file.path(config$rdir, "mantel.R"))
source(file.path(config$rdir, "variogram.R"))
source(file.path(config$rdir, "qstats.R"))
source(file.path(config$rdir, "taxonomy.R"))
source(file.path(config$rdir, "plate_check.R"))
source(file.path(config$rdir, "output_functions.R"))

datasets <- read_csv(config$dataset, col_types = "cccccccicccccccccccicc")
regions <- read_csv(config$regions, col_types = "cccciiiiic")

# Data frame for mapping between targets for guild assignment
guilds_meta <- tibble(
  preguild_taxa = c(
    "preguild_taxa_long_fungi_PHYLOTAX",
    "preguild_taxa_short_fungi_Consensus",
    "preguild_taxa_short_fungi_PHYLOTAX.Cons"
  ),
  amplicon = c("Long", "Short", "Short"),
  algorithm = c("PHYLOTAX", "Consensus", "PHYLOTAX+Cons")
) %>%
  mutate_at("preguild_taxa", syms)

# Data frame for mapping between targets for building phyloseq experiment objects
physeq_meta <-
  tidyr::crossing(
    dplyr::select(datasets, "seq_run", "tech", "dataset", "amplicon"),
    guild = c("fungi", "ecm") #, "ecm2", "ecm3")
  ) %>%
  dplyr::filter(
    seq_run != "is_057", # couldn't be demultiplexed
  ) %>%
  dplyr::select(-dataset) %>%
  left_join(select(guilds_meta, amplicon, algorithm), by = "amplicon") %>%
  mutate(
    fungi = glue::glue("fungi_{algorithm}"),
    ecm = glue::glue("ecm_{algorithm}"),
    ecm2 = glue::glue("ecm2_{algorithm}"),
    ecm3 = glue::glue("ecm3_{algorithm}")
  ) %>%
  mutate_at(vars(starts_with("ecm"), starts_with("fungi")), compose(syms, make.names))

latexfiles <- c("transect_paper.tex", "transect_supplement.tex")
latexdirs <- file.path(c("transect_paper_files", "transect_supplement_files"), "figure-latex")

# Drake plan
plan2 <- drake_plan(
  # targets which are imported from first plan
  file_meta = target(
    transform = map(primer_ID = !!unique(itsx_meta$primer_ID)),
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
  taxon_table = target(trigger = trigger(mode = "blacklist")),
  taxon = target(
    transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID),
    trigger = trigger(mode = "blacklist")
  ),
  # raxml_decipher_LSU = target(trigger = trigger(mode = "blacklist")),
  # raxml_decipher_long = target(trigger = trigger(mode = "blacklist")),
  raxml_decipher_unconst_long = target(trigger = trigger(mode = "blacklist")),
  # raxml_epa_full = target(trigger = trigger(mode = "blacklist")),
  # raxml_infernal_32S = target(trigger = trigger(mode = "blacklist")),
  big_seq_table = target(
    transform = map(region = !!(region_meta$region), .tag_in = step, .id = region),
    trigger = trigger(mode = "blacklist")
  ),
  chimeras = target(
    transform = map(.data = !!(dada_meta[,c("seq_run", "region")]), .id = c(seq_run, region)),
    trigger = trigger(mode = "blacklist")
  ),
  allchimeras = target(
    transform = map(region = !!(region_meta$region), .id = region),
    trigger = trigger(mode = "blacklist")
  ),
  qstats = target(trigger = trigger(mode = "blacklist")),
  
  # big_fasta ----
  # write the big_seq_table as a fasta files so that they can be clustered by
  # VSEARCH.
  big_fasta = target(
    write_big_fasta(big_seq_table,
                    file_out(big_fasta_file)),
    transform = map(.data = !!region_meta, .tag_in = step, .id = region)
  ),
  
  # Create labels for the tree(s) which show the assigned taxonomy.
  taxon_labels = make_taxon_labels(taxon_table),
  
  # funguild_db ----
  # Download the FUNGuild database
  funguild_db = FUNGuildR::get_funguild_db(),
  
  # platemap ----
  # Read the map between plate locations and samples
  platemap = target(
    read_platemap(file_in(!!config$platemap), !!config$platemap_sheet),
    format = "fst"),
  
  # write out tables
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
  
  # tree_decipher_LSU = target(
  #   raxml_decipher_LSU$bipartitions,
  #   transform = map(
  #     treename = "decipher_LSU",
  #     group = "euk",
  #     .tag_in = step, .tag_out = c(euktree, tree), .id = FALSE)
  # ),
  # 
  # tree_decipher_LSU_long = target(
  #   raxml_decipher_long$bipartitions,
  #   transform = map(
  #     treename = "decipher_LSU_long",
  #     group = "euk",
  #     .tag_in = step, .tag_out = c(euktree, tree), .id = FALSE)
  # ),
  
  tree_decipher_unconst_long = target(
    raxml_decipher_unconst_long$bipartitions,
    transform = map(
      treename = "decipher_unconst_long",
      group = "euk",
      .tag_in = step, .tag_out = c(euktree, tree), .id = FALSE)
  ),
  
  # tree_infernal_32S = target(
  #   raxml_infernal_32S$bipartitions,
  #   transform = map(
  #     treename = "infernal_32S",
  #     group = "euk",
  #     .tag_in = step, .tag_out = c(euktree, tree), .id = FALSE)
  # ),
  
  # tree_epa_mafft_full = target(
  #   raxml_epa_mafft_full$bestTree,
  #   transform = map(treename = "epa_mafft_full", .tag_out = tree)
  # ),
  
  labeled_tree = target(
    relabel_tree(
      tree = tree,
      old = taxon_labels$label,
      new = taxon_labels$tip_label,
      chimeras = allchimeras_ITS2
    ) %T>%
      castor::write_tree(file_out(!!glue::glue(
        "data/trees/{treename}_{group}.tree",
        treename = treename,
        group = group
      ))),
    transform = map(tree, tree, group, .tag_in = step, .id = c(treename, group))
  ),
  
  phylo_labeled_tree = target(
    relabel_tree(
      tree = tree,
      old = phylotaxon_labels$label,
      new = phylotaxon_labels$tip_label,
      chimeras = allchimeras_ITS2
    ) %T>%
      castor::write_tree(file_out(!!glue::glue(
        "data/trees/{treename}_{group}_phylo.tree",
        treename = treename,
        group = group
      ))),
    transform = map(tree, phylotaxon_labels, treename, group, .tag_in = step, .id = c(treename, group))
  ),
  
  # find the most abundant sequence which we are confident is not a fungus
  best_nonfungus =
    taxon_table %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank == "kingdom",
                  taxon != "Fungi",
                  n_tot == 6,
                  n_diff == 1) %>%
    dplyr::filter(n_reads == max(n_reads)) %$%
    label[1],
  
  # root the tree outside the fungi.  This is almost certainly not the correct
  # root for Eukaryotes (inside Chlorophyta), but it ensures that
  # the fungi can be indentified.
  rooted_tree = target(
    ape::root.phylo(euktree, best_nonfungus),
    transform = map(euktree, .tag_in = step, .id = treename)
  ),
  
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
  
  # Extract the minimally inclusive clade including all confidently identified
  # fungi (except Tulasnella).
  fungi_tree = target(
    ape::getMRCA(rooted_tree, tip = intersect(surefungi, rooted_tree$tip.label)) %>%
      ape::extract.clade(phy = rooted_tree) %>%
      ape::drop.tip(phy = ., intersect(.$tip.label, allchimeras_ITS2)),
    transform = map(rooted_tree, group = "fungi", .tag_in = step, .tag_out = tree, .id = treename)
  ),
  
  phylotaxon = target(
    phylotax(tree = tree, taxa = taxon_table),
    transform = map(tree, .tag_in = step,
                    taxname = "long",
                    algorithm = "PHYLOTAX",
                    .tag_out = consensus_taxa, .id = c(treename, group))
  ),
  
  phylotaxon_labels = target(
    make_taxon_labels(phylotaxon$tip_taxa),
    transform = map(phylotaxon, .tag_in = step, .id = c(treename, group))
  ),
  
  strict_taxon_short_euk = target(
    taxon_table %>%
      group_by(label, method, region, reference) %>%
      arrange(rank) %>%
      filter(cumall(n_diff == 1)) %>%
      ungroup() %>%
      mutate(
        name = "consensus",
        region = "All",
        reference = "All",
        ref_region = "All",
        method = "consensus",
        confidence = NA
      ) %>%
      unique() %>%
      list(tip_taxa = .),
    
    transform = map(
      taxname = "short",
      group = "euk",
      algorithm = "Consensus",
      .tag_out = c(consensus_taxa, strict_taxon),
      .id = FALSE
    )
  ),
  
  strict_taxon_short_fungi = target(
    group_by(strict_taxon_short_euk$tip_taxa, label) %>%
      filter("Fungi" %in% taxon | any(endsWith(taxon, "mycota"))) %>%
      list(tip_taxa = .),
    transform = map(
      taxname = "short",
      group = "fungi",
      algorithm = "Consensus",
      .tag_out = c(consensus_taxa, strict_taxon),
      .id = FALSE
    )
  ),
  
  best_taxon_short_euk = target(
    strict_taxon_short_euk$tip_taxa %>%
      filter(!label %in% tree_decipher_unconst_long$tip.label) %>%
      bind_rows(
        phylotaxon_decipher_unconst_long_euk$tip_taxa %>%
          filter(method == "phylotax")
      ) %>%
      mutate(method = "phylotax+c") %>%
      list(tip_taxa = .),
    transform = map(group = "euk",
                    taxname = "short",
                    algorithm = "PHYLOTAX+Cons",
                    .tag_out = consensus_taxa, .id = FALSE)
  ),
  
  best_taxon_short_fungi = target(
    strict_taxon_short_fungi$tip_taxa %>%
      filter(!label %in% fungi_tree_decipher_unconst_long$tip.label) %>%
      bind_rows(
        phylotaxon_decipher_unconst_long_fungi$tip_taxa %>%
          filter(method == "phylotax")
      ) %>%
      mutate(method = "phylotax+c") %>%
      list(tip_taxa = .),
    transform = map(group = "fungi",
                    taxname = "short",
                    algorithm = "PHYLOTAX+Cons",
                    .tag_out = consensus_taxa, .id = FALSE)
  ),
  strict_taxon_ITS_euk = target(
    taxon_table %>%
      filter(ref_region == "ITS") %>%
      group_by(label, rank) %>%
      mutate(n_diff = n_distinct(taxon)) %>%
      group_by(label, method, region, reference) %>%
      arrange(rank) %>%
      filter(cumall(n_diff == 1)) %>%
      ungroup() %>%
      mutate(
        name = "consensus",
        region = "ITS",
        reference = "All",
        ref_region = "ITS",
        method = "consensus",
        confidence = NA
      ) %>%
      unique() %>%
      list(tip_taxa = .),
    transform = map(
      taxname = "ITS",
      group = "euk",
      algorithm = "Consensus",
      .tag_out = c(consensus_taxa, strict_taxon),
      .id = FALSE
    )
  ),
  
  strict_taxon_ITS_fungi = target(
    group_by(strict_taxon_ITS_euk$tip_taxa, label) %>%
      filter("Fungi" %in% taxon | any(endsWith(taxon, "mycota"))) %>%
      list(tip_taxa = .),
    transform = map(
      taxname = "ITS",
      group = "fungi",
      algorithm = "Consensus",
      .tag_out = c(consensus_taxa, strict_taxon),
      .id = FALSE
    )
  ),
  
  best_taxon_ITS_euk = target(
    strict_taxon_ITS_euk$tip_taxa %>%
      filter(!label %in% tree_decipher_unconst_long$tip.label) %>%
      bind_rows(
        phylotaxon_decipher_unconst_long_euk$tip_taxa %>%
          filter(method == "phylotax")
      ) %>%
      mutate(method = "phylotax+c") %>%
      list(tip_taxa = .),
    transform = map(group = "euk",
                    taxname = "ITS",
                    algorithm = "PHYLOTAX+Cons",
                    .tag_out = consensus_taxa, .id = FALSE)
  ),
  
  best_taxon_ITS_fungi = target(
    strict_taxon_ITS_fungi$tip_taxa %>%
      filter(!label %in% fungi_tree_decipher_unconst_long$tip.label) %>%
      bind_rows(
        phylotaxon_decipher_unconst_long_fungi$tip_taxa %>%
          filter(method == "phylotax")
      ) %>%
      mutate(method = "phylotax+c") %>%
      list(tip_taxa = .),
    transform = map(group = "fungi",
                    taxname = "ITS",
                    algorithm = "PHYLOTAX+Cons",
                    .tag_out = consensus_taxa, .id = FALSE)
  ),
  
  preguild_taxa = target(
    consensus_taxa$tip_taxa %>%
      dplyr::group_by(label, rank) %>%
      dplyr::filter((!"ITS" %in% region) | region != "Short") %>%
      dplyr::mutate(
        n_tot = dplyr::n(),
        n_diff = dplyr::n_distinct(taxon, na.rm = TRUE),
        n_method = dplyr::n_distinct(method, na.rm = TRUE),
        n_reference = dplyr::n_distinct(reference, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(label, rank, n_diff, taxon) %>%
      dplyr::mutate_at("rank", rank_factor) %>%
      dplyr::arrange(rank) %>%
      dplyr::group_by(label) %>%
      dplyr::filter(dplyr::cumall(n_diff == 1)) %>%
      unique() %>%
      dplyr::ungroup() %>%
      tidyr::spread(., key = rank, value = taxon) %>%
      dplyr::mutate(
        kingdom = dplyr::if_else(
          !is.na(phylum) & endsWith(phylum, "mycota"),
          "Fungi",
          kingdom
        )
      ) %>%
      dplyr::left_join(
        .,
        tidyr::gather(., key = "rank", value = "taxon", kingdom:genus) %>%
          dplyr::group_by(label) %>%
          dplyr::filter(!is.na(taxon)) %>%
          dplyr::summarize(Taxonomy = paste(taxon, collapse = ";")),
        by = "label"
      ),
    transform = map(consensus_taxa, .tag_in = step, .id = c(taxname, group, algorithm))
  ),
  
  guilds = target(
      FUNGuildR::funguild_assign(preguild_taxa, funguild_db),
    transform = map(.data = !!guilds_meta, .tag_in = step, .id = algorithm)
  ),
  
  fungi = target(
    dplyr::filter(guilds, kingdom == "Fungi"),
    transform = map(guilds, .tag_in = step, .id = algorithm)
  ),
  
  ecm = target(
    dplyr::filter(guilds, grepl("Ectomycorrhizal", guild)),
    transform = map(guilds, .tag_in = step, .id = algorithm)
  ),
  
  ecm2 = target(
    dplyr::filter(ecm, !taxon %in% c("Peziza", "Pezizaceae", "Pyronemataceae")),
    transform = map(ecm, .tag_in = step, .id = algorithm)
  ),
  
  ecm3 = target(
    dplyr::filter(ecm, confidenceRanking != "Possible"),
    transform = map(ecm, .tag_in = step, .id = algorithm)
  ),
  
  proto_physeq = assemble_physeq(
    platemap = platemap,
    datasets = datasets,
    seqtable = relabel_seqtable(big_seq_table_ITS2),
    tree = NULL,
    chimeras = allchimeras_ITS2
  ),
  
  physeq = target(
    {
      physeq <- proto_physeq %>%
        phyloseq::prune_taxa(taxa = !phyloseq::taxa_names(.) %in% pos_control)
      if (amplicon == "Long") {
        phyloseq::phy_tree(physeq) <- fungi_tree_decipher_unconst_long
      }
      physeq <- physeq %>%
        phyloseq::subset_samples(sample_type == "Sample") %>%
        phyloseq::subset_samples(buffer == "Xpedition") %>%
        phyloseq::subset_samples(site == "Gan") %>%
        phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 100) %>%
        phyloseq::prune_samples(samples = phyloseq::sample_data(.)[["tech"]] == tech) %>%
        phyloseq::prune_samples(samples = phyloseq::sample_data(.)[["amplicon"]] == amplicon)
      
      if (guild == "ecm") {
        physeq <- phyloseq::prune_taxa(ecm$label, physeq) %>%
          phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0)
        
      } else if (guild == "ecm2") {
        physeq <- phyloseq::prune_taxa(ecm2$label, physeq) %>%
          phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0)
        
      } else if (guild == "ecm3") {
        physeq <- phyloseq::prune_taxa(ecm3$label, physeq) %>%
          phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0)
        
      } else if (guild == "fungi") {
        physeq <- phyloseq::prune_taxa(fungi$label, physeq) %>%
          phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0)
      }
      phyloseq::transform_sample_counts(physeq, function(x) x / sum(x))
    },
    transform = map(.data = !!physeq_meta, .tag_in = step, .id = c(guild, tech, amplicon, algorithm))
  ),
  
  correlog = target(
    correlog(physeq, metric, timelag),
    transform = cross(physeq, metric = c("bray", "wunifrac"), timelag = c(0, 1),
                      .tag_in = step, .id = c(guild, metric, tech, amplicon, algorithm, timelag))
  ),
  
  variog = target(
    variog(physeq, metric, breaks = c(1:25 - 0.5, 31000)),
    transform = cross(physeq, metric = c("bray", "wunifrac"),
                      .tag_in = step,
                      .id = c(guild, metric, tech, amplicon, algorithm))
  ),
  
  variofit2 = target(
    variog %>%
      as_variogram() %>%
      stats::nls(
        gamma ~ (1 - sill) - (1 - sill) * (1 - nugget) * exp(dist/range*log(0.05)),
        data = .,
        start = list(
          nugget = min(.$gamma),
          sill = 1 - max(.$gamma),
          range = 30
        ),
        upper = list(nugget = 0.99, sill = 0.99, range = Inf),
        lower = list(nugget = 0, sill = 0, range = 0.001),
        algorithm = "port",
        weights = pmin(.$np/.$dist, 1),
        control = list(warnOnly = TRUE, maxiter = 1000, minFactor = 1/4096)
      ) %>% {
        list(
          pars = as.list(.$m$getPars()),
          summary = variog_summary(.),
          predict = variog_predict(
            .,
            newdata = expand_grid(dist = seq(1, 25, 0.5), timelag = 0:1)
          )
        )
      },
    transform = map(variog, .tag_in = step,.id = c(guild, metric, tech, amplicon, algorithm))
  ),
  
  variogST = target(
    variogST(physeq, metric, breaks = c(1:25 - 0.5, 30000)),
    transform = cross(physeq, metric = c("bray", "wunifrac"),
                      .tag_in = step,
                      .id = c(guild, metric, tech, amplicon, algorithm))
  ),
  
  variofitST2 = target(
    as_variogramST(variogST) %>%
      stats::nls(
        gamma ~ (1 - sill) - (1 - sill) * (1 - nugget) * exp((dist/range + timelag/timerange)*log(0.05)),
        # add the parameters from the spatial-only fit.
        data = .,
        start = c(
          list(timerange = 2),
          variofit2$pars
        ),
        lower = c(
          list(timerange = 0.001),
          variofit2$pars
        ),
        algorithm = "port",
        weights = pmin(.$np/.$dist, 1),
        control = list(warnOnly = TRUE, maxiter = 1000, minFactor = 1/4096)
      ) %>% {
        list(
          pars = as.list(.$m$getPars()),
          summary = variog_summary(.),
          predict = variog_predict(
            .,
            newdata = expand_grid(dist = seq(1, 25, 0.5), timelag = 0:1)
          )
        )
      },
    transform = map(variogST, variofit2, .tag_in = step, .id = c(guild, metric, tech, amplicon, algorithm))
  ),
  
  ##### Tables and Figures #####################################################
  # Targets below this point are intended to compile data in useful formats    #
  # for making the tables, figures, and in-line results in the paper.          #
  ##############################################################################
  
  # list of all sequences by region and sequencing run, along with number of
  # reads and sequence length
  reads_table = target(
    list(big_seq_table) %>%
      purrr::map(tibble::as_tibble, rownames = "sample") %>%
      purrr::map_dfr(tidyr::pivot_longer, -"sample", names_to = "seq", values_to = "reads") %>%
      dplyr::filter(reads > 0) %>%
      tidyr::extract(
        "sample",
        c("seq_run", "plate", "well", "region"),
        regex = "([piS][bsH][-_]\\d{3,4})_(\\d{3})_([A-H]1?\\d)_(.+)"
      ) %>%
      dplyr::group_by(seq_run, region, seq) %>%
      dplyr::summarize(reads = sum(reads)) %>%
      dplyr::mutate(length = nchar(seq)) %>%
      dplyr::ungroup(),
    transform = combine(big_seq_table)
  ),
  
  # Summary information for each region
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
  
  # Table of reads per ASV for each sequencing run
  asv_table = big_seq_table_ITS2 %>%
    tibble::as_tibble(rownames = "filename") %>%
    tidyr::gather(key = "seq", value = "reads", -1) %>%
    dplyr::filter(reads >= 1) %>%
    tidyr::extract(col = "filename",
                   into = c("seq.run", "plate", "well", "dir", "region"),
                   regex = "([a-zA-Z]+[-_]\\d+)_(\\d+)_([A-H]1?[0-9])([fr]?)_([:alnum:]+).+") %>%
    dplyr::group_by(seq.run, seq) %>%
    dplyr::summarize(reads = sum(reads)) %>%
    tidyr::spread(key = seq.run, value = reads, fill = 0) %>% 
    dplyr::ungroup(),
  
  # Table of reads per OTU for each sequencing run
  otu_table = {
    
    big_fasta_ITS2
    read_tsv(
      file_in("data/clusters/ITS2.table"),
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
      dplyr::filter(reads >= 1) %>%
      dplyr::group_by(seq_run, seq) %>%
      dplyr::summarize(reads = sum(reads)) %>%
      tidyr::spread(key = seq_run, value = reads, fill = 0)
  },
  
  # Find the central sequence for each 97% OTU
  otu_map =
    read_tsv(
      file_in(here::here("data/clusters/ITS2.uc")),
      col_names = paste0("V", 1:10),
      col_types = "ciidfccccc",
      na = c("", "NA", "*")
      ) %>%
    filter(V1 == "H") %>%
    select(identity = V4, ASVseq = V9, OTUseq = V10) %>%
    mutate_all(str_replace, ";.*", "") %>%
    unique(),
  
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
  
  readcounts = parsed_qstat %>% filter(is.na(read) | read != "_R2", stat == "n"),
  
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
  
  asvcounts =
    pivot_longer(asv_table, -1, names_to = "seq_run", values_to = "reads") %>%
    filter(reads > 0) %>%
    mutate_at("seq", chartr, old = "T", new = "U") %>%
    inner_join(
      select(allseqs, seq = ITS2, everything()) %>%
        select(-hash, -nread) %>%
        mutate(ITS2 = seq) %>%
        pivot_longer(-1, names_to = "region", values_to = "consensus") %>%
        filter(!is.na(consensus)) %>%
        select(-consensus) %>%
        unique(),
      by = "seq"
    ) %>%
    group_by(seq_run, region) %>%
    summarize(reads = sum(reads, na.rm = TRUE), ASVs = n()) %>%
    select(seq_run, reads, step = region, ASVs),
  
  bioinf_table = bind_rows(
    enframe(rawcounts, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "Raw"
      ),
    enframe(demuxcounts, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "Trim"
      ),
    enframe(filtercounts_full, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "Filter (full)"
      ),
    enframe(regioncounts, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "LSUx"
      ),
    enframe(filtercounts_ITS2, name = "seq_run", value = "reads") %>%
      mutate(
        # reads = as.integer(gsub(" ", "", reads)),
        step = "Filter (ITS2)"
      ),
    filter(region_table, region == "ITS2") %>%
      select(seq_run, reads, step = region, ASVs),
    pivot_longer(asv_table, -1, names_to = "seq_run", values_to = "reads") %>%
      filter(reads > 0) %>%
      mutate_at("seq", chartr, old = "T", new = "U") %>%
      inner_join(
        select(allseqs, seq = ITS2, long, short, ITS, LSU) %>%
          pivot_longer(-1, names_to = "region", values_to = "consensus") %>%
          filter(!is.na(consensus)) %>%
          select(-consensus) %>%
          unique(),
        by = "seq"
      ) %>%
      group_by(seq_run, region) %>%
      summarize(reads = sum(reads, na.rm = TRUE), ASVs = n()) %>%
      select(seq_run, reads, step = region, ASVs),
    # pivot_longer(otu_table, -1, names_to = "seq_run", values_to = "reads") %>%
    #   filter(reads > 0) %>%
    #   mutate_at("seq", chartr, old = "T", new = "U") %>%
    #   inner_join(
    #     select(readd(allseqs, cache = cache), seq = ITS2, long, short, ITS, LSU) %>%
    #       pivot_longer(-1, names_to = "region", values_to = "consensus") %>%
    #       filter(!is.na(consensus)) %>%
    #       select(-consensus) %>%
    #       unique(),
    #     by = "seq"
    #   ) %>%
    #   group_by(seq_run, region) %>%
    #   summarize(reads = sum(reads, na.rm = TRUE), OTUs = n()) %>%
    #   select(seq_run, reads, step = region, OTUs)
  ) %>%
    left_join(select(datasets, seq_run, tech, amplicon), by = "seq_run") %>%
    select(tech, amplicon, step, reads, ASVs) %>%
    pivot_longer(c("reads", "ASVs"), names_to = "type", values_to = "count") %>%
    filter(!is.na(count)) %>%
    mutate(
      step = factor(step, levels = c("Raw", "Trim", "Filter (full)", "LSUx", "Filter (ITS2)", "ITS2",
                                     "short", "ITS", "LSU", "long")),
      tech = factor(tech, levels = c("PacBio", "Ion Torrent", "Illumina")),
      amplicon = factor(stringr::str_to_title(amplicon), levels = c("Long", "Short")),
      type = factor(type, c("ASVs", "reads"))
    ) %>%
    arrange(tech, amplicon, type) %>%
    filter(ifelse(amplicon == "Short", !step  %in% c("ITS", "LSU", "long"), step != "short")) %>%
    pivot_wider(
      id_cols = "step",
      names_from = c("tech", "amplicon", "type"),
      values_from = c("count")
    ) %>%
    arrange(step) %>%
    column_to_rownames("step"),
  
  bioinf_header =
    tibble(name = names(bioinf_table)) %>%
    separate(name, into = c("tech", "amplicon", "type"), sep = "_") %>%
    mutate(
      tech = paste(tech, plyr::mapvalues(tech, datasets$tech, datasets$machine, FALSE))
    ),
  
  venn_ASV = venndata(asv_table, ASVs, cols = c("pb_500", "pb_483", "SH-2257", "is_057")),
  
  vennplot_ASV = vennplot_data(venn_ASV, ASVs),
  
  venn_OTU = venndata(otu_table, OTUs, cols = c("pb_500", "pb_483", "SH-2257", "is_057")),
  
  vennplot_OTU = vennplot_data(venn_OTU, OTUs),
  
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
  
  taxon_reads = {
    out <- 
      # Get the number of reads per sequencing run for each ASV
      proto_physeq %>%
      phyloseq::merge_samples(group = "seq_run") %>%
      phyloseq::otu_table() %>%
      t() %>%
      as("matrix") %>%
      as_tibble(rownames = "label") %>%
      # Normalize reads to a fraction of the sequencing run.
      mutate_if(is.numeric, ~./sum(.)) %>%
      # Pivot so that each ASV/sequencing run combination is its own line, and remove 0's
      pivot_longer(-1, names_to = "seq_run", values_to = "reads") %>%
      filter(reads > 0) %>%
      # Add the amplicon and technology used for each sequencing run
      left_join(select(datasets, seq_run, amplicon, tech), by = "seq_run") %>%
      # put the amplicon in title case for display
      mutate_at("amplicon", stringr::str_to_title) %>%
      # Make all valid combinations of reference, region, and method
      expand_grid(
        tibble(
          reference = c("unite", "warcup", "rdp_train"),
          region = c("ITS", "ITS", "LSU")
        ),
        method = c("dada2", "idtaxa", "sintax")
      ) %>%
      filter(
        ifelse(amplicon == "Short", region == "ITS", TRUE),
        ifelse(amplicon == "Long", region %in% c("ITS", "LSU"), TRUE)
      ) %>%
      # Add the identifications from the taxon table.
      # This will leave NAs where a method did not make an identification at a rank,
      # Including full NA rows where no identification was made at all.
      left_join(
        taxon_table %>%
          mutate(region = ifelse(region == "short", "ITS", region)) %>%
          mutate_at("taxon", na_if, "NA") %>%
          select(label, method, reference, rank, taxon, region) %>%
          pivot_wider(
            names_from = "rank",
            values_from = "taxon"
          ) %>%
          mutate(
            family = ifelse(is.na(genus), family, coalesce(family, "incertae_sedis")),
            order = ifelse(is.na(family), order, coalesce(order, "incertae_sedis")),
            class = ifelse(is.na(order), class, coalesce(class, "incertae_sedis")),
            phylum = ifelse(class == "Leotiomycetes", "Ascomycota", phylum)
          )
        ,
        by = c("label", "reference", "region", "method")
      ) %>%
      # Warcup doesn't explicitly identify fungi.
      mutate(
        kingdom = ifelse(!is.na(phylum) & reference == "warcup", "Fungi", kingdom)
      )
    
    # Add strict consensus and PHYLOTAX identifications
    bind_rows(
      out,
      phylotaxon_decipher_unconst_long_euk$tip_taxa %>%
        filter(method == "phylotax") %>%
        select(method, label, rank, taxon) %>%
        pivot_wider(names_from = "rank", values_from = "taxon") %>%
        left_join(
          out %>%
            filter(seq_run == "pb_500") %>%
            select(seq_run, label, reads, amplicon, tech) %>%
            unique(),
          .,
          by = "label"
        ) %>%
        mutate_at("kingdom", na_if, "NA") %>%
        mutate(
          method = "PHYLO",
          reference = "All",
          kingdom = ifelse(endsWith(phylum, "mycota"), "Fungi", kingdom)
        ),
      strict_taxon_short_euk$tip_taxa %>%
        select(method, label, rank, taxon, region) %>%
        pivot_wider(names_from = "rank", values_from = "taxon") %>%
        inner_join(
          out %>%
            filter(amplicon == "Short") %>%
            select(seq_run, label, reads, amplicon, tech) %>%
            unique(),
          .,
          by = "label"
        ) %>%
        mutate_at("kingdom", na_if, "NA") %>%
        mutate(
          kingdom = ifelse(endsWith(phylum, "mycota"), "Fungi", kingdom)
        ),
      strict_taxon_ITS_euk$tip_taxa %>%
        select(method, label, rank, taxon, region) %>%
        pivot_wider(names_from = "rank", values_from = "taxon") %>%
        inner_join(
          out %>%
            select(seq_run, label, reads, amplicon, tech) %>%
            unique(),
          .,
          by = "label"
        ) %>%
        mutate_at("kingdom", na_if, "NA") %>%
        mutate(
          kingdom = ifelse(endsWith(phylum, "mycota"), "Fungi", kingdom)
        ),
      best_taxon_short_euk$tip_taxa %>%
        select(method, label, rank, taxon) %>%
        pivot_wider(names_from = "rank", values_from = "taxon") %>%
        inner_join(
          out %>%
            filter(amplicon == "Short") %>%
            select(seq_run, label, reads, amplicon, tech) %>%
            unique(),
          .,
          by = "label"
        ) %>%
        mutate_at("kingdom", na_if, "NA") %>%
        mutate(
          kingdom = ifelse(endsWith(phylum, "mycota"), "Fungi", kingdom)
        )
    ) %>%
      mutate(
        family = ifelse(is.na(genus), family, coalesce(family, "incertae_sedis")),
        order = ifelse(is.na(family), order, coalesce(order, "incertae_sedis")),
        class = ifelse(is.na(order), class, coalesce(class, "incertae_sedis")),
        phylum = ifelse(class == "Leotiomycetes", "Ascomycota", phylum) # Not sure what happened here...
      ) %>%
      mutate_at(
        "method",
        factor,
        levels = c("PHYLO", "phylotax+c", "consensus", "dada2", "sintax", "idtaxa"),
        labels = c("PHYLOTAX", "PHYLOTAX", "Consensus", "RDPC", "SINTAX", "IDTAXA")
      ) %>%
      mutate_at(
        "reference",
        replace_na,
        "All"
      ) %>%
      mutate_at(
        "reference",
        factor,
        levels = c("All", "unite", "warcup", "rdp_train"),
        labels = c("All", "Unite", "Warcup", "RDP")
      ) %>%
      mutate_at(
        "reference",
        replace_na,
        "All"
      ) %>%
      mutate_at("tech", factor, levels = c("PacBio", "Illumina", "Ion Torrent")) %>%
      rename(Algorithm = method) %>%
      mutate(
        Taxonomy = paste(kingdom, phylum, class, order, family, genus, sep = ";") %>%
          str_replace_all(";NA", "")
      ) %>%
      bind_cols(
        FUNGuildR::funguild_assign(., funguild_db) %>%
          mutate(
            ECM = grepl("Ectomycorrhizal", guild) %>%
              ifelse(paste(confidenceRanking, "ECM"),
                     ifelse(is.na(guild), NA, "non-ECM")) %>%
              ifelse(is.na(family), NA, .)) %>%
          mutate_at(
            "ECM",
            factor,
            levels = c(NA_character_, "non-ECM", "Possible ECM",
                       "Probable ECM", "Highly Probable ECM"),
            exclude = NULL,
            ordered = TRUE
          ) %>%
          select(ECM)
      ) %>%
      select(-Taxonomy) %>%
      left_join(select(otu_map, label = ASVseq, OTU = OTUseq), by = "label")
  },
  
  tax_chart =
    taxon_reads %>%
    pivot_longer(
      names_to = "rank",
      cols = c("kingdom", "phylum", "class", "order", "family", "genus"),
      values_to = "taxon"
    ) %>%
    filter(ifelse(Algorithm == "Consensus", region == "All", TRUE)) %>%
    group_by(amplicon, tech, Algorithm, reference, rank, ID = !is.na(taxon)) %>%
    summarize(reads = sum(reads), ASVs = n()) %>%
    group_by(amplicon, tech, Algorithm, reference, rank) %>%
    mutate(ASVs = ASVs / sum(ASVs)) %>%
    ungroup() %>%
    filter(ID) %>%
    mutate_at("rank", factor,
              levels = c("kingdom", "phylum", "class", "order", "family", "genus")
    ) %>%
    rename(Reads = reads) %>%
    mutate_at("reference", as.character) %>%
    mutate_at("reference", replace_na, "All") %>%
    mutate(
      reference = map2(
        reference,
        amplicon,
        ~ if (.x == "All") {
          c("Unite", "Warcup", if (.y == "Long") "RDP" else NULL)
          } else {.}
      )
    ) %>%
    unnest(reference) %>%
    mutate_at("reference", factor, levels = c("Unite", "Warcup", "RDP")) %>%
    pivot_longer(cols = c("Reads", "ASVs"), names_to = "type", values_to = "frac"),
  
  tax_chart_plot =
    tax_chart %>%
    ggplot(aes(x = rank, y = frac, ymax = frac, group = Algorithm, fill = Algorithm)) +
    ggnomics::facet_nested(
      tech + amplicon ~ type + reference,
      resect = unit(1, "mm"),
      nest_line = TRUE,
      bleed = FALSE
    ) +
    # geom_ribbon(ymin = 0, alpha = 0.5) +
    # geom_line(aes(color = Algorithm)) +
    geom_col(position = "dodge") +
    # scale_fill_brewer(type = "qual", palette = 2) +
    ylab("Fraction of reads assigned") +
    xlab(NULL) +
    theme(axis.text.x = element_text(vjust = 0.5, angle = 90),
          legend.direction = "horizontal",
          legend.position = "bottom",
          strip.background = element_blank()),
  
  agaricus =
    taxon_table %>%
    filter(rank == "genus", taxon == "Agaricus", n_diff == 1) %$%
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
        metric = plyr::mapvalues(metric, c("bray", "wunifrac"), c("Bray-Curtis", "W-UNIFRAC")),
        amplicon = factor(
          amplicon,
          levels = c("Short", "Long")
        ),
        algorithm = factor(
          algorithm,
          levels = c("Consensus", "PHYLOTAX", "PHYLOTAX+Cons"),
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
  
  taxdata = {
    ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
    out <- taxon_reads %>%
      filter(region == "ITS", Algorithm == "Consensus") %>%
      mutate(domain = "Root", ASVs = 1) %>%
      select(-region, -seq_run, -label, -OTU) %>%
      mutate_at(
        ranks,
        replace_na,
        "?"
      ) %>%
      group_by_at(vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
      mutate(ASVs = ASVs / sum(ASVs))
    for (i in seq_along(ranks)) {
      out <- out %>%
        group_by_at(rev(ranks)[i:length(ranks)]) %>%
        group_map(
          function(x, y, i) {
            d <- group_by_at(x, vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
              summarize_at(c("reads", "ASVs"), sum)
            d2 <- group_by_at(d, "amplicon") %>% summarize_at(c("reads", "ASVs"), mean)
            reads <- select(d2, "amplicon", "reads") %>% deframe()
            ASVs <- select(d2, "amplicon", "ASVs") %>% deframe()
            if (max(c((d$reads), (d$ASVs))) < 0.01) {
              if (
                (isTRUE(abs(log10(reads["Long"]) - log10(reads["Short"])) < 1) &
                 isTRUE(abs(log10(ASVs["Long"]) - log10(ASVs["Short"])) < 1)) |
                max(c((d$reads), (d$ASVs))) < 0.01
              ) {
                mutate_at(x, rev(ranks)[i], ~"*")
              } else if (i > 1) {
                mutate_at(x, rev(ranks)[i - 1], ~"*") %>%
                mutate_at(rev(ranks)[i], ~paste0("(",., ")"))
              } else {
                mutate_at(x, rev(ranks)[i], ~paste0("(",., ")"))
              }
            } else {
              x
            }
          },
          i = i,
          keep = TRUE
        ) %>%
        bind_rows()
    }
    out <- group_by_at(out, vars(-one_of("reads", "ASVs", "ECM"))) %>%
      summarize(reads = sum(reads), ASVs = sum(ASVs), ECM = min(ECM)) %>%
      pivot_wider(
        names_from = c("tech", "amplicon", "reference", "Algorithm"),
        values_from = c("reads", "ASVs"),
        values_fill = list(reads = 0, ASVs = 0),
        names_repair = "universal"
      ) %>% 
      taxa::parse_tax_data(
        class_cols = ranks
      ) %>%
      taxa::filter_taxa(
        taxon_names != "incertae_sedis",
        drop_obs = TRUE,
        reassign_obs = FALSE,
        reassign_taxa = TRUE
      )
    
    read_cols <- keep(names(out$data$tax_data), startsWith, "reads_")
    asv_cols <- keep(names(out$data$tax_data), startsWith, "ASVs_")
    
    out$data$asv_table <- metacoder::calc_obs_props(out, "tax_data", cols = asv_cols)
    out$data$read_table <- metacoder::calc_obs_props(out, "tax_data", cols = read_cols)
    out$data$tax_asv <- metacoder::calc_taxon_abund(out, "asv_table", cols = asv_cols)
    out$data$tax_read <- metacoder::calc_taxon_abund(out, "read_table", cols = read_cols)
    
    child_of_unknown <-  unlist(
      out$supertaxa(
        recursive = FALSE,
        value = "taxon_names",
        na = TRUE
      )
    ) %in% c("?", "*")
    unknown_taxa <- out$taxon_names() %in% c("?", "*", "incertae_sedis")
    out <- out$filter_taxa(
      !(child_of_unknown & unknown_taxa) | taxon_names == "Root",
      reassign_obs = FALSE
    )
    dangling_taxa <- map_lgl(
      out$subtaxa(recursive = FALSE)[unlist(out$supertaxa(recursive = FALSE, na = TRUE))],
      ~all(out$taxon_names()[.] %in% c("?", "*") | is.na(out$taxon_names()[.]))
    )
    out <- out$filter_taxa(
      !dangling_taxa | taxon_names == "Root",
      reassign_obs = FALSE
    )
    
    retaxa_table <- tibble::tibble(
      classification = out$classifications(),
      id = out$taxon_ids()
    ) %>%
      group_by(classification) %>%
      mutate(newid = first(id)) %>%
      filter(id != newid)
    for (n in names(out$data)) {
      if ("taxon_id" %in% names(out$data[[n]]))
        out$data[[n]]$taxon_id <- plyr::mapvalues(
          out$data[[n]]$taxon_id,
          retaxa_table$id,
          retaxa_table$newid,
          FALSE
        )
    }
    for (n in c("asv_table", "read_table", "tax_asv", "tax_read")) {
      out$data[[n]] <- out$data[[n]] %>%
        group_by(taxon_id) %>%
        summarize_all(sum)
    }
    out <- out$filter_taxa(!taxon_ids %in% retaxa_table$id, drop_obs = TRUE,
                           reassign_obs = FALSE, reassign_taxa = FALSE)
    
    out$data$tax_read$mean = rowMeans(out$data$tax_read[,-1])
    out$data$tax_read$max = do.call(pmax, out$data$tax_read[,-1])
    out$data$tax_asv$mean = rowMeans(out$data$tax_asv[,-1])
    out$data$tax_asv$mean = do.call(pmax, out$data$tax_asv[,-1])
    
    out$data$diff_read <- metacoder::compare_groups(
      out,
      data = "tax_read",
      cols = c(read_cols, "mean"),
      groups = c(read_cols, "mean"),
      combinations = lapply(read_cols, c, "mean")
    )
    out$data$diff_asv <- metacoder::compare_groups(
      out,
      data = "tax_read",
      cols = c(read_cols, "mean"),
      groups = c(read_cols, "mean"),
      combinations = lapply(read_cols, c, "mean")
    )
    out
  },
  taxdata_ECM = {
    ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
    out <- taxon_reads %>%
      filter(region == "ITS", Algorithm == "Consensus") %>%
      filter(!is.na(as.character(ECM)), ECM != "non-ECM") %>%
      mutate(ASVs = 1) %>%
      select(-region, -seq_run, -label, -OTU) %>%
      mutate_at(
        ranks,
        replace_na,
        "?"
      ) %>%
      group_by_at(vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
      mutate(ASVs = ASVs / sum(ASVs))
    for (i in seq_along(ranks)) {
      out <- out %>%
        group_by_at(rev(ranks)[i:length(ranks)]) %>%
        group_map(
          function(x, y, i) {
            d <- group_by_at(x, vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
              summarize_at(c("reads", "ASVs"), sum)
            d2 <- group_by_at(d, "amplicon") %>% summarize_at(c("reads", "ASVs"), mean)
            reads <- select(d2, "amplicon", "reads") %>% deframe()
            ASVs <- select(d2, "amplicon", "ASVs") %>% deframe()
            if (max(c((d$reads), (d$ASVs))) < 0.01) {
              if (
                isTRUE(abs(log10(reads["Long"]) - log10(reads["Short"])) < 1) &
                isTRUE(abs(log10(ASVs["Long"]) - log10(ASVs["Short"])) < 1)
              ) {
                mutate_at(x, rev(ranks)[i], ~"*")
              } else if (i > 1) {
                mutate_at(x, rev(ranks)[i - 1], ~"*") %>%
                mutate_at(rev(ranks)[i], ~paste0("(",., ")"))
              } else {
                mutate_at(x, rev(ranks)[i], ~paste0("(",., ")"))
              }
            } else {
              x
            }
          },
          i = i,
          keep = TRUE
        ) %>%
        bind_rows()
    }
    out <- group_by_at(out, vars(-one_of("reads", "ASVs", "ECM"))) %>%
      summarize(reads = sum(reads), ASVs = sum(ASVs), ECM = min(ECM)) %>%
      pivot_wider(
        names_from = c("tech", "amplicon", "reference", "Algorithm"),
        values_from = c("reads", "ASVs"),
        values_fill = list(reads = 0, ASVs = 0),
        names_repair = "universal"
      ) %>% 
      taxa::parse_tax_data(
        class_cols = ranks
      ) %>%
      taxa::filter_taxa(
        taxon_names != "incertae_sedis",
        drop_obs = TRUE,
        reassign_obs = FALSE,
        reassign_taxa = TRUE
      )
    
    read_cols <- keep(names(out$data$tax_data), startsWith, "reads_")
    asv_cols <- keep(names(out$data$tax_data), startsWith, "ASVs_")
    
    out$data$asv_table <- metacoder::calc_obs_props(out, "tax_data", cols = asv_cols)
    out$data$read_table <- metacoder::calc_obs_props(out, "tax_data", cols = read_cols)
    out$data$tax_asv <- metacoder::calc_taxon_abund(out, "asv_table", cols = asv_cols)
    out$data$tax_read <- metacoder::calc_taxon_abund(out, "read_table", cols = read_cols)
    
    child_of_unknown <-  unlist(
      out$supertaxa(
        recursive = FALSE,
        value = "taxon_names",
        na = TRUE
      )
    ) %in% c("?", "*")
    out <- out$filter_taxa(
      !child_of_unknown | taxon_names == "Fungi",
      reassign_obs = FALSE
    )
    dangling_taxa <- map_lgl(
      out$subtaxa(recursive = FALSE)[unlist(out$supertaxa(recursive = FALSE, na = TRUE))],
      ~all(out$taxon_names()[.] %in% c("?", "*") | is.na(out$taxon_names()[.]))
    )
    out <- out$filter_taxa(
      !dangling_taxa | taxon_names == "Fungi",
      reassign_obs = FALSE
    )
    
    retaxa_table <- tibble::tibble(
      classification = out$classifications(),
      id = out$taxon_ids()
    ) %>%
      group_by(classification) %>%
      mutate(newid = first(id)) %>%
      filter(id != newid)
    for (n in names(out$data)) {
      if ("taxon_id" %in% names(out$data[[n]]))
        out$data[[n]]$taxon_id <- plyr::mapvalues(
          out$data[[n]]$taxon_id,
          retaxa_table$id,
          retaxa_table$newid,
          FALSE
        )
    }
    for (n in c("asv_table", "read_table", "tax_asv", "tax_read")) {
      out$data[[n]] <- out$data[[n]] %>%
        group_by(taxon_id) %>%
        summarize_all(sum)
    }
    out <- out$filter_taxa(!taxon_ids %in% retaxa_table$id, drop_obs = TRUE,
                           reassign_obs = FALSE, reassign_taxa = FALSE)
    
    out$data$tax_read$mean = rowMeans(out$data$tax_read[,-1])
    out$data$tax_read$max = do.call(pmax, out$data$tax_read[,-1])
    out$data$tax_asv$mean = rowMeans(out$data$tax_asv[,-1])
    out$data$tax_asv$mean = do.call(pmax, out$data$tax_asv[,-1])
    
    out$data$diff_read <- metacoder::compare_groups(
      out,
      data = "tax_read",
      cols = c(read_cols, "mean"),
      groups = c(read_cols, "mean"),
      combinations = lapply(read_cols, c, "mean")
    )
    out$data$diff_asv <- metacoder::compare_groups(
      out,
      data = "tax_read",
      cols = c(read_cols, "mean"),
      groups = c(read_cols, "mean"),
      combinations = lapply(read_cols, c, "mean")
    )
    out
  },
  
  # data for supplementary figures ####
  phylotax_fig = 
    treeio::as.treedata(
      ape::read.tree(text = "(A:1,((B:1,C:1):1,((E:1,F:1):1,D:1):1):1);")
    ) %>%
    tidytree::full_join(
      tibble(label = LETTERS[1:6],
             tax1 = c("unk", "Tax1", "Tax2", "Tax2", "unk", "unk"),
             tax2 = c("unk", "Tax2", "Tax2", "Tax1", "unk", "Tax1")),
      by = "label"
    ) %>%
    ggtree::ggtree() +
    ggtree::geom_tiplab(x = 4, align = 2) +
    xlim(0, 6) +
    ggtree::geom_tiplab(
      aes(color = tax1, label = tax1),
      x = 4.7,
      hjust = 0.5
    ) +
    ggtree::geom_tiplab(
      aes(color = tax2, label = tax2),
      x = 5.7,
      hjust = 0.5
    ) +
    ggtree::geom_tiplab(x = 5.2, label = "/", hjust = 0.5) +
    ggtree::geom_nodelab(aes(label = node - 6), geom = "label") +
    ggtree::geom_hilight(node = 10, fill = "steelblue", extendto = 4) +
    ggtree::geom_hilight(node = 9, fill = "tomato", extendto = 4) +
    scale_color_manual(
      values = c(Tax1 = "steelblue", Tax2 = "tomato", unk = "grey50")
    ) +
    scale_y_reverse(),
  
  length_table = 
    pivot_longer(asv_table, -1, names_to = "seq_run", values_to = "reads") %>%
    filter(reads > 0) %>%
    mutate_at("seq", chartr, old = "T", new = "U") %>%
    inner_join(
       select(allseqs, seq = ITS2, short, long) %>%
           mutate(ITS2 = seq) %>%
           pivot_longer(-1, names_to = "region", values_to = "consensus") %>%
           filter(!is.na(consensus)) %>%
           unique(),
       by = "seq"
      ) %>%
    mutate(length = nchar(consensus)) %>%
    group_by(seq_run, region, length) %>%
    summarize_at("reads", sum) %>%
    group_by(seq_run, region) %>%
    arrange(length) %>%
    mutate(reads = reads/sum(reads), ecdf = cumsum(reads)) %>%
    dplyr::left_join(dplyr::select(datasets, seq_run, tech, amplicon)) %>%
    dplyr::mutate_at("amplicon", factor, levels = c("Short", "Long")) %>%
    dplyr::mutate_at("tech", factor, levels = c("PacBio", "Illumina", "Ion Torrent")) %>%
    arrange(tech, amplicon) %>%
    mutate(group = paste(tech, amplicon) %>% factor(levels = unique(.))),
  
  full_dens_short = length_table %>%
    filter(amplicon == "Short", region == "short") %>%
    ggplot(aes(length, weight = reads, color = group, group = group)) +
    stat_density(bw = 0.5, geom = "line", position = "identity") +
    scale_color_strategy(guide = "none") +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
    ggnomics::facet_nested(rows = vars(tech), scales = "free_y") +
    theme(strip.background = element_blank()),
  
  full_dens_long = length_table %>%
    filter(amplicon == "Long", region == "long") %>%
    ggplot(aes(length, weight = reads, color = group, group = group)) +
    stat_density(bw = 0.5, geom = "line", position = "identity") +
    scale_color_strategy(guide = "none") +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
    ggnomics::facet_nested(rows = vars(tech), scales = "free_y") +
    theme(strip.background = element_blank()),
  
  full_ecdf_short = length_table %>%
    filter(amplicon == "Short", region == "short") %>%
    ggplot(aes(x = length, y = ecdf, color = group, group = group)) +
    geom_step() +
    scale_color_strategy(name = NULL) +
    scale_y_continuous(name = "ECDF") +
    scale_x_continuous(name = NULL) +
    facet_grid( amplicon ~ ., scales = "free_x") +
    theme(
      legend.position = c(1, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      strip.background = element_blank()
    ),
  
  full_ecdf_long = length_table %>%
    filter(amplicon == "Long", region == "long") %>%
    ggplot(aes(x = length, y = ecdf, color = group, group = group)) +
    geom_step() +
    scale_color_strategy(name = NULL) +
    scale_y_continuous(name = "ECDF") +
    scale_x_continuous(name = "Length (bp)") +
    facet_grid( amplicon ~ ., scales = "free_x") +
    theme(
      legend.position = c(1, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      strip.background = element_blank()
    ),
  
  its2_dens = length_table %>%
    filter(region == "ITS2") %>%
    ggplot(aes(length, weight = reads, color = group, group = group)) +
    stat_density(bw = 0.5, geom = "line", position = "identity") +
    scale_color_strategy(guide = "none") +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
    ggnomics::facet_nested(rows = vars(tech, amplicon), scales = "free_y",
                           nest_line = TRUE) +
    theme(strip.background = element_blank()),
  
  its2_ecdf = length_table %>%
    filter(region == "ITS2") %>%
    mutate(tech = "", amplicon = "") %>%
    ggplot(aes(x = length, y = ecdf, color = group, group = group)) +
    geom_step() +
    scale_color_strategy(name = NULL) +
    scale_y_continuous(name = "ECDF") +
    scale_x_continuous(name = "ITS2 Length (bp)") +
    facet_grid(rows = vars(tech, amplicon), scales = "free_y") +
    theme(
      strip.background = element_blank(),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_blank()
    ),
    
  
  region_length_fig = {
    baseregions <-
      c("ITS1", "5_8S", "ITS2", "LSU1", "D1", "LSU2", "D2", "LSU3", "D3", "LSU4")
    
    region_limits <-
      dplyr::filter(regions, region %in% baseregions) %>%
      dplyr::mutate_at("region", factor, levels = baseregions) %>%
      dplyr::mutate_at("region", fct_recode, "5.8S" = "5_8S") %>%
      tidyr::pivot_longer(
        cols = c("min_length", "max_length"),
        names_to = "limit",
        names_pattern = "(min|max)_length",
        values_to = "length"
      )
    
    reads_table %>%
      dplyr::filter(region %in% baseregions, seq_run == "pb_500") %>%
      dplyr::mutate_at("region", factor, levels = baseregions) %>%
      dplyr::mutate_at("region", fct_recode, "5.8S" = "5_8S") %>%
      ggplot(aes(length, weight = reads)) +
      stat_density(bw = 0.5, color = "#8da0cb", geom = "line") +
      geom_vline(data = region_limits,
                 aes(xintercept = length),
                 linetype = "dashed", alpha = 0.5) +
      facet_wrap(~region,ncol = 1, strip.position = "right", scales = "free_y") +
      scale_y_continuous(breaks = NULL, labels = NULL, name = "Read density") +
      scale_x_continuous(limits = c(0, 500), name = "Length (bp)")
  },
  
  heattree = {
    set.seed(3)
    theme_update(panel.border = element_blank())
    taxdata %>%
      metacoder::heat_tree(
        layout = "davidson-harel",
        # initial_layout = "reingold-tilford",
        node_size = rowMeans(.$data$tax_asv[,-1]),
        node_color = rowMeans(.$data$tax_asv[,-1]),
        node_label = replace_na(taxon_names, "unknown"),
        node_size_range = c(0.002, .035),
        node_size_axis_label = "ASV count",
        node_color_axis_label = "ASV count",
        node_label_size_range = c(0.01, 0.03),
        edge_size = rowMeans(.$data$tax_read[,-1]),
        edge_size_trans = "linear",
        edge_size_axis_label = "Read abundance",
        edge_size_range = c(0.001, 0.03),
        edge_color = rowMeans(.$data$tax_read[,-1]),
        edge_color_axis_label = "Read abundance",
        # aspect_ratio = 3/2,
        make_node_legend = FALSE,
        make_edge_legend = FALSE,
        output_file = file_out("temp/heattree.pdf")
      )
  },
  
  heattree_length_compare = {
    taxdata2 <- taxdata
    val_cols = names(taxdata2$data$tax_read) %>%
      keep(str_detect, "reads_.+_.+_.+_Consensus")
    
    taxdata2$filter_obs(
      "tax_read",
      rowSums(select_at(taxdata2$data$tax_read, val_cols)) > 0,
      drop_taxa = TRUE,
      drop_obs = TRUE,
      supertaxa = TRUE,
      reassign_obs = FALSE) %>%
      invisible()
    
    taxdata2$data$diff_read <- metacoder::compare_groups(
      taxdata2,
      data = "tax_read",
      cols = val_cols,
      groups = val_cols %>%
        str_replace("reads_.+_(.+)_.+_.+", "\\1"),
      func = function(abund1, abund2) {
        list(read_ratio = log10(mean(abund1) / mean(abund2)))
      }
    )
    taxdata2$data$diff_asv <- metacoder::compare_groups(
      taxdata2,
      data = "tax_asv",
      cols = names(taxdata2$data$tax_asv) %>% keep(startsWith, "ASVs_"),
      groups = names(taxdata2$data$tax_asv) %>%
        keep(startsWith, "ASVs_") %>%
        str_replace("ASVs_.+_(.+)_.+_.+", "\\1"),
      func = function(abund1, abund2) {
        list(asv_ratio = log10(mean(abund1) / mean(abund2)))
      }
    )
    
    set.seed(3)
    theme_update(panel.border = element_blank())
    taxdata2 %>%
      metacoder::heat_tree(
        layout = "davidson-harel",
        # initial_layout = "reingold-tilford",
        
        node_size = rowMeans(.$data$tax_asv[,-1]),
        node_size_range = c(0.002, .035),
        node_size_axis_label = "ASV richness",
        
        node_color = asv_ratio,
        node_color_range = metacoder::diverging_palette(),
        node_color_trans = "linear",
        node_color_axis_label = "Log richness ratio",
        node_color_interval = c(-1.5, 1.5),
        
        edge_size = rowMeans(.$data$tax_read[,-1]),
        edge_size_range = c(0.001, 0.03),
        edge_size_trans = "linear",
        edge_size_axis_label = "Read abundance",
        
        edge_color = read_ratio,
        edge_color_range = metacoder::diverging_palette(),
        edge_color_trans = "linear",
        edge_color_axis_label = "Log abundance ratio",
        edge_color_interval = c(-3, 3),
        
        node_label = ifelse(
          pmax(abs(.$data$diff_read$read_ratio) > 0.5,
               abs(.$data$diff_asv$asv_ratio) > 0.25),
          taxon_names,
          ""
        ),
        # aspect_ratio = 3/2,
        node_label_size_range = c(0.010, 0.03),
        output_file = file_out("temp/heattree_amplicons.pdf")
        # title = unique(paste(.$data$diff_read$treatment_1, "vs.", .$data$diff_read$treatment_2))
      )
  },
  
  ecm_heattree = {
    taxdata2 <- taxdata_ECM
    
    val_cols <- names(taxdata2$data$tax_read) %>%
      keep(str_detect, "reads_.+_.+_.+_Consensus")
    
    taxdata2$filter_obs(
      "tax_read",
      rowSums(select_at(taxdata2$data$tax_read, val_cols)) > 0,
      drop_taxa = TRUE,
      drop_obs = TRUE,
      supertaxa = TRUE,
      reassign_obs = FALSE) %>%
      invisible()
    
    taxdata2$data$diff_read <- metacoder::compare_groups(
      taxdata2,
      data = "tax_read",
      cols = val_cols,
      groups = val_cols %>%
        str_replace("reads_.+_(.+)_.+_Consensus", "\\1"),
      func = function(abund1, abund2) {
        list(read_ratio = log10(mean(abund1) / mean(abund2)))
      }
    )
    taxdata2$data$diff_asv <- metacoder::compare_groups(
      taxdata2,
      data = "tax_asv",
      cols = names(taxdata2$data$tax_asv) %>% keep(startsWith, "ASVs_"),
      groups = names(taxdata2$data$tax_asv) %>%
        keep(startsWith, "ASVs_") %>%
        str_replace("ASVs_.+_(.+)_.+_.+", "\\1"),
      func = function(abund1, abund2) {
        list(ASV_ratio = log10(mean(abund1) / mean(abund2)))
      }
    )
    
    taxdata2 <- taxdata2$filter_taxa(!is.nan(ASV_ratio), !is.nan(read_ratio))
    
    set.seed(3)
    theme_update(panel.border = element_blank())
    taxdata2 %>%
      metacoder::heat_tree(
        layout = "davidson-harel",
        initial_layout = "reingold-tilford",
        
        node_size = rowMeans(.$data$tax_asv[,-1]),
        node_size_range = c(0.002, .035),
        node_size_axis_label = "ASV richness",
        
        node_color = ASV_ratio,
        node_color_range = metacoder::diverging_palette(),
        node_color_trans = "linear",
        node_color_axis_label = "Log richness ratio",
        node_color_interval = c(-0.9, 0.9),
        
        edge_size = rowMeans(.$data$tax_read[,-1]),
        edge_size_range = c(0.001, 0.03),
        edge_size_trans = "linear",
        edge_size_axis_label = "Read abundance",
        
        edge_color = read_ratio,
        edge_color_range = metacoder::diverging_palette(),
        edge_color_trans = "linear",
        edge_color_axis_label = "Log abundance ratio",
        edge_color_interval = c(-1.8, 1.8),
        
        node_label = taxon_names,
        node_label_size_range = c(0.02, 0.03),
        # title = unique(paste(.$data$diff_read$treatment_1, "vs.", .$data$diff_read$treatment_2)),
        aspect_ratio = 3/2,
        output_file = file_out("temp/ecm_heattree.pdf")
        
      )
  },
  
  buffer_compare_physeq = {
    physeq <- proto_physeq %>%
      phyloseq::`sample_data<-`(inset2(
        phyloseq::sample_data(.),
        "qual",
        value = fct_collapse(phyloseq::sample_data(.)$qual, "-" = c("", "*"))
      )) %>%
      phyloseq::subset_samples(!is.na(site)) %>%
      phyloseq::merge_samples(group = glue::glue_data(
        phyloseq::sample_data(.),
        "{year}_{buffer}_{tech}_{amplicon}"
      ))
    
    ranks <- c("root", "kingdom", "phylum", "class", "order", "family", "genus")
    taxon_reads %>%
      mutate(root = "Root") %>%
      group_by(label) %>%
      arrange(amplicon, reference, Algorithm, .by_group = TRUE) %>%
      summarize_at(ranks, first) %>%
      left_join(
        tibble(label = phyloseq::taxa_names(physeq)),
        .,
        by = "label"
      ) %>% {
        split(select(., -label), .$label)
      } %>%
      lapply(unlist) %>%
      phyloseq::build_tax_table() %>%
      phyloseq::`tax_table<-`(physeq, .) %>%
      phyloseq::tax_glom(taxrank = "family", NArm = FALSE)
  },
  
  buffer_compare_taxmap = {
    taxmap <- metacoder::parse_phyloseq(buffer_compare_physeq)
    taxmap$data$otu_table <- metacoder::calc_obs_props(taxmap, "otu_table")
    taxmap$data$read_count <- metacoder::calc_taxon_abund(taxmap, "otu_table")
    taxmap
  },
  
  buffer_compare_long = {
    taxmap <- buffer_compare_taxmap
    value_cols <- names(taxmap$data$read_count) %>% keep(endsWith, "Long")
    taxmap$data$read_diff <- metacoder::compare_groups(
      taxmap,
      "read_count",
      cols = value_cols,
      groups = value_cols %>%
        str_match("(2015|2016)_(Xpedition|LifeGuard)") %>% {
          paste(.[,3], .[,2])
        },
      func = function(abund1, abund2) {
        list(read_ratio = log10(mean(abund1) / mean(abund2)) %>% ifelse(is.nan(.), 0, .))
      },
      combinations = list(
        c("Xpedition 2016", "Xpedition 2015"),
        c("Xpedition 2016", "LifeGuard 2016"),
        c("Xpedition 2015", "LifeGuard 2016")
      )
    )
    theme_update(panel.border = element_blank())
    taxa::filter_taxa(
      taxmap,
      taxmap$data$read_count %>% select_at(value_cols) %>% apply(MARGIN = 1, FUN = max) > 0.001,
      drop_obs = TRUE,
      reassign_obs = FALSE
    ) %>%
      metacoder::heat_tree_matrix(
        "read_diff",
        node_label = taxon_names,
        node_color = read_ratio,
        node_size = rowMeans(.$data$read_count[-1]),
        node_size_range = c(0.0001, 0.03),
        node_size_axis_label = "Read abundance",
        node_label_size_range = c(0.005, 0.03),
        node_color_range = metacoder::diverging_palette(),
        node_color_axis_label = "Log abundance ratio",
        node_color_interval = c(-3, 3),
        node_color_trans = "linear",
        # edge_size = rowMeans(.$data$read_count[-1]),
        # edge_size_axis_label = "Read abundance",
        # edge_size_range = c(0.001, 0.03),
        # edge_size_trans = "linear",
        layout = "davidson-harel",
        initial_layout = "reingold-tilford",
        output_file = file_out("temp/compare_long.pdf")
      )
  },
  
  buffer_compare_short = {
    taxmap <- buffer_compare_taxmap
    value_cols <- names(taxmap$data$read_count) %>% keep(endsWith, "Short")
    taxmap$data$read_diff <- metacoder::compare_groups(
      taxmap,
      "read_count",
      cols = value_cols,
      groups = value_cols %>%
        str_match("(2015|2016)_(Xpedition|LifeGuard)") %>% {
          paste(.[,3], .[,2])
        },
      func = function(abund1, abund2) {
        list(read_ratio = log10(mean(abund1) / mean(abund2)) %>% ifelse(is.nan(.), 0, .))
      },
      combinations = list(
        c("Xpedition 2016", "Xpedition 2015"),
        c("Xpedition 2016", "LifeGuard 2016"),
        c("Xpedition 2015", "LifeGuard 2016")
      )
    )
    theme_update(panel.border = element_blank())
    taxa::filter_taxa(
      taxmap,
      taxmap$data$read_count %>% select_at(value_cols) %>% apply(MARGIN = 1, FUN = max) > 0.001,
      drop_obs = TRUE,
      reassign_obs = FALSE
    ) %>%
      metacoder::heat_tree_matrix(
        "read_diff",
        node_label = taxon_names,
        node_color = read_ratio,
        node_size = rowMeans(.$data$read_count[-1]),
        node_size_range = c(0.0001, 0.03),
        node_size_axis_label = "Read abundance",
        node_label_size_range = c(0.005, 0.03),
        node_color_range = metacoder::diverging_palette(),
        node_color_axis_label = "Log abundance ratio",
        node_color_interval = c(-3, 3),
        node_color_trans = "linear",
        # edge_size = rowMeans(.$data$read_count[-1]),
        # edge_size_axis_label = "Read abundance",
        # edge_size_range = c(0.001, 0.03),
        # edge_size_trans = "linear",
        layout = "davidson-harel",
        initial_layout = "reingold-tilford",
        output_file = file_out("temp/compare_short.pdf")
      )
  },
  
  # Generate ASV table for comparisons between buffer solutions
  buffer_asv_table = 
    proto_physeq %>%
    # No Ion Torrent samples, because they are not completely demultiplexed
    phyloseq::subset_samples(tech != "Ion Torrent") %>%
    # Don't use QC samples
    phyloseq::subset_samples(sample_type == "Sample") %>%
    # Only include fungi
    phyloseq::prune_taxa(taxa = fungi_PHYLOTAX.Cons$label) %>%
    # Only use samples with more than 100 reads
    phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 100) %>%
    # Only use samples where there reads from all three tech×amplicon combinations
    phyloseq::prune_samples(
      samples = phyloseq::sample_data(.) %>%
        as("data.frame") %>%
        as_tibble(rownames = "sample") %>%
        group_by(site, x, tech, amplicon) %>%
        filter(n() >= 3) %>%
        pull(sample) %>%
        unique()
    ),
  
  # Create distance matrix for the buffer solution comparison
  buffer_asv_bc_dist =
    buffer_asv_table %>%
    # Normalize reads to 1
    phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%
    # Bray-Curtis distance
    phyloseq::distance("bray"),
  
  # Calculate PERMANOVA
  buffer_adonis = target(
    vegan::adonis2(
      buffer_asv_bc_dist ~ paste(site, x) + paste(tech, amplicon) + buffer + year,
      data = as(phyloseq::sample_data(buffer_asv_table), "data.frame"),
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
      data = as(phyloseq::sample_data(buffer_asv_table), "data.frame")
    ),
  
  # Generate ASV table for comparisons between sequencing methods
  tech_asv_table =
    proto_physeq %>%
    # Don't use Ion Torrent because it didn't multiplex properly
    phyloseq::subset_samples(tech != "Ion Torrent") %>%
    # Only use the Xpedition buffer
    phyloseq::subset_samples(buffer == "Xpedition") %>%
    # Don't use QC sample
    phyloseq::subset_samples(sample_type == "Sample") %>%
    # Only fungi
    phyloseq::prune_taxa(taxa = fungi_PHYLOTAX.Cons$label),
  
  # Add taxonomy to the ASV table and cluster at the class level
  tech_class_table =
    phyloseq::`tax_table<-`(
      tech_asv_table,
      phyloseq::tax_table(
        fungi_PHYLOTAX.Cons %>%
          select(label, kingdom:genus) %>%
          column_to_rownames("label") %>%
          as.matrix()
      )
    ) %>%
    phyloseq::tax_glom(taxrank = "class") %>%
    phyloseq::`taxa_names<-`(phyloseq::tax_table(.)[,"class"]) %>%
    # Only samples with at least 100 reads
    phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 100) %>%
    # Only soil samples with results from all three tech/amplicon combinations
    phyloseq::prune_samples(
      samples = phyloseq::sample_data(.) %>%
        as("data.frame") %>%
        as_tibble(rownames = "sample") %>%
        group_by(site, x, year) %>%
        filter(n() >= 3) %>%
        pull(sample) %>%
        unique()
    ) %>%
    # Normalize reads
    phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%
    # Take only classes which represent at least 1% of reads in at least one sample.
    phyloseq::prune_taxa(
      taxa = (
        phyloseq::otu_table(.) %>%
          t() %>%
          as.data.frame() %>%
          do.call(what = pmax)
      ) > 0.01
    ),
  
  # Sample data for easy reference
  tech_class_sample_data =
    as(phyloseq::sample_data(tech_class_table), "data.frame"),
  
  # Calculate Bray-Curtis distance for comparing technologies
  tech_class_bc_dist = phyloseq::distance(tech_class_table, "bray"),
  
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
      phyloseq::otu_table(tech_class_table) ~ Condition(paste(site, x, year)),
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
  tech_class_pcoa_plot = 
    left_join(
      vegan::scores(tech_class_pcoa, display = "sites") %>%
        as_tibble(rownames = "sample"),
      tech_class_sample_data %>%
        as_tibble(rownames = "sample"),
      by = "sample"
    ) %>%
    mutate(
      group = paste(tech, "--", amplicon) %>%
        factor(
          levels = c("Illumina -- Short", "PacBio -- Short", "PacBio -- Long"),
          ordered = TRUE
        )
    ) %>%
    ggplot(aes(x = MDS1, y = MDS2, color = group, pch = group)) +
    geom_vline(xintercept = 0, color = "gray80") +
    geom_hline(yintercept = 0, color = "gray80") +
    geom_point() +
    geom_text(
      aes(x = MDS1, y = MDS2, label = abbrev),
      data = tech_class_pcoa_scores,
      inherit.aes = FALSE
    ) +
    scale_color_discrete(name = NULL) +
    scale_discrete_manual(aesthetics = "pch", name = NULL, values = 1:3) +
    coord_equal(),
  
  # Generate ASV table for comparisons between sequencing methods for ECM
  tech_ecm_asv_table =
    tech_asv_table %>%
    # Only ECM
    phyloseq::prune_taxa(taxa = ecm_PHYLOTAX.Cons$label),
  
  # Add taxonomy to the ECM ASV table and cluster at the family level
  tech_ecm_fam_table =
    phyloseq::`tax_table<-`(
      tech_asv_table,
      phyloseq::tax_table(
        ecm_PHYLOTAX.Cons %>%
          select(label, kingdom:genus) %>%
          column_to_rownames("label") %>%
          as.matrix()
      )
    ) %>%
    phyloseq::tax_glom(taxrank = "family") %>%
    phyloseq::`taxa_names<-`(phyloseq::tax_table(.)[,"family"]) %>%
    # Only samples with at least 100 reads
    phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 100) %>%
    # Only soil samples with results from all three tech/amplicon combinations
    phyloseq::prune_samples(
      samples = phyloseq::sample_data(.) %>%
        as("data.frame") %>%
        as_tibble(rownames = "sample") %>%
        group_by(site, x, year) %>%
        filter(n() >= 3) %>%
        pull(sample) %>%
        unique()
    ) %>%
    # Normalize reads
    phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%
    # Take only classes which represent at least 1% of reads in at least one sample.
    phyloseq::prune_taxa(
      taxa = (
        phyloseq::otu_table(.) %>%
          t() %>%
          as.data.frame() %>%
          do.call(what = pmax)
      ) > 0.01
    ),
  
  # Sample data for easy reference
  tech_ecm_fam_sample_data =
    as(phyloseq::sample_data(tech_ecm_fam_table), "data.frame"),
  
  # Calculate Bray-Curtis distance for comparing technologies
  tech_ecm_fam_bc_dist = phyloseq::distance(tech_ecm_fam_table, "bray"),
  
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
      phyloseq::otu_table(tech_ecm_fam_table) ~ Condition(paste(site, x, year)),
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
    ggplot(aes(x = MDS2, y = MDS1, color = group, pch = group)) +
    geom_vline(xintercept = 0, color = "gray80") +
    geom_hline(yintercept = 0, color = "gray80") +
    geom_point() +
    geom_text(
      aes(x = MDS2 * 3, y = MDS1 * 3, label = abbrev),
      data = tech_ecm_fam_pcoa_scores,
      inherit.aes = FALSE
    ) +
    scale_color_discrete(name = NULL) +
    scale_discrete_manual(aesthetics = "pch", name = NULL, values = 1:3) +
    coord_equal(),
  
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
    set_names(names(.) %>% str_extract("(PacBio|Illumina|Ion.Torrent)_(Long|Short)")) %>%
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
    set_names(names(.) %>% str_extract("(PacBio|Illumina|Ion.Torrent)_(Long|Short)")) %>%
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
  
  # short_length_glm = 
  #   short_length_table %>%
  #   glm(reads ~ seq_run + seq_run:length - 1, family = "poisson", data = .),
  
  length_model_data =
    reads_table %>%
    filter(region == "short") %>%
    mutate_at("seq", tzara::seqhash) %>%
    mutate_at("seq", factor) %>%
    tidyr::complete(seq_run, nesting(seq, length), fill = list(reads = 0)),
  
  # length_glmm =
  #   MASS::glmmPQL(
  #     fixed = reads ~ seq_run * length,
  #     random = ~ 1 | seq,
  #     data = length_model_data,
  #     family = "poisson"
  #   ),
  
  # length_glmer =
  #   lme4::glmer(
  #     formula = reads ~ seq_run * length + (1 | seq),
  #     data = length_model_data,
  #     family = "poisson"
  #   ),
  
  # confint_glmer = 
  #   lme4::confint.merMod(length_glmer),
  
  # length_glmernb =
  #   length_model_data %>%
  #   lme4::glmer.nb(
  #     formula = reads ~ seq_run * length + (1 | seq),
  #     data = .
  #   ),
  
  # confint_glmernb = 
  #   lme4::confint.merMod(length_glmernb),
  
  # length_mcmcglmm =
  #   MCMCglmm::MCMCglmm(
  #     fixed = reads ~ trait:seq_run + trait:seq_run:length - 1,
  #     random = ~ idh(trait):seq,
  #     data = length_model_data,
  #     family = "zipoisson",
  #     rcov = ~us(trait):units
  #   ),
  # 
  # length_mcmcglmm2 =
  #   MCMCglmm::MCMCglmm(
  #     fixed = reads ~ seq + seq_run + seq_run:length - 1,
  #     data = short_otu_table,
  #     family = "poisson",
  #     nitt = 103000,
  #     thin = 100,
  #     burnin = 3000
  #   ),
  
  # write *.aux and *.tex files.  Remove extra bibliographies from *.tex files
  supplement_tex =
    withr::with_dir(
      "writing",
      {
        file_in(!!file.path("writing", "all.bib"))
        knitr_in(!!file.path("writing", "transect_supplement.Rmd"))
        file_in(!!file.path("writing", "preamble_supplement.tex"))
        file_in(!!file.path("scripts", "output_functions.R"))
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
      file_in(!!file.path("scripts", "output_functions.R"))
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
  
  supplement_diff_tex = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_supplement_original.tex"))
      file_in(!!file.path("writing", "transect_supplement.tex"))
      file_out(!!file.path("writing", "transect_supplement_diff.tex"))
      tinytex::tlmgr_install("latexdiff")
      system2(
        command = "latexdiff",
        args = list(
          "transect_supplement_original.tex",
          "transect_supplement.tex"
        ),
        stdout = TRUE
      ) %>%
        stringr::str_replace_all(
          "transect_paper.tex",
          "transect_paper_diff.tex"
        ) %>%
        writeLines("transect_supplement_diff.tex")
    }
  ),
  
  article_diff_tex = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_paper_original.tex"))
      file_in(!!file.path("writing", "transect_paper.tex"))
      file_out(!!file.path("writing", "transect_paper_diff.tex"))
      tinytex::tlmgr_install("latexdiff")
      system2(
        command = "latexdiff",
        args = list(
          "transect_paper_original.tex",
          "transect_paper.tex"
        ),
        stdout = TRUE
      ) %>%
        stringr::str_replace_all(
          "transect_supplement.tex",
          "transect_supplement_diff.tex"
        ) %>%
        writeLines("transect_paper_diff.tex")
    }
  ),
  
  supplement_pdf = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_supplement.tex"))
      file_in(!!file.path("writing", "transect_paper.tex"))
      file_out(!!file.path("writing", "transect_supplement.pdf"))
      tinytex::lualatex("transect_supplement.tex", bib_engine = "biber")
    }
  ),
  
  supplement_diff_pdf = withr::with_dir(
    "writing",
    {
      file_in(!!file.path("writing", "transect_supplement_diff.tex"))
      file_in(!!file.path("writing", "transect_paper_diff.tex"))
      file_out(!!file.path("writing", "transect_supplement_diff.pdf"))
      tinytex::lualatex("transect_supplement_diff.tex", bib_engine = "biber")
    }
  ),
  
  article_pdf = withr::with_dir(
    "writing", 
    {
      file_in(!!file.path("writing", "transect_paper.tex"))
      file_in(!!file.path("writing", "transect_supplement.tex"))
      file_out(!!file.path("writing", "transect_paper.pdf"))
      tinytex::lualatex("transect_paper.tex", bib_engine = "biber")
    }
  ),
  
  article_diff_pdf = withr::with_dir(
    "writing", 
    {
      file_in(!!file.path("writing", "transect_paper_diff.tex"))
      file_in(!!file.path("writing", "transect_supplement_diff.tex"))
      file_out(!!file.path("writing", "transect_paper_diff.pdf"))
      tinytex::lualatex("transect_paper_diff.tex", bib_engine = "biber")
    }
  ),
  
  tree_fig = {
    
    charscale <- 0.015
    rankgap <- 0.015
    tip_taxa <- 
      phylotaxon_decipher_unconst_long_fungi$tip_taxa %>%
      filter(rank != "kingdom") %>%
      group_by(label, rank) %>%
      filter(if (any("phylotax" == method)) method == "phylotax" else TRUE) %>%
      group_by(label, rank) %>%
      summarize(
        taxon =  table(taxon) %>%
          paste0(names(.), collapse = "/") %>%
          gsub(pattern = "(.+/.+)", replacement = "<\\1>") %>%
          gsub(pattern = "^\\d+", replacement = "")
      ) %>% inner_join(
        phylo_labeled_tree_decipher_unconst_long_fungi$tip.label %>%
          gsub(pattern = "\"([a-f0-9]{8}).*", replacement = "\\1") %>%
          enframe(name = "tip", value = "label"),
        by = "label"
      ) %>%
      mutate(
        depth = ape::node.depth.edgelength(
          phylo_labeled_tree_decipher_unconst_long_fungi
        )[tip]
      ) %>%
      group_by(label) %>%
      arrange(desc(rank), .by_group = TRUE) %>%
      mutate(width = nchar(taxon) * charscale + rankgap,
             offset = cumsum(lag(width, default = 0)))
    clade_annot <- list()
    ranks <- c("genus", "family", "order", "class", "phylum")
    for (i in seq_along(ranks)) {
      rank_nodes <- phylotaxon_decipher_unconst_long_fungi$node_taxa %>%
        filter(rank == ranks[i]) %>%
        group_by(taxon) %>%
        mutate(color = ifelse(n() > 1, "tomato", "black")) %>%
        ungroup() %>%
        mutate(tip = phangorn::Descendants(
          phylo_labeled_tree_decipher_unconst_long_fungi,
          node = node,
          type = "tips"
        )) %>%
        unnest(tip) %>%
        left_join(tip_taxa, by = c("rank", "tip", "taxon")) %>%
        group_by(node, taxon, color) %>%
        summarize(offset = max(depth + offset) - max(depth)) %>%
        rename(label = taxon) %>%
        bind_rows(
          tip_taxa %>%
            ungroup() %>%
            filter(rank == ranks[i], startsWith(taxon, "<")) %>%
            select(label = taxon, node = tip, offset) %>%
            mutate(color = "steelblue")
        )
      clade_annot <- c(
        clade_annot,
        pmap(rank_nodes,
             ggtree::geom_cladelabel,
             align = FALSE,
             fontsize = 1.5,
             extend = 0.3)
        )
    }
    (phylo_labeled_tree_decipher_unconst_long_fungi %>%
        ggtree::ggtree(
          aes(size = if_else(
            !is.na(as.integer(label)) & as.integer(label) > 90,
            0.5,
            0.2
          )),
          lineend = "butt",
          linejoin = "mitre"
        ) %>%
        reduce(clade_annot, `+`, .init = .) + 
        scale_size_identity() +
        xlim(0, max(tip_taxa$offset + tip_taxa$width + tip_taxa$depth))) %T>%
      ggsave(
        filename = file_out(!!file.path("output", "supp_file_3.pdf")),
        plot = .,
        device = "pdf",
        width = 8,
        height = 20,
        units = "in"
      )
  },
  
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
  
  #### prepare submission materials for ENA ####
  soil_samples =
    readd(proto_physeq) %>%
      phyloseq::sample_data() %>%
      as("data.frame") %>%
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
    ),
    ENA_soil_samples =  {
      template_in <- file_in("config/ENA_soil_sample_template.tsv")
      template_names <- names(read_tsv(template_in, skip = 2))
      tsv_out <- file_out("output/ENA/ENA_soil_samples.tsv")
      if (!dir.exists(dirname(tsv_out))) dir.create(dirname(tsv_out), recursive = TRUE)
      file.copy(template_in, tsv_out, overwrite = TRUE)
      soil_samples %>%
        assertr::verify(template_names %in% names(.)) %>%
        select_at(template_names) %>%
        write_tsv(tsv_out, append = TRUE, col_names = FALSE)
  },
  control_samples =
    readd(proto_physeq) %>%
    phyloseq::sample_data() %>%
    as("data.frame") %>%
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
    ),
  ENA_control_samples = {
    template_in <- file_in("config/ENA_control_sample_template.tsv")
    template_names <- names(read_tsv(template_in, skip = 2))
    tsv_out <- file_out("output/ENA/ENA_control_samples.tsv")
    if (!dir.exists(dirname(tsv_out))) dir.create(dirname(tsv_out), recursive = TRUE)
    file.copy(template_in, tsv_out, overwrite = TRUE)
    control_samples %>%
      assertr::verify(template_names %in% names(.)) %>%
      select_at(template_names) %>%
      write_tsv(tsv_out, append = TRUE, col_names = FALSE)
  },
  blank_samples =
    readd(proto_physeq) %>%
    phyloseq::sample_data() %>%
    as("data.frame") %>%
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
    ),
  ENA_blank_samples = {
    template_in <- file_in("config/ENA_control_sample_template.tsv")
    template_names <- names(read_tsv(template_in, skip = 2))
    tsv_out <- file_out("output/ENA/ENA_blank_samples.tsv")
    if (!dir.exists(dirname(tsv_out))) dir.create(dirname(tsv_out), recursive = TRUE)
    file.copy(template_in, tsv_out, overwrite = TRUE)
    blank_samples %>%
      assertr::verify(template_names %in% names(.)) %>%
      select_at(template_names) %>%
      write_tsv(tsv_out, append = TRUE, col_names = FALSE)
  },
  pacbio_ion_manifest = 
    dplyr::inner_join(
      dplyr::bind_rows(
        soil_samples,
        control_samples,
        blank_samples
      ) %>% select(sample_alias, "project name", plate, well, amplicon, sample_type),
      dplyr::bind_rows(
        file_meta_gits7its4 %>%
          mutate(plate = ifelse(tech == "Ion Torrent", list(c("001", "002")), as.list(plate))) %>%
          unnest(plate),
        readd(file_meta_its1lr5)
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
          d$DESCRIPTION <- paste("also contains samples from", d$SAMPLE[2], "due to incomplete multiplexing")
          d$SAMPLE <- d$SAMPLE[1]
          d$NAME <- d$NAME[1]
          d$LIBRARY_SOURCE <- d$LIBRARY_SOURCE[1]
        }
        if (length(d$NAME) > 1) print(d)
        sink(file.path("output", "ENA", paste0(d$NAME, ".manifest")), append = FALSE)
        iwalk(d, ~ for (x in .x) {cat(.y, x); cat("\n")})
        sink()
      }
    ),
  
  illumina_manifest = dplyr::inner_join(
    dplyr::bind_rows(
      soil_samples,
      control_samples,
      blank_samples
    ) %>%
      select(sample_alias, "project name", plate, well, amplicon, sample_type),
    illumina_group_SH.2257,
    by = c("plate", "well", "amplicon")
  ) %>% 
    dplyr::transmute(
      STUDY = `project name`,
      SAMPLE = sample_alias,
      NAME = paste(sample_alias, tech, direction, sep = "-") %>% str_replace("-$", "") %>% str_replace(" ", ""),
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
    ),
  read_accnos =
    list.files(path = file_in("output/ENA/reads"), pattern = "receipt.xml", recursive = TRUE) %>%
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
    tidyr::pivot_wider(names_from = "tech", values_from = "accno"),
  pretaxid =
    preguild_taxa_short_euk_PHYLOTAX.Cons %>%
    dplyr::select(kingdom:genus, Taxonomy, label) %>%
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
    ),
  taxid = target(
    do.call(lookup_ncbi_taxon, pretaxid),
    dynamic = map(pretaxid),
    retries = 2
  ),
  taxid_all = bind_rows(readd(taxid)) %T>%
    write_csv(file_out("output/taxa.csv")) %T>% {
      dplyr::select(., label, targetTaxonomy, targetRank) %>%
        write_csv(file_out("output/taxa_target.csv"))
    } %T>% {
      dplyr::select(., label, Taxonomy, rank) %>%
        write_csv(file_out("output/taxa_ncbi.csv"))
    },
  pretaxid_env =
    preguild_taxa_short_euk_PHYLOTAX.Cons %>%
    dplyr::select(kingdom:order, label) %>%
    dplyr::group_by_at(vars(kingdom:order)) %>%
    dplyr::summarize(label = paste(label, collapse = ",")) %>% 
    tidyr::pivot_longer(cols = kingdom:order, names_to = "rank", values_to = "taxon") %>%
    dplyr::filter(complete.cases(.), !stringr::str_detect(taxon, "\\d")) %>%
    dplyr::group_by(label) %>%
    dplyr::summarize(
      Taxonomy = paste(taxon, collapse = ";"),
      taxon = dplyr::last(taxon),
      rank = dplyr::last(rank)
    ) %>%
    dplyr::group_by(Taxonomy, taxon, rank) %>%
    dplyr::summarize(label = paste(label, collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_at(
      "taxon",
      stringi::stri_replace_all_regex,
      c("Fungi", "Alveolata", "Gregarinasina", "Neogregarinorida", "NA", "Angiospermae", "Stramenopila",
        "Kickxellomycota", "Zoopagomycota", "Monoblepharomycota", "Morteriellomycota", "Bryopsida",
        "Chromerida", "Gigasporales", "Paraglomeromycetes", "Ichthyosporia"),
      c("fungus", "alveolate", "gregarine", "gregarine", "eukaryote", "Magnoliophyta", "Stramenopile",
        "Kickxellomycotina", "Zoopagomycotina", "Monoblepharidomycetes", "Morteriellomycotina", "Bryophyta",
        "Colpodellida", "Diversisporales", "Paraglomerales", "Ichthyosporea"),
      vectorize_all = FALSE
    ) %>%
    dplyr::mutate(taxon = ifelse(
      grepl("(Streptophyta|Metazoa)", Taxonomy),
      paste(taxon, "environmental sample"),
      paste("uncultured", taxon)
    )),
  taxid_env = target(
    do.call(lookup_ncbi_taxon, pretaxid_env),
    dynamic = map(pretaxid_env),
    retries = 2
  ),
  taxid_env_all = bind_rows(readd(taxid_env)) %T>%
    write_csv(file_out("output/env_taxa.csv")) %T>% {
      dplyr::select(., targetTaxonomy) %>%
        write_csv(file_out("output/env_taxa_target.csv"))
    } %T>% {
      dplyr::select(., Taxonomy) %>%
        write_csv(file_out("output/env_taxa_ncbi.csv"))
    },
  trace = TRUE
) %>%
  dplyr::filter(ifelse(amplicon == '"Short"', metric != '"wunifrac"', TRUE) %|% TRUE)

saveRDS(plan2, "data/plan/drake_light.rds")
njobs <- max(local_cpus() %/% 2, 1)
parallelism <- if (njobs > 1) "clustermq" else "loop"
options(clustermq.scheduler = "multicore")
#if (interactive()) {
cache_dir <- ".drake"
cache <- drake_cache(cache_dir)
# dconfig <- drake_config(plan2, cache = cache)
# vis_drake_graph(dconfig)
make(
  plan2,
  cache = cache,
  parallelism = parallelism,
  jobs = njobs,
  # targets = "variofit2_ecm_bray_PacBio_Short_Consensus",
  lazy_load = FALSE,
  memory_strategy = "autoclean",
  garbage_collection = TRUE,
  cache_log_file = "drake_light_cache.csv",
  console_log_file = file.path(normalizePath("logs"), "drake_light.log")
)
#}
