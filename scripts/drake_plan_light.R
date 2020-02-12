remove(list = ls())
if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("data/plan/drake.Rdata")
}
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

source(file.path(config$rdir, "dada.R"))
source(file.path(config$rdir, "mantel.R"))
source(file.path(config$rdir, "variogram.R"))
source(file.path(config$rdir, "qstats.R"))
source(file.path(config$rdir, "taxonomy.R"))
source(file.path(config$rdir, "plate_check.R"))
source(file.path(config$rdir, "output_functions.R"))

choosevars <- function(d, g, .data) {
  .data %>%
    filter(type == g$type) %>%
    select(type, x = !!paste(g$x_var, "-", g$x_amplicon), y = !!paste(g$y_var, "-", g$y_amplicon)) %>%
    filter(x > 0 | y > 0) %>%
    mutate(
      x_var = g$x_var, 
      x_amplicon = g$x_amplicon,
      y_var = g$y_var,
      y_amplicon = g$y_amplicon
    )
}

datasets <- read_csv(config$dataset, col_types = "cccccccicccccicc")
regions <- read_csv(config$regions, col_types = "cccciiiiic")

# combo_meta <- tibble(
#   group = c("fungi", "euk"),
#   taxname = "short",
#   reftree = c("fungi_tree_decipher_unconst_long", "tree_decipher_unconst_long"),
#   conf_taxon = glue::glue("conf_taxon_short_{group}"),
#   phylotaxon = glue::glue("phylotaxon_unconst_long_{group}")
# ) %>%
#   mutate_at(c("reftree", "conf_taxon", "phylotaxon"), syms)

guilds_meta <- tibble(
  consensus_taxa = c("phylotaxon_decipher_unconst_long_fungi", "conf_taxon_short_fungi", "best_taxon_short_fungi"),
  amplicon = c("Long", "Short", "Short"),
  algorithm = c("PHYLOTAX", "Consensus", "PHYLOTAX+Cons")
) %>%
  mutate_at("consensus_taxa", syms)

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
    ecm = glue::glue("ecm_{algorithm}"),
    ecm2 = glue::glue("ecm2_{algorithm}"),
    ecm3 = glue::glue("ecm3_{algorithm}")
  ) %>%
  mutate_at(vars(starts_with("ecm")), compose(syms, make.names))

plan2 <- drake_plan(
  # targets which are imported from first plan
  allseqs = target(trigger = trigger(mode = "blacklist")),
  taxon_table = target(trigger = trigger(mode = "blacklist")),
  taxon = target(
    transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID),
    trigger = trigger(mode = "blacklist")
  ),
  raxml_decipher_LSU = target(trigger = trigger(mode = "blacklist")),
  raxml_decipher_long = target(trigger = trigger(mode = "blacklist")),
  raxml_decipher_unconst_long = target(trigger = trigger(mode = "blacklist")),
  raxml_epa_full = target(trigger = trigger(mode = "blacklist")),
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
  qstats_n = target(trigger = trigger(mode = "blacklist")),
  qstats_length = target(trigger = trigger(mode = "blacklist")),
  qstats_eexp = target(trigger = trigger(mode = "blacklist")),
  qstats_erate = target(trigger = trigger(mode = "blacklist")),
  qstats_pnoerr = target(trigger = trigger(mode = "blacklist")),
  qstats_minq = target(trigger = trigger(mode = "blacklist")),
  
  
  # big_fasta ----
  # write the big_seq_table as a fasta files so that they can be clustered by
  # VSEARCH.
  big_fasta = target(
    write_big_fasta(big_seq_table,
                    file_out(big_fasta_file)),
    transform = map(.data = !!region_meta, .tag_in = step, .id = region)
  ),
  
  taxon_labels = make_taxon_labels(taxon_table),
  
  # funguild_db ----
  # Download the FUNGuild database
  funguild_db = FUNGuildR::get_funguild_db(),
  
  # guilds_table ----
  # Assign ecological guilds to the ASVs based on taxonomy.
  # guilds_table = target(
  #   FUNGuildR::funguild_assign(taxon_table, db = funguild_db),
  #   format = "fst"),
  
  # platemap ----
  # Read the map between plate locations and samples
  platemap = target(
    read_platemap(file_in(!!config$platemap), !!config$platemap_sheet),
    format = "fst"),
  
  tree_decipher_LSU = target(
    raxml_decipher_LSU$bipartitions,
    transform = map(
      treename = "decipher_LSU",
      group = "euk",
      .tag_in = step, .tag_out = c(euktree, tree), .id = FALSE)
  ),
  
  tree_decipher_LSU_long = target(
    raxml_decipher_long$bipartitions,
    transform = map(
      treename = "decipher_LSU_long",
      group = "euk",
      .tag_in = step, .tag_out = c(euktree, tree), .id = FALSE)
  ),
  
  tree_decipher_unconst_long = target(
    raxml_decipher_unconst_long$bipartitions,
    transform = map(
      treename = "decipher_unconst_long",
      group = "euk",
      .tag_in = step, .tag_out = c(euktree, tree), .id = FALSE)
  ),
  
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
  
  # labeled_fungi_tree = target(
  #   relabel_tree(
  #     tree = fungi_tree,
  #     old = taxon_labels$label,
  #     new = taxon_labels$tip_label
  #   ) %T>%
  #     castor::write_tree(
  #       file_out(!!paste0("data/trees/fungi_", treename, ".tree"))
  #     ),
  #   transform = map(fungi_tree, treename, .id = treename)
  # ),
  # 
  # phylolabeled_fungi_tree = target(
  #   relabel_tree(
  #     tree = fungi_tree,
  #     old = phylotaxon_labels$label,
  #     new = phylotaxon_labels$tip_label
  #   ) %T>%
  #     castor::write_tree(
  #       file_out(!!paste0("data/trees/fungi_", treename, "_phylo.tree"))
  #     ),
  #   transform = map(fungi_tree, phylotaxon_labels, treename, .id = treename)
  # ),
  # 
  phylotaxon = target(
    phylotax(tree = tree, taxa = taxon_table),
    transform = map(tree, .tag_in = step, 
                    .tag_out = consensus_taxa, .id = c(treename, group))
  ),
  
  phylotaxon_labels = target(
    make_taxon_labels(phylotaxon$tip_taxa),
    transform = map(phylotaxon, .tag_in = step, .id = c(treename, group))
  ),
  
  conf_taxon_short_euk = target(
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
      .tag_out = c(consensus_taxa, conf_taxon),
      .id = FALSE
    )
  ),
  
  conf_taxon_short_fungi = target(
    group_by(conf_taxon_short_euk$tip_taxa, label) %>%
      filter("Fungi" %in% taxon | any(endsWith(taxon, "mycota"))) %>%
      list(tip_taxa = .),
    transform = map(
      taxname = "short",
      group = "fungi",
      .tag_out = c(consensus_taxa, conf_taxon),
      .id = FALSE
    )
  ),
  
  best_taxon_short_euk = target(
    conf_taxon_short_euk$tip_taxa %>%
      filter(!label %in% tree_decipher_unconst_long$tip.label) %>%
      bind_rows(
        phylotaxon_decipher_unconst_long_euk$tip_taxa %>%
          filter(method == "phylotax")
      ) %>%
      mutate(method = "phylotax+c") %>%
      list(tip_taxa = .),
    transform = map(group = "euk",
                    taxname = "short",
                    .tag_out = consensus_taxa, .id = FALSE)
  ),
  
  best_taxon_short_fungi = target(
    conf_taxon_short_fungi$tip_taxa %>%
      filter(!label %in% fungi_tree_decipher_unconst_long$tip.label) %>%
      bind_rows(
        phylotaxon_decipher_unconst_long_fungi$tip_taxa %>%
          filter(method == "phylotax")
      ) %>%
      mutate(method = "phylotax+c") %>%
      list(tip_taxa = .),
    transform = map(group = "fungi",
                    taxname = "short",
                    .tag_out = consensus_taxa, .id = FALSE)
  ),
  
  guilds = target(
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
          dplyr::summarize(Taxonomy = paste(taxon, collapse = ";")),
        by = "label"
      ) %>%
      FUNGuildR::funguild_assign(funguild_db),
    transform = map(.data = !!guilds_meta, .tag_in = step, .id = algorithm)
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
      physeq <- proto_physeq
      if (amplicon == "Long") {
        phyloseq::phy_tree(physeq) <- fungi_tree_decipher_unconst_long
      }
      physeq <- physeq %>%
        phyloseq::subset_samples(sample_type == "Sample") %>%
        phyloseq::subset_samples(buffer == "Xpedition") %>%
        phyloseq::subset_samples(site == "Gan") %>%
        phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 100 ) %>%
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
        
      }
      physeq
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
  
  # variofit = target(
  #   gstat::fit.variogram(
  #     as_variogram(variog),
  #     gstat::vgm(variog$gamma[22], "Exp", 3, variog$gamma[22]/2),
  #     fit.method = 1
  #   ),
  #   transform = map(variog, .tag_in = step,.id = c(guild, metric, tech, amplicon, algorithm))
  # ),
  
  variofit2 = target(
    variog %>%
      # filter(!is.na(bin)) %>%
      as_variogram() %>%
      nls(
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
      ),
    transform = map(variog, .tag_in = step,.id = c(guild, metric, tech, amplicon, algorithm))
  ),
  
  variogST = target(
    variogST(physeq, metric, breaks = c(1:25 - 0.5, 30000)),
    transform = cross(physeq, metric = c("bray", "wunifrac"),
                      .tag_in = step,
                      .id = c(guild, metric, tech, amplicon, algorithm))
  ),
  
  # variofitST = target(
  #   gstat::fit.StVariogram(
  #     as_variogramST(variogST) %>%
  #       inset(,"np", ifelse(.$dist > 30, 1, .$np)),
  #     gstat::vgmST(
  #       "metric",
  #       joint = gstat::vgm(max(variogST$gamma), "Exp", 5, quantile(variogST$gamma, 0.25)),
  #       stAni = 1
  #     ),
  #     fit.method = 1
  #   ),
  #   transform = map(variogST, .tag_in = step, .id = c(guild, metric, tech, amplicon, algorithm))
  # ),
  
  variofitST2 = target(
    variogST %>%
      # filter(!is.na(bin)) %>%
      as_variogramST() %>% {
        nls(
          gamma ~ (1 - sill) - (1 - sill) * (1 - nugget) * exp((dist/range + timelag/timerange)*log(0.05)),
          # add the parameters from the spatial-only fit.
          data = .,
          start = c(
            list(timerange = 2),
            as.list(variofit2$m$getPars())
          ),
          lower = c(
            list(timerange = 0.001),
            as.list(variofit2$m$getPars())
          ),
          algorithm = "port",
          weights = pmin(.$np/.$dist, 1),
          control = list(warnOnly = TRUE, maxiter = 1000, minFactor = 1/4096)
        )
      },
    transform = map(variogST, variofit2, .tag_in = step, .id = c(guild, metric, tech, amplicon, algorithm))
  ),
  
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
  
  asv_table = big_seq_table_ITS2 %>%
    tibble::as_tibble(rownames = "filename") %>%
    tidyr::gather(key = "seq", value = "reads", -1) %>%
    dplyr::filter(reads >= 1) %>%
    tidyr::extract(col = "filename",
                   into = c("seq.run", "plate", "well", "dir", "region"),
                   regex = "([a-zA-Z]+[-_]\\d+)_(\\d+)_([A-H]1?[0-9])([fr]?)_([:alnum:]+).+") %>%
    dplyr::group_by(seq.run, seq) %>%
    dplyr::summarize(reads = sum(reads)) %>%
    tidyr::spread(key = seq.run, value = reads, fill = 0),
  
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
  
  demuxlength =  parse_qstat(qstats_length) %>%
    filter(is.na(read) | read != "_R2") %>%
    filter(is.na(step), !is.na(well), ) %>%
    group_by(seq_run) %>%
    summarize(length = reldist::wtd.quantile(length, weight = nreads)) %>%
    deframe(),
  
  demuxqual = parse_qstat(qstats_erate) %>%
    filter(is.na(read) | read != "_R2") %>%
    filter(is.na(step), !is.na(well)) %>%
    group_by(seq_run) %>%
    summarize(qual = round(weighted.mean(erate, w = nreads, na.rm = TRUE), 4)) %>%
    deframe(),
  
  readcounts = parse_qstat(qstats_n) %>% filter(is.na(read) | read != "_R2"),
  
  rawcounts = readcounts %>%
    filter(is.na(step), is.na(well)) %>%
    group_by(seq_run) %>%
    summarize(nreads =  prettyNum(sum(nreads), big.mark = " ")) %>%
    deframe(),
  
  demuxcounts = readcounts %>%
    filter(is.na(step), !is.na(well)) %>%
    group_by(seq_run) %>%
    summarize(nreads = prettyNum(sum(nreads), big.mark = " ")) %>%
    deframe(),
  
  filtercounts_full = readcounts %>%
    filter(step == "filter", is.na(region) | region %in% c("long", "short")) %>%
    group_by(seq_run, region) %>%
    summarize(nreads = prettyNum(sum(nreads), big.mark = " ")) %>%
    select(seq_run, nreads) %>%
    deframe(), 
  
  filtercounts_ITS2 = readcounts %>%
    filter(step == "filter", region == "ITS2") %>%
    group_by(seq_run, region) %>%
    summarize(nreads = prettyNum(sum(nreads), big.mark = " ")) %>%
    select(seq_run, nreads) %>%
    deframe(),
    
  regioncounts = readcounts %>%
    filter(step == "lsux") %>%
    group_by(seq_run) %>%
    summarize(nreads = prettyNum(
      sum(nreads) %/% if (seq_run[1] == "pb_500") 14 else 4,
      big.mark = " ")
    ) %>%
    deframe(),
  
  bioinf_table = bind_rows(
    enframe(rawcounts, name = "seq_run", value = "reads") %>%
      mutate(
        reads = as.integer(gsub(" ", "", reads)),
        step = "Raw"
      ),
    enframe(demuxcounts, name = "seq_run", value = "reads") %>%
      mutate(
        reads = as.integer(gsub(" ", "", reads)),
        step = "Trim"
      ),
    enframe(regioncounts, name = "seq_run", value = "reads") %>%
      mutate(
        reads = as.integer(gsub(" ", "", reads)),
        step = "LSUx"
      ),
    enframe(filtercounts_full, name = "seq_run", value = "reads") %>%
      mutate(
        reads = as.integer(gsub(" ", "", reads)),
        step = "Filter (full)"
      ),
    enframe(filtercounts_ITS2, name = "seq_run", value = "reads") %>%
      mutate(
        reads = as.integer(gsub(" ", "", reads)),
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
    mutate(
      step = ifelse(tech == "Ion Torrent" & step == "Trim", "Raw", step),
    ) %>%
    select(tech, amplicon, step, reads, ASVs) %>%
    pivot_longer(c("reads", "ASVs"), names_to = "type", values_to = "count") %>%
    filter(!is.na(count)) %>%
    mutate(
      step = factor(step, levels = c("Raw", "Trim", "LSUx", "Filter (full)", "Filter (ITS2)", "ITS2",
                                     "short", "ITS", "LSU", "long")),
      tech = factor(tech, levels = c("PacBio", "Ion Torrent", "Illumina")),
      amplicon = factor(stringr::str_to_title(amplicon), levels = c("Long", "Short")),
      type = factor(type, c("ASVs", "reads"))
    ) %>%
    arrange(tech, amplicon, type) %>%
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
      proto_physeq %>%
      phyloseq::merge_samples(group = "seq_run") %>%
      phyloseq::otu_table() %>%
      t() %>%
      {.@.Data} %>%
      as_tibble(rownames = "label") %>%
      mutate_if(is.numeric, ~./sum(.)) %>%
      pivot_longer(-1, names_to = "seq_run", values_to = "reads") %>%
      filter(reads > 0) %>%
      left_join(select(datasets, seq_run, amplicon, tech), by = "seq_run") %>%
      mutate_at("amplicon", stringr::str_to_title) %>%
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
      mutate(
        kingdom = ifelse(!is.na(phylum) & reference == "warcup", "Fungi", kingdom)
      )
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
      conf_taxon_short_euk$tip_taxa %>%
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
      select(-Taxonomy)
  },
  
  tax_chart =
    taxon_reads %>%
    pivot_longer(
      names_to = "rank",
      cols = c("kingdom", "phylum", "class", "order", "family", "genus"),
      values_to = "taxon"
    ) %>%
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
    pivot_longer(cols = c("Reads", "ASVs"), names_to = "type", values_to = "frac") %>%
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
  
  nonpos_control_physeq = 
    phyloseq::subset_samples(
      proto_physeq,
      sample_type != "Pos" &
        !(seq_run == "is_057" & well == "E7")
    ),
  
  pos_control =
    agaricus_reads %>%
    as.matrix() %>%
    colSums() %>%
    which.max() %>%
    names(),
  
  pos_control_reads =
    nonpos_control_physeq %>%
    phyloseq::prune_taxa(taxa = pos_control) %>%
    phyloseq::otu_table() %>%
    rowSums(),
  
  all_reads =
    nonpos_control_physeq %>%
    phyloseq::otu_table() %>%
    as.matrix() %>%
    rowSums(),
  
  pos_control_data =
    nonpos_control_physeq %>%
    phyloseq::sample_data() %>%
    as.data.frame() %>%
    mutate(
      pc_reads = pos_control_reads,
      all_reads = all_reads,
      pc_frac = pc_reads/all_reads
    ) %>%
    select(seq_run, well, plate, pc_frac, pc_reads, all_reads),
  
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
  
  variog_data = target(
    ignore(plan2) %>%
    filter(step == "correlog") %>%
    select("correlog", timelag, guild, metric, tech, amplicon, algorithm) %>%
    mutate_at("correlog", lapply, readd, character_only = TRUE, cache = ignore(cache)) %>%
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
      cache = ignore(cache)
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
    ),
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
    mutate_at("variofit", map,
              ~ mutate(.y, gamma = predict(.x, newdata = .y)),
              expand_grid(dist = seq(1, 25, 0.5), timelag = 0:1)
    ) %>%
    unnest("variofit") %>%
    filter(mantel_1 | timelag == "0", mantel_0) %>%
    mutate_at("timelag", factor),
  
  taxdata = {
    ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
    out <- taxon_reads %>%
      mutate(domain = "Root", ASVs = 1) %>%
      select(-region, -seq_run, -label) %>%
      mutate_at(
        ranks,
        replace_na,
        "?"
      ) %>%
      group_by_at(vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
      mutate(ASVs = ASVs / sum(ASVs))
    for (i in seq_along(ranks)) {
      out <- out %>%
        group_by_at(rev(ranks)[i]) %>%
        group_map(
          function(x, y, i) {
            d <- group_by_at(x, vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
              summarize_at(c("reads", "ASVs"), sum)
            if (max(c((d$reads), (d$ASVs))) < 0.01) {
              mutate_at(x, rev(ranks)[i], ~"*")
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
    unknown_taxa <- out$taxon_names() %in% c("?", "*")
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
      filter(!is.na(as.character(ECM)), ECM != "non-ECM") %>%
      mutate(ASVs = 1) %>%
      select(-region, -seq_run, -label) %>%
      mutate_at(
        ranks,
        replace_na,
        "?"
      ) %>%
      group_by_at(vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
      mutate(ASVs = ASVs / sum(ASVs))
    for (i in seq_along(ranks)) {
      out <- out %>%
        group_by_at(rev(ranks)[i]) %>%
        group_map(
          function(x, y, i) {
            d <- group_by_at(x, vars(-one_of(ranks, "reads", "ASVs", "ECM"))) %>%
              summarize_at(c("reads", "ASVs"), sum)
            if (max(c((d$reads), (d$ASVs))) < 0.01) {
              mutate_at(x, rev(ranks)[i], ~"*")
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
    mutate(group = paste(tech, "-", amplicon) %>% factor(levels = unique(.))),
  
  full_dens_short = length_table %>%
    filter(amplicon == "Short", region == "short") %>%
    ggplot(aes(length, weight = reads, color = tech, group = tech)) +
    stat_density(bw = 0.5, geom = "line", position = "identity") +
    scale_color_brewer(type = "qual", guide = "none", palette = "Set2") +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
    ggnomics::facet_nested(rows = vars(tech), scales = "free_y") +
    theme(strip.background = element_blank()),
  
  full_dens_long = length_table %>%
    filter(amplicon == "Long", region == "long") %>%
    ggplot(aes(length, weight = reads, color = tech, group = tech)) +
    stat_density(bw = 0.5, geom = "line", position = "identity") +
    scale_color_brewer(type = "qual", guide = "none", palette = "Set2") +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
    ggnomics::facet_nested(rows = vars(tech), scales = "free_y") +
    theme(strip.background = element_blank()),
  
  full_ecdf_short = length_table %>%
    filter(amplicon == "Short", region == "short") %>%
    ggplot(aes(x = length, y = ecdf, color = tech, group = tech)) +
    geom_step() +
    scale_color_brewer(type = "qual",
                       name = NULL, palette = "Set2") +
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
    ggplot(aes(x = length, y = ecdf, color = tech, group = tech)) +
    geom_step() +
    scale_color_brewer(type = "qual",
                       name = NULL, palette = "Set2", direction = -1) +
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
    scale_color_brewer(type = "qual", guide = "none", palette = "Set2", direction = -1) +
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
    scale_color_brewer(type = "qual",
                       name = NULL, palette = "Set2", direction = -1) +
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
  
  ecm_heattree = {
    taxdata2 <- taxdata_ECM
    
    val_cols <- names(taxdata2$data$tax_read) %>%
      keep(str_detect, "reads_.+_.+_Unite_.+")
    
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
        str_replace("reads_.+_(.+)_Unite_.+", "\\1"),
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
    taxdata2 %>%
      metacoder::heat_tree(
        layout = "davidson-harel",
        initial_layout = "reingold-tilford",
        
        node_size = rowMeans(.$data$tax_read[,-1]),
        node_size_range = c(0.002, .035),
        node_size_axis_label = "Read abundance",
        
        node_color = read_ratio,
        node_color_range = metacoder::diverging_palette(),
        node_color_trans = "linear",
        node_color_axis_label = "Read abundance ratio",
        node_color_interval = c(-3, 3),
        
        edge_size = rowMeans(.$data$tax_asv[,-1]),
        edge_size_range = c(0.001, 0.03),
        edge_size_trans = "linear",
        edge_size_axis_label = "ASV count",
        
        edge_color = ASV_ratio,
        edge_color_range = metacoder::diverging_palette(),
        edge_color_trans = "linear",
        edge_color_axis_label = "ASV count ratio",
        edge_color_interval = c(-0.75, 0.75),
        
        node_label = taxon_names,
        node_label_size_range = c(0.02, 0.03),
        # title = unique(paste(.$data$diff_read$treatment_1, "vs.", .$data$diff_read$treatment_2)),
        aspect_ratio = 3/2
      )
  },
  trace = TRUE
) %>%
  dplyr::filter(ifelse(amplicon == '"Short"', metric != '"wunifrac"', TRUE) %|% TRUE)

saveRDS(plan2, "data/plan/drake_light.rds")
options(clustermq.scheduler = "multicore")
if (interactive()) {
  cache <- drake_cache(".light")
  # dconfig <- drake_config(plan2, cache = cache)
  # vis_drake_graph(dconfig)
  make(plan2, cache = cache)#, parallelism = "clustermq", jobs = local_cpus() - 1)
}
