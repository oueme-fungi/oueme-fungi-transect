if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}
library(drake)
library(magrittr)

source(file.path(config$rdir, "variogram.R"))

guild_metric_meta <- tidyr::crossing(
  metric = c("bray", "wunifrac"),
  guild = c("fungi", "ecm")
)

physeq_meta <-
  tidyr::crossing(
    dplyr::select(datasets, "seq_run", "tech", "dataset"),
    guild_metric_meta
  ) %>%
  dplyr::mutate(
    gm = glue::glue("{guild}_{metric}"),
    amplicon = stringr::str_extract(dataset, "^[a-z]+"),
    physeq_all = glue::glue("physeq_all_{gm}")
  ) %>%
  dplyr::filter(
    seq_run != "is_057", # couldn't be demultiplexed
    !(metric == "wunifrac" & amplicon == "short") # no good tree
  ) %>%
  dplyr::mutate_at("physeq_all", syms) %>%
  dplyr::select(-dataset)

correlog_meta <- tidyr::crossing(
  physeq_meta,
  timelag = c("0", "1")
) %>%
  dplyr::mutate(
    gmta = glue::glue("{gm}_{tech}_{amplicon}"),
    dist = glue::glue("dist_{gmta}"),
    dist_spatial = glue::glue("dist_spatial_{timelag}_{gmta}")
  ) %>%
  dplyr::mutate_at(c("dist", "dist_spatial"), make.names) %>%
  dplyr::mutate_at(c("dist", "dist_spatial"), rlang::syms) %>%
  dplyr::select(-guild, -metric, -tech, -amplicon)

plan2 <- drake_plan(
  # targets which are imported from first plan
  taxon_table = target(trigger = trigger(mode = "blacklist")),
  taxon = target(
    transform = map(.data = !!taxonomy_meta, .tag_in = step, .id = tax_ID),
    trigger = trigger(mode = "blacklist")
  ),
  raxml_decipher_LSU = target(trigger = trigger(mode = "blacklist")),
  raxml_decipher_long = target(trigger = trigger(mode = "blacklist")),
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
  
  # big_fasta ----
  # write the big_seq_table as a fasta files so that they can be clustered by
  # VSEARCH.
  big_fasta = target(
    write_big_fasta(big_seq_table,
                    file_out(big_fasta_file)),
    transform = map(.data = !!region_meta, .id = region)
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
    transform = map(outname = "decipher_LSU", group = "euk", .tag_out = c(euktree, tree))
  ),
  
  tree_decipher_LSU_long = target(
    raxml_decipher_long$bipartitions,
    transform = map(outname = "decipher_LSU_long", group = "euk", .tag_out = c(euktree, tree))
  ),
  
  # tree_epa_mafft_full = target(
  #   raxml_epa_mafft_full$bestTree,
  #   transform = map(outname = "epa_mafft_full", .tag_out = tree)
  # ),
  
  labeled_tree = target(
    relabel_tree(
      tree = tree,
      old = taxon_labels$label,
      new = taxon_labels$tip_label,
      chimeras = allchimeras_ITS2
    ) %T>%
      castor::write_tree(file_out(!!glue::glue(
        "data/trees/{outname}_{group}.tree",
        outname = outname,
        group = group
      ))),
    transform = map(tree, outname, group, .id = c(outname, group))
  ),
  
  phylo_labeled_tree = target(
    relabel_tree(
      tree = tree,
      old = phylotaxon_labels$label,
      new = phylotaxon_labels$tip_label,
      chimeras = allchimeras_ITS2
    ) %T>%
      castor::write_tree(file_out(!!glue::glue(
        "data/trees/{outname}_{group}_phylo.tree",
        outname = outname,
        group = group
      ))),
    transform = map(tree, phylotaxon_labels, outname, group, .id = c(outname, group))
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
    ape::root(tree, best_nonfungus),
    transform = map(euktree, .id = outname)
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
    transform = map(rooted_tree, group = "fungi", .tag_out = tree, .id = outname)
  ),
  
  # labeled_fungi_tree = target(
  #   relabel_tree(
  #     tree = fungi_tree,
  #     old = taxon_labels$label,
  #     new = taxon_labels$tip_label
  #   ) %T>%
  #     castor::write_tree(
  #       file_out(!!paste0("data/trees/fungi_", outname, ".tree"))
  #     ),
  #   transform = map(fungi_tree, outname, .id = outname)
  # ),
  # 
  # phylolabeled_fungi_tree = target(
  #   relabel_tree(
  #     tree = fungi_tree,
  #     old = phylotaxon_labels$label,
  #     new = phylotaxon_labels$tip_label
  #   ) %T>%
  #     castor::write_tree(
  #       file_out(!!paste0("data/trees/fungi_", outname, "_phylo.tree"))
  #     ),
  #   transform = map(fungi_tree, phylotaxon_labels, outname, .id = outname)
  # ),
  # 
  phylotaxon = target(
    phylotax(tree = tree, taxa = taxon_table),
    transform = map(tree, .id = c(outname, group))
  ),
  
  phylotaxon_labels = target(
    make_taxon_labels(phylotaxon$tip_taxa),
    transform = map(phylotaxon, .id = c(outname, group))
  ),
  
  guilds = target(
    phylotaxon$tip_taxa %>%
      dplyr::group_by(label, rank) %>%
      dplyr::filter((!"ITS" %in% region) | region != "short") %>%
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
    transform = map(phylotaxon, .id = c(outname, group))
  ),
  
  ecm = target(
      dplyr::filter(guilds, grepl("Ectomycorrhizal", guild)),
      transform = map(guilds, .id = c(outname, group))
  ),
  
  physeq = target(
    assemble_physeq(
      platemap = platemap,
      datasets = datasets,
      seqtable = relabel_seqtable(big_seq_table_ITS2),
      tree = if (metric == "wunifrac") fungi_tree_decipher_LSU_long else NULL,
      chimeras = allchimeras_ITS2
    ) %>%
      phyloseq::prune_samples(
        samples = phyloseq::sample_data(.)$sample_type == "Sample" &
          (phyloseq::sample_data(.)$qual == "X" |
             phyloseq::sample_data(.)$year == "2015") &
          rowSums(phyloseq::otu_table(.)) > 0 &
          phyloseq::sample_data(.)[["tech"]] == tech &
          phyloseq::sample_data(.)[["amplicon"]] == amplicon
      ) %>%
      phyloseq::prune_taxa(
        taxa = if (guild == "ecm") ecm_decipher_LSU_long_fungi$label
        else phyloseq::taxa_names(.)
      ) %>%
      phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0),
    transform = map(.data = !!physeq_meta, .id = c(guild, metric, tech, amplicon))
  ),
  
  dist = target(
    phyloseq::distance(physeq, metric),
    transform = map(physeq, metric, .id = c(gm, tech, amplicon))
  ),
  
  dist_spatial_0 = target(
    phyloseq::sample_data(physeq) %>%
    with(x + 30000 * as.integer(site) + 100000 * as.integer(year)) %>%
    dist(),
    transform = map(
      physeq,
      .id = c(gm, tech, amplicon),
      .tag_out = "dist_spatial"
      )
  ),
  
  dist_spatial_1 = target(
    dist_spatial_0 + ifelse(dist_spatial_0 > 50000, -100000, 100000),
    transform = map(
      dist_spatial_0,
      .id = c(gm, tech, amplicon),
      .tag_out = "dist_spatial"
      )
  ),
  
  correlog = target(
    vegan::mantel.correlog(
      dist,
      dist_spatial,
      break.pts = 0:13 - 0.5,
      cutoff = FALSE
    ),
    transform = map(.data = !!correlog_meta, .id = c(gmta, timelag))
  ),
  
  variog = target(
    variogram_dist(
      eco_dist = dist,
      sp_dist = dist_spatial,
      breaks = c(1:21 - 0.5, 31000)
    ),
    transform = map(.data = !!dplyr::filter(correlog_meta, timelag == 0),
                    .id = c(gmta))
  ),
  
  variofit = target(
    gstat::fit.variogram(
      variog,
      gstat::vgm(variog$gamma[22], "Exp", 3, variog$gamma[22]/2),
      fit.method = 1
    ),
    transform = map(variog, .id = c(gmta))
  ),
  trace = TRUE
)

saveRDS(plan2, "data/plan/drake_light.rds")

if (interactive()) {
  cache <- drake_cache(".light")
  dconfig <- drake_config(plan2, cache = cache)
  vis_drake_graph(dconfig)
  make(plan2, cache = cache)
}
# 