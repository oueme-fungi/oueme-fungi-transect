if (exists("snakemake")) {
  snakemake@source(".Rprofile", echo = FALSE)
  load(snakemake@input[["drakedata"]])
} else {
  load("drake.Rdata")
}
library(drake)
library(magrittr)

correlog_meta <- tidyr::crossing(
  dplyr::select(datasets, "dataset"),
  metric = c("bray", "unifrac"),
  timelag = c("0", "1")
) %>%
  dplyr::mutate(
    dist = glue::glue("dist_{metric}_{dataset}"),
    dist_spatial = glue::glue("dist_spatial_{timelag}_{dataset}")
  ) %>%
  dplyr::mutate_at(c("dist", "dist_spatial"), make.names) %>%
  dplyr::mutate_at(c("dist", "dist_spatial"), rlang::syms) %>%
  dplyr::rename(dataset2 = dataset)

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
  
  labeled_tree_decipher_LSU =
    relabel_tree(
      tree = raxml_decipher_LSU$bipartitions,
      old = taxon_labels$label,
      new = taxon_labels$tip_label,
      chimeras = allchimeras_ITS2
    ) %T>%
    castor::write_tree(file_out("data/trees/decipher_LSU.tree")),
  
  labeled_tree_decipher_long =
    relabel_tree(
      tree = raxml_decipher_long$bipartitions,
      old = taxon_labels$label,
      new = taxon_labels$tip_label,
      chimeras = allchimeras_ITS2
    ) %T>%
    castor::write_tree(file_out("data/trees/decipher_long.tree")),
  
  labeled_tree_epa_full =
    relabel_tree(
      tree = raxml_epa_full$bestTree,
      old = taxon_labels$label,
      new = taxon_labels$tip_label,
      chimeras = allchimeras_ITS2
    ) %T>%
    castor::write_tree(file_out("data/trees/epa_full.tree")),
  
  phylo_labeled_tree_epa_full =
    relabel_tree(
      tree = raxml_epa_full$bestTree,
      old = phylotaxon_labels$label,
      new = phylotaxon_labels$tip_label,
      chimeras = allchimeras_ITS2
    ) %T>%
    castor::write_tree(file_out("data/trees/epa_full_phylo.tree")),
  
  
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
  
  # find the most abundant sequence which we are confident is a plant
  best_nonfungus =
    taxon_table %>%
    dplyr::ungroup() %>%
    dplyr::filter(rank == "kingdom",
                  taxon != "Fungi",
                  n_tot == 6,
                  n_diff == 1) %>%
    dplyr::filter(n_reads == max(n_reads)) %$%
    label[1],
  
  # root the tree outside the fungi.  This isn't accurate, but it ensures that
  # the fungi can be indentified.
  tree_epa_full =
    ape::root(raxml_epa_full$bestTree,
              best_nonfungus),
  
  # Extract the minimally inclusive clade including all confidently identified
  # fungi (except Tulasnella).
  fungi_tree_epa_full =
    ape::getMRCA(
      phy = tree_epa_full,
      tip = intersect(surefungi, tree_epa_full$tip.label)
    ) %>%
    ape::extract.clade(phy = tree_epa_full) %>%
    ape::drop.tip(phy = ., intersect(.$tip.label, allchimeras_ITS2)),
  
  labeled_fungi_tree_epa_full =
    relabel_tree(
      tree = fungi_tree_epa_full,
      old = taxon_labels$label,
      new = taxon_labels$tip_label
    ) %T>%
    castor::write_tree(file_out("data/trees/fungi_epa_full.tree")),
  
  phylolabeled_fungi_tree_epa_full =
    relabel_tree(
      tree = fungi_tree_epa_full,
      old = phylotaxon_labels$label,
      new = phylotaxon_labels$tip_label
    ) %T>%
    castor::write_tree(file_out("data/trees/fungi_epa_full_phylo.tree")),
  
  taxon_phy = phylotax(
    tree = fungi_tree_epa_full,
    taxa = taxon_table
  ),
  
  taxon_phy_long = phylotax(
    tree = ape::keep.tip(
      fungi_tree_epa_full,
      intersection(fungi_tree_epa_full)),
    taxa = taxon_table
  ),
  
  phylotaxon_labels = make_taxon_labels(taxon_phy$tip_taxa),
  
  physeq_all = assemble_physeq(
    platemap,
    datasets,
    relabel_seqtable(big_seq_table_ITS2),
    fungi_tree_epa_full,
    chimeras = allchimeras_ITS2
  ),
  
  physeq = target(
    physeq_all %>%
    phyloseq::prune_samples(
      samples = phyloseq::sample_data(.)[["dataset"]] == dataset & 
        phyloseq::sample_data(.)$sample_type == "Sample" &
        (phyloseq::sample_data(.)$qual == "X" |
           phyloseq::sample_data(.)$year == "2015") &
        rowSums(phyloseq::otu_table(.)) > 0
    ),
    transform = map(.data = !!datasets, .id = dataset)
  ),
  
  dist = target(
    phyloseq::distance(physeq, metric),
    transform = cross(
      physeq,
      metric = c("unifrac", "bray"),
      .id = c(metric, dataset)
    )
  ),
  
  dist_spatial_0 = target(
    phyloseq::sample_data(physeq) %>%
    with(x + 30000 * as.integer(site) + 100000 * as.integer(year)) %>%
    dist(),
    transform = map(physeq, .id = dataset, .tag_out = "dist_spatial")
  ),
  
  dist_spatial_1 = target(
    dist_spatial_0 + ifelse(dist_spatial_0 > 50000, -100000, 100000),
    transform = map(dist_spatial_0, .id = dataset, .tag_out = "dist_spatial")
  ),
  
  correlog = target(
    vegan::mantel.correlog(
      dist,
      dist_spatial,
      break.pts = 0:13 - 0.5,
      cutoff = FALSE
    ),
    transform = map(.data = !!correlog_meta, .id = c(metric, timelag, dataset2))
  ),
  trace = TRUE
)

if (interactive()) {
  cache <- drake_cache(".light")
  dconfig <- drake_config(plan2, cache = cache)
  vis_drake_graph(dconfig)
  make(plan2, cache = cache)
}
# 