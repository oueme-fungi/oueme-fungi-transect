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

source(file.path(config$rdir, "variogram.R"))

physeq_meta <-
  tidyr::crossing(
    dplyr::select(datasets, "seq_run", "tech", "dataset", "amplicon"),
    guild = c("fungi", "ecm") #, "ecm2", "ecm3")
  ) %>%
  dplyr::filter(
    seq_run != "is_057", # couldn't be demultiplexed
  ) %>%
  dplyr::select(-dataset)

plan2 <- drake_plan(
  # targets which are imported from first plan
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
    transform = map(outname = "decipher_LSU", group = "euk", .tag_in = step, .tag_out = c(euktree, tree))
  ),
  
  tree_decipher_LSU_long = target(
    raxml_decipher_long$bipartitions,
    transform = map(outname = "decipher_LSU_long", group = "euk", .tag_in = step, .tag_out = c(euktree, tree))
  ),
  
  tree_decipher_unconst_long = target(
    raxml_decipher_unconst_long$bipartitions,
    transform = map(outname = "decipher_unconst_long", group = "euk", .tag_in = step, .tag_out = c(euktree, tree))
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
    transform = map(tree, outname, group, .tag_in = step, .id = c(outname, group))
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
    transform = map(tree, phylotaxon_labels, outname, group, .tag_in = step, .id = c(outname, group))
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
    transform = map(euktree, .tag_in = step, .id = outname)
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
    transform = map(rooted_tree, group = "fungi", .tag_in = step, .tag_out = tree, .id = outname)
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
    transform = map(tree, .tag_in = step, .id = c(outname, group))
  ),
  
  phylotaxon_labels = target(
    make_taxon_labels(phylotaxon$tip_taxa),
    transform = map(phylotaxon, .tag_in = step, .id = c(outname, group))
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
    transform = map(phylotaxon, .tag_in = step, .id = c(outname, group))
  ),
  
  ecm = target(
      dplyr::filter(guilds, grepl("Ectomycorrhizal", guild)),
      transform = map(guilds, .tag_in = step, .id = c(outname, group))
  ),
  
  ecm2 = target(
    dplyr::filter(ecm, !taxon %in% c("Peziza", "Pezizaceae", "Pyronemataceae")),
    transform = map(ecm, .tag_in = step, .id = c(outname, group))
  ),
  
  ecm3 = target(
    dplyr::filter(ecm, confidenceRanking != "Possible"),
    transform = map(ecm, .tag_in = step, .id = c(outname, group))
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
        phyloseq::subset_samples(qual == "X" | year == "2015") %>%
        phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0 ) %>%
        phyloseq::prune_samples(samples = phyloseq::sample_data(.)[["tech"]] == tech) %>%
        phyloseq::prune_samples(samples = phyloseq::sample_data(.)[["amplicon"]] == amplicon)
      
      if (guild == "ecm") {
        physeq <- phyloseq::prune_taxa(ecm_decipher_unconst_long_fungi$label, physeq) %>%
          phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0)
          
      } else if (guild == "ecm2") {
        physeq <- phyloseq::prune_taxa(ecm2_decipher_unconst_long_fungi$label, physeq) %>%
          phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0)
        
      } else if (guild == "ecm3") {
        physeq <- phyloseq::prune_taxa(ecm3_decipher_unconst_long_fungi$label, physeq) %>%
          phyloseq::prune_samples(samples = rowSums(phyloseq::otu_table(.)) > 0)
        
      }
      physeq
    },
    transform = map(.data = !!physeq_meta, .tag_in = step, .id = c(guild, tech, amplicon))
  ),
  
  correlog = target(
    correlog(physeq, metric, timelag),
    transform = cross(physeq, metric = c("bray", "wunifrac"), timelag = c(0, 1),
                    .tag_in = step, .id = c(guild, metric, tech, amplicon, timelag))
  ),
  
  variog = target(
    variog(physeq, metric, breaks = c(1:21 - 0.5, 31000)),
    transform = cross(physeq, metric = c("bray", "wunifrac"),
                    .tag_in = step,
                    .id = c(guild, metric, tech, amplicon))
  ),
  
  variofit = target(
    gstat::fit.variogram(
      as_variogram(variog),
      gstat::vgm(variog$gamma[22], "Exp", 3, variog$gamma[22]/2),
      fit.method = 1
    ),
    transform = map(variog, .tag_in = step,.id = c(guild, metric, tech, amplicon))
  ),
  
  variogST = target(
    variogST(physeq, metric, breaks = c(1:26 - 0.5, 30000)),
    transform = cross(physeq, metric = c("bray", "wunifrac"),
                      .tag_in = step,
                      .id = c(guild, metric, tech, amplicon))
  ),
  
  variofitST = target(
    gstat::fit.StVariogram(
      as_variogramST(variogST) %>%
        inset(,"np", ifelse(.$dist > 30, 1, .$np)),
      gstat::vgmST(
        "metric",
        joint = gstat::vgm(max(variogST$gamma), "Exp", 5, quantile(variogST$gamma, 0.25)),
        stAni = 1
      ),
      fit.method = 1
    ),
    transform = map(variogST, .tag_in = step, .id = c(guild, metric, tech, amplicon))
  ),
  trace = TRUE
) %>%
  dplyr::filter(ifelse(amplicon == '"Short"', metric != '"wunifrac"', TRUE) %|% TRUE)

saveRDS(plan2, "data/plan/drake_light.rds")

if (interactive()) {
  cache <- drake_cache(".light")
  dconfig <- drake_config(plan2, cache = cache)
  vis_drake_graph(dconfig)
  make(plan2, cache = cache)
}
