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
source(file.path(config$rdir, "qstats.R"))

choosevars <- function(d, g, .data) {
  .data %>%
    filter(type == g$type) %>%
    select(type, x = !!g$x_var, y = !!g$y_var) %>%
    filter(x > 0 | y > 0) %>%
    mutate(x_var = g$x_var, y_var = g$y_var)
}

datasets <- read_csv(config$dataset, col_types = "cccccccicccccicc") %>%
  mutate_at("amplicon", stringr::str_to_lower) %>%
  filter(tech != "Illumina")
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
      if (amplicon == "long") {
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
  
  reads_table = target(
    list(big_seq_table) %>%
      purrr::map(tibble::as_tibble, rownames = "sample") %>%
      purrr::map_dfr(tidyr::pivot_longer, -"sample", names_to = "seq", values_to = "reads") %>%
      dplyr::filter(reads > 0) %>%
      tidyr::extract(
        "sample",
        c("seq_run", "plate", "well", "region"),
        regex = "([pi][bs]_\\d{3})_(\\d{3})_([A-H]1?\\d)_(.+)"
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
                   regex = "([a-z]+_\\d+)_(\\d+)_([A-H]1?[0-9])([fr]?)_([:alnum:]+).+") %>%
    dplyr::group_by(seq.run, seq) %>%
    dplyr::summarize(reads = sum(reads)) %>%
    tidyr::spread(key = seq.run, value = reads, fill = 0),
  
  otu_table = read_tsv(
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
      regex = "([a-z]{2}_\\d{3})_(\\d{3})([A-H]1?\\d)"
    ) %>%
    tidyr::gather(key = "seq", value = "reads", -(1:3)) %>%
    dplyr::filter(reads >= 1) %>%
    dplyr::group_by(seq_run, seq) %>%
    dplyr::summarize(reads = sum(reads)) %>%
    tidyr::spread(key = seq_run, value = reads, fill = 0),
  
  demuxlength =  parse_qstat(qstats_length) %>%
    filter(is.na(step), !is.na(well)) %>%
    group_by(seq_run) %>%
    summarize(length = reldist::wtd.quantile(length, weight = nreads)) %>%
    deframe(),
  
  demuxqual = parse_qstat(qstats_erate) %>%
    filter(is.na(step), !is.na(well)) %>%
    group_by(seq_run) %>%
    summarize(qual = round(weighted.mean(erate, w = nreads, na.rm = TRUE), 4)) %>%
    deframe(),
  
  readcounts = parse_qstat(qstats_n),
  
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
    filter(region_table, region == "ITS2") %>%
      select(seq_run, reads, step = region, ASVs),
    pivot_longer(asv_table, -1, names_to = "seq_run", values_to = "reads") %>%
      filter(reads > 0) %>%
      mutate_at("seq", chartr, old = "T", new = "U") %>%
      inner_join(
        select(readd(allseqs, cache = cache), seq = ITS2, long, short, ITS, LSU) %>%
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
      step = factor(step, levels = c("Raw", "Trim", "LSUx", "ITS2",
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
  
  venn_ASV =
    asv_table %>%
    select(seq, pb_500, pb_483, is_057) %>%
    mutate_if(is.numeric, list(found = ~. > 0)) %>%
    group_by_at(vars(ends_with("_found"))) %>%
    summarize_if(is.numeric, list(ASVs = ~sum(. > 0), reads = sum)) %>%
    ungroup() %>%
    mutate_at(vars(ends_with("_found")), as.integer) %>%
    tidyr::unite(col = "set", 1:3, sep = "") %>%
    pivot_longer(
      -1,
      names_to = c("seq_run", "what"),
      values_to = "value",
      names_pattern = "([[:alpha:]]+[-._]\\d+)_(.+)"
    ) %>%
    pivot_wider(names_from = "what", values_from = "value") %>%
    group_by(seq_run) %>%
    mutate_if(is.numeric, list(frac = ~./sum(.))) %>%
    mutate_at(
      vars(ends_with("_frac")),
      formatC,
      digits = 2,
      format = "f",
      drop0trailing = FALSE
    ) %>%
    mutate_at(
      c("reads", "ASVs"),
      format
    ) %>%
    mutate(
      ASVs_frac = if_else(grepl(" +0", ASVs), "", ASVs_frac),
      ASVs = if_else(grepl(" +0", ASVs), "", ASVs),
      reads_frac = if_else(grepl(" +0", reads), "", reads_frac),
      reads = if_else(grepl(" +0", reads), "", reads)
    ) %>%
    select(set, seq_run, ASVs, ASVs_frac, reads, reads_frac) %>% {
      left_join(
        group_by(., set) %>%
          summarize(ASVs = as.integer(first(ASVs[ASVs != ""]))),
        
        select(., set, seq_run, ASV_frac = ASVs_frac, reads, frac = reads_frac) %>%
          pivot_wider(
            names_from = "seq_run",
            values_from = c("ASV_frac", "reads", "frac")
          ),
        by = "set"
      )
    } %>%
    mutate(n = str_count(set, "1")) %>%
    arrange(n, desc(set)) %>%
    select(set, ASVs, ends_with("pb_500"), ends_with("pb_483"), ends_with("is_057")) %>%
    column_to_rownames("set"),
  
  venn_OTU =
    otu_table %>%
    select(seq, pb_500, pb_483, is_057) %>%
    mutate_if(is.numeric, list(found = ~. > 0)) %>%
    group_by_at(vars(ends_with("_found"))) %>%
    summarize_if(is.numeric, list(OTUs = ~sum(. > 0), reads = sum)) %>%
    ungroup() %>%
    mutate_at(vars(ends_with("_found")), as.integer) %>%
    tidyr::unite(col = "set", 1:3, sep = "") %>%
    pivot_longer(
      -1,
      names_to = c("seq_run", "what"),
      values_to = "value",
      names_pattern = "([[:alpha:]]+[-._]\\d+)_(.+)"
    ) %>%
    pivot_wider(names_from = "what", values_from = "value") %>%
    group_by(seq_run) %>%
    mutate_if(is.numeric, list(frac = ~./sum(.))) %>%
    mutate_at(
      vars(ends_with("_frac")),
      formatC,
      digits = 2,
      format = "f",
      drop0trailing = FALSE
    ) %>%
    mutate_at(
      c("reads", "OTUs"),
      format
    ) %>%
    mutate(
      OTUs_frac = if_else(grepl(" +0", OTUs), "", OTUs_frac),
      OTUs = if_else(grepl(" +0", OTUs), "", OTUs),
      reads_frac = if_else(grepl(" +0", reads), "", reads_frac),
      reads = if_else(grepl(" +0", reads), "", reads)
    ) %>%
    select(set, seq_run, OTUs, OTUs_frac, reads, reads_frac) %>%  {
      left_join(
        group_by(., set) %>%
          summarize(OTUs = as.integer(first(OTUs[OTUs != ""]))),
        
        select(., set, seq_run, OTU_frac = OTUs_frac, reads, frac = reads_frac) %>%
          pivot_wider(
            names_from = "seq_run",
            values_from = c("OTU_frac", "reads", "frac")
          ),
        by = "set"
      )
    } %>%
    mutate(n = str_count(set, "1")) %>%
    arrange(n, desc(set)) %>%
    select(set, OTUs, ends_with("pb_500"), ends_with("pb_483"), ends_with("is_057")) %>%
    column_to_rownames("set"),
  
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
  
  multi_table = crossing(
    x_var = factor(glue::glue_data(datasets, "{tech} {machine} - {amplicon}"), ordered = TRUE),
    y_var = factor(glue::glue_data(datasets, "{tech} {machine} - {amplicon}"), ordered = TRUE),
    type = c("ASV", "OTU")
  ) %>%
    filter(
      x_var < y_var,
      x_var %in% names(big_table),
      y_var %in% names(big_table)
    ) %>%
    mutate_all(as.character) %>%
    group_by(x_var, y_var, type) %>%
    group_map(choosevars, .data = big_table) %>%
    bind_rows() %>%
    group_by(x_var, y_var, type),
  
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
        kingdom = ifelse(endsWith(phylum, "mycota"), "Fungi", kingdom)) %>%
      bind_rows(out) %>%
    mutate_at(
      "method",
      factor,
              levels = c("PHYLO", "dada2", "sintax", "idtaxa"),
              labels = c("PHYLOTAX", "RDPC", "SINTAX", "IDTAXA")
    ) %>%
    mutate_at(
      "reference",
      factor,
      levels = c("All", "unite", "warcup", "rdp_train"),
      labels = c("All", "Unite", "Warcup", "RDP")
    ) %>%
      mutate_at("tech", factor, levels = c("PacBio", "Ion Torrent")) %>%
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
                mutate_at(x, rev(ranks)[i], ~"...")
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
      ) %in% c("?", "...")
      unknown_taxa <- out$taxon_names() %in% c("?", "...")
      out <- out$filter_taxa(
        !(child_of_unknown & unknown_taxa) | taxon_names == "Root",
        reassign_obs = FALSE
      )
      dangling_taxa <- map_lgl(
        out$subtaxa(recursive = FALSE)[unlist(out$supertaxa(recursive = FALSE, na = TRUE))],
        ~all(out$taxon_names()[.] %in% c("?", "...") | is.na(out$taxon_names()[.]))
      )
      out <- out$filter_taxa(
        !dangling_taxa | taxon_names == "Root",
        reassign_obs = FALSE
      )
      
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
              mutate_at(x, rev(ranks)[i], ~"...")
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
    ) %in% c("?", "...")
    out <- out$filter_taxa(
      !child_of_unknown | taxon_names == "Fungi",
      reassign_obs = FALSE
    )
    dangling_taxa <- map_lgl(
      out$subtaxa(recursive = FALSE)[unlist(out$supertaxa(recursive = FALSE, na = TRUE))],
      ~all(out$taxon_names()[.] %in% c("?", "...") | is.na(out$taxon_names()[.]))
    )
    out <- out$filter_taxa(
      !dangling_taxa | taxon_names == "Fungi",
      reassign_obs = FALSE
    )
    
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
  trace = TRUE
) %>%
  dplyr::filter(ifelse(stringr::str_to_lower(amplicon) == '"short"', metric != '"wunifrac"', TRUE) %|% TRUE)

saveRDS(plan2, "data/plan/drake_light.rds")

if (interactive()) {
  cache <- drake_cache(".light")
  # dconfig <- drake_config(plan2, cache = cache)
  # vis_drake_graph(dconfig)
  make(plan2, cache = cache)
}
