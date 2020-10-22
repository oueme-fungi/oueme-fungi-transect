# functions which draw a figure or generate data used to draw a figure
# author Brendan Furneaux

# draw the rDNA, primers, and amplicons (Fig1a)
draw_amplicon_plot <- function(amplicons_file) {
  readr::read_csv(amplicons_file) %>%
    dplyr::mutate(xmid = (x + xend) / 2,
                  label = tidyr::replace_na(label, "")) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend, color = color,
               linetype = linetype, size = size)) +
    geom_segment(data = ~dplyr::filter(., !arrow)) +
    geom_segment(
      data = ~dplyr::filter(., arrow),
      arrow = arrow(length = unit(0.1, "inches"))
    ) +
    geom_text(aes(x = xmid, y = y_label, label = label, color = labelcolor,
                  size = size_label)) +
    scale_linetype_identity() +
    scale_color_identity() +
    scale_size_identity() +
    ylim(c(-19, 24)) +
    coord_equal(ratio = 10) +
    theme(axis.line = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()
    )
}

# list of all sequences by region and sequencing run, along with number of
# reads and sequence length
# [tibble] columns: seq_run, region, seq, reads, length
compile_reads_table <- function(...) {
  purrr::map(list(...), tibble::as_tibble, rownames = "sample") %>%
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
    dplyr::ungroup()
}

# plot of ASV richness vs. denoising success vs. length
draw_length_richness_plot <- function(region_table, readcounts) {
  region_table %>%
    dplyr::filter(seq_run == "pb_500") %>%
    mutate_at(
      "reads",
      divide_by,
      sum(filter(readcounts, seq_run == "pb_500", step == "demux")$nreads)
    ) %>%
    dplyr::mutate_at("region", fct_recode, "5.8S" = "5_8S") %>%
    ggplot(aes(reads, ASVs, color = length_med, label = region)) +
    geom_point() +
    ggrepel::geom_text_repel() +
    scale_color_gradientn(
      colors = c("blue", "red", "green1"),
      limits = c(30, 2000),
      breaks = c(50, 150, 400, 1000),
      trans = "log10",
      name = "median length"
    ) +
    scale_x_continuous(name = "Fraction of reads mapped",limits = c(0, 1)) +
    ylab("ASV richness")
}

# Table of reads per ASV for each sequencing run
# [tibble] columns: seq, is_057, pb_483, pb_500, SH-2257
compile_asv_table <- function(seq_table) {
  tibble::as_tibble(seq_table, rownames = "filename") %>%
  tidyr::gather(key = "seq", value = "reads", -1) %>%
  dplyr::filter(reads >= 1) %>%
  tidyr::extract(col = "filename",
                 into = c("seq.run", "plate", "well", "dir", "region"),
                 regex = "([a-zA-Z]+[-_]\\d+)_(\\d+)_([A-H]1?[0-9])([fr]?)_([:alnum:]+).+") %>%
  dplyr::group_by(seq.run, seq) %>%
  dplyr::summarize(reads = sum(reads)) %>%
  tidyr::spread(key = seq.run, value = reads, fill = 0) %>%
  dplyr::ungroup()
}

# Table of reads per OTU for each sequencing run
# [tibble] columns: seq, is_057, pb_483, pb_500, SH-2257
compile_otu_table <- function(infile, ...) {
  read_tsv(
    infile,
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
}


# Find the central sequence for each 97% OTU
# [tibble] identity, ASVseq, OTUseq
# ASVseq and OTUseq are ASV sequence hashes.
generate_otu_map <- function(infile) {
  read_tsv(
    infile,
    col_names = paste0("V", 1:10),
    col_types = "ciidfccccc",
    na = c("", "NA", "*")
  ) %>%
  #filter(V1 == "H") %>%
  select(identity = V4, ASVseq = V9, OTUseq = V10) %>%
  mutate_all(str_replace, ";.*", "") %>%
  mutate(OTUseq = coalesce(OTUseq, ASVseq),
         identity = tidyr::replace_na(identity, 100)) %>%
  unique()
}

count_asvs <- function(asv_table, allseqs) {
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
    select(seq_run, reads, step = region, ASVs)
}

draw_phylotax_figure <- function() {
  treeio::as.treedata(
    ape::read.tree(text = "(A:1,((B:1,C:1):1,((E:1,F:1):1,D:1):1):1);")
  ) %>%
    tidytree::full_join(
      tibble::tibble(label = LETTERS[1:6],
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
    scale_y_reverse()
}


compile_length_table <- function(asv_table, allseqs, datasets) {
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
    dplyr::left_join(
      dplyr::select(datasets, seq_run, tech, amplicon),
      by = "seq_run"
    ) %>%
    dplyr::mutate_at("amplicon", factor, levels = c("Short", "Long")) %>%
    dplyr::mutate_at("tech", factor, levels = c("PacBio", "Illumina", "Ion Torrent")) %>%
    arrange(tech, amplicon) %>%
    mutate(group = paste(tech, amplicon) %>% factor(levels = unique(.)))
}

length_density_plot <- function(length_table, nesting, ...) {
    filter(length_table, ...) %>%
    ggplot(aes(length, weight = reads, color = group, group = group)) +
    stat_density(bw = 0.5, geom = "line", position = "identity") +
    scale_color_strategy(guide = "none") +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
    ggnomics::facet_nested(nesting, scales = "free_y", nest_line = TRUE) +
    theme(strip.background = element_blank())
}

full_length_ecdf_plot <- function(length_table, amplicon, xlabel = NULL) {
    filter(length_table,
           .data$amplicon == !!amplicon, region == tolower(!!amplicon)) %>%
    ggplot(aes(x = length, y = ecdf, color = group, group = group)) +
    geom_step() +
    scale_color_strategy(name = NULL) +
    scale_y_continuous(name = "ECDF") +
    scale_x_continuous(name = xlabel) +
    facet_grid( amplicon ~ ., scales = "free_x") +
    theme(
      legend.position = c(1, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_blank(),
      strip.background = element_blank()
    )
}

its2_length_ecdf_plot <- function(length_table) {
    filter(length_table, region == "ITS2") %>%
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
    )
}

draw_region_length_figure <- function(reads_table, regions) {
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
}

make_buffer_compare_physeq <- function(proto_physeq, taxon_reads) {
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
}

generate_buffer_asv_physeq <- function(proto_physeq, fungi) {
    # No Ion Torrent samples, because they are not completely demultiplexed
    phyloseq::subset_samples(proto_physeq, tech != "Ion Torrent") %>%
    # Don't use QC samples
    phyloseq::subset_samples(sample_type == "Sample") %>%
    # Only include fungi
    phyloseq::prune_taxa(taxa = fungi) %>%
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
    )

}

generate_tech_asv_physeq <- function(proto_physeq, fungi) {
    # Don't use Ion Torrent because it didn't multiplex properly
    phyloseq::subset_samples(proto_physeq, tech != "Ion Torrent") %>%
    # Only use the Xpedition buffer
    phyloseq::subset_samples(buffer == "Xpedition") %>%
    # Don't use QC sample
    phyloseq::subset_samples(sample_type == "Sample") %>%
    # Only fungi
    phyloseq::prune_taxa(taxa = fungi)
}

generate_tech_class_physeq <- function(tech_asv_physeq, fungi) {
  phyloseq::`tax_table<-`(
    tech_asv_physeq,
    phyloseq::tax_table(
      fungi %>%
        select(label, kingdom:genus) %>%
        unique() %>%
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
  )
}

draw_tech_class_pcoa_plot <- function(pcoa, sample_data, scores) {
  left_join(
    vegan::scores(pcoa, display = "sites") %>%
      as_tibble(rownames = "sample"),
    sample_data %>%
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
    data = scores,
    inherit.aes = FALSE
  ) +
  scale_color_discrete(name = NULL) +
  scale_discrete_manual(aesthetics = "pch", name = NULL, values = 1:3) +
  coord_equal() +
  xlab("PCoA1") +
  ylab("PCoA2")
}

# Add taxonomy to the ECM ASV table and cluster at the family level
generate_tech_ecm_fam_physeq <- function(tech_asv_physeq, ecm_combined)
  phyloseq::`tax_table<-`(
    tech_asv_physeq,
    phyloseq::tax_table(
      ecm_combined %>%
        select(label, kingdom:genus) %>%
        unique() %>%
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
  )

tree_figure <- function(phylotax, outfile, charscale = 0.015, rankgap = 0.015) {
  # tip_labels <-
  #   dplyr::filter(phylotax$retained, label %in% phylotax$tree$tip.label) %>%
  #   phylotax::make_taxon_labels(abbrev = TRUE)
  tip_taxa <-
    bind_rows(phylotax$assigned, phylotax$retained) %>%
    filter(rank != "kingdom", label %in% phylotax$tree$tip.label) %>%
    group_by(label, rank) %>%
    filter(if (any("PHYLOTAX" == method)) method == "PHYLOTAX" else TRUE) %>%
    group_by(label, rank) %>%
    summarize(
      taxon =  table(taxon) %>%
        paste0(names(.), collapse = "/") %>%
        gsub(pattern = "(.+/.+)", replacement = "<\\1>") %>%
        gsub(pattern = "^\\d+", replacement = "")
    ) %>%
    mutate(
      tip = match(label, phylotax$tree$tip.label),
      depth = ape::node.depth.edgelength(phylotax$tree)[tip]
    ) %>%
    group_by(label) %>%
    arrange(desc(rank), .by_group = TRUE) %>%
    mutate(width = nchar(taxon) * charscale + rankgap,
           offset = cumsum(lag(width, default = 0))) %>%
    ungroup()
  clade_annot <- list()
  ranks <- c("genus", "family", "order", "class", "phylum")
  for (i in seq_along(ranks)) {
    rank_tips <- select(tip_taxa, tip, rank, taxon) %>%
      filter(rank == ranks[[i]])
    rank_nodes <- phylotax$node_assigned %>%
      filter(rank == ranks[i]) %>%
      select(node, rank, taxon) %>%
      unique() %>%
      mutate(tip = phangorn::Descendants(
        phylotax$tree,
        node = node,
        type = "tips"
      )) %>%
      unnest(tip) %>%
      bind_rows(
        anti_join(filter(rank_tips, !startsWith(taxon, "<")), ., by = "tip")
      ) %>%
      mutate(node = coalesce(node, tip)) %>%
      group_by(taxon) %>%
      mutate(color = ifelse(n_distinct(node) > 1, "tomato", "black")) %>%
      ungroup() %>%
      left_join(tip_taxa, by = c("rank", "tip", "taxon")) %>%
      group_by(node, taxon, color) %>%
      summarize(offset = max(depth + offset) - max(depth)) %>%
      rename(label = taxon) %>%
      bind_rows(
        tip_taxa %>%
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
  (phylotax$tree %>%
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
      filename = outfile,
      plot = .,
      device = "pdf",
      width = 8,
      height = 20,
      units = "in"
    )
}

