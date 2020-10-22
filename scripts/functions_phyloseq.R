# functions to build phyloseq objects from subsets of other data
# author Brendan Furneaux

# rename the rows and columns of a sequence abundance table
# convert rownames from a fileame to {Seq.Run}{Plate}{Well}
# convert colnames from the full sequence to an 8-character hash.
# also pools forward and reverse reads from same well.

relabel_seqtable <- function(seqtable) {
  seqtable %>%
    # set the new column names
    magrittr::set_colnames(tzara::seqhash(chartr("T", "U", colnames(.)))) %>%
    # convert to a tibble for easier name column operations
    tibble::as_tibble(rownames = "file") %>%
    # parse the filename
    tidyr::extract(
      col = "file",
      into = c("seq_run", "plate", "well", "direction", "region"),
      regex = "([a-zA-Z]+[_-]\\d+)_(\\d+)_([A-H]1?\\d)([fr]?)_([a-zA-Z0-9]+).*") %>%
    # create the ID
    tidyr::unite("ID", seq_run, plate, well, sep = "") %>%
    # remove unnecessary columns
    dplyr::select(-direction, -region) %>%
    # pool all reads from the same well
    dplyr::group_by(ID) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::ungroup() %>%
    # convert back to a matrix
    tibble::column_to_rownames("ID") %>%
    as.matrix()
}

assemble_sample_data <- function(platemap, datasets) {
  platemap %>%
    dplyr::mutate_at("primer_pair", tolower) %>%
    dplyr::left_join(
      datasets %>%
        dplyr::select(dataset, seq_run, tech, forward, reverse, amplicon) %>%
        dplyr::mutate(
          primer_pair = paste(
            stringr::str_replace(forward, "_tag.*$", ""),
            stringr::str_replace(reverse, "_tag.*$", ""),
            sep = "_")
        ),
      by = "primer_pair") %>%
    dplyr::mutate_at(c("year", "site", "qual", "plate", "well", "primer_pair"),
                     factor) %>%
    tidyr::unite("ID", seq_run, plate, well, sep = "", remove = FALSE) %>%
    tibble::column_to_rownames("ID") %>%
    phyloseq::sample_data()
}

assemble_otu_table <- function(samp, seqtable, drop) {
  asvs <- colnames(seqtable)
  asvs <- setdiff(asvs, drop)
  miss_samples <- setdiff(row.names(samp), row.names(seqtable))
  rbind(
    seqtable[,asvs, drop = FALSE],
    matrix(0, nrow = length(miss_samples), ncol = length(asvs),
           dimnames = list(miss_samples, asvs))
  ) %>%
    phyloseq::otu_table(taxa_are_rows = FALSE)
}

build_physeq <- function(taxa, physeq_otu_table, sample_data, amplicon, tech,
                         guild, fungi, ecm, ecm2, ecm3) {
  {
    physeq <-
      phylotax::phylotax_to_phyloseq(
        phylotax = taxa,
        otu_table = physeq_otu_table,
        sample_data,
        use_tree = !is.null(taxa$tree) && amplicon == "Long"
      ) %>%
      phyloseq::subset_samples(sample_type == "Sample") %>%
      phyloseq::subset_samples(buffer == "Xpedition") %>%
      phyloseq::subset_samples(site == "Gan") %>%
      phyloseq::prune_samples(
        samples = rowSums(phyloseq::otu_table(.)) > 100
      ) %>%
      phyloseq::prune_samples(
        samples = phyloseq::sample_data(.)[["tech"]] == tech
      ) %>%
      phyloseq::prune_samples(
        samples = phyloseq::sample_data(.)[["amplicon"]] == amplicon
      )

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
  }
}
