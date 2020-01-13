
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
      regex = "([a-z]+_\\d+)_(\\d+)_([A-H]1?\\d)([fr]?)_([a-zA-Z0-9]+).*") %>%
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
    as.matrix
}

assemble_physeq <- function(platemap, datasets, seqtable, tree = NULL, chimeras) {
  samp <- platemap %>%
    dplyr::mutate_at("primer_pair", tolower) %>%
    dplyr::left_join(
      datasets %>%
        dplyr::select(dataset, seq_run, tech, forward, reverse) %>%
        dplyr::mutate(
          amplicon = stringr::str_extract(dataset, "^[a-z]+"),
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
  asvs <- colnames(seqtable)
  if (!is.null(tree)) asvs <- intersect(asvs, tree$tip.label)
  asvs <- setdiff(asvs, chimeras)
  asvtab <- seqtable[,asvs, drop = FALSE] %>%
    phyloseq::otu_table(taxa_are_rows = FALSE)
  
  if (is.null(tree)) {
    phyloseq::phyloseq(samp, asvtab)
  } else {
    phyloseq::phyloseq(samp, asvtab, ape::keep.tip(tree, asvs))
  }
}


max_cophenetic <- function(tree) {
  max_coph <- max_depth <- max_length <- numeric(length(tree$tip.label) + tree$Nnode)
  # tips have no depth, length, or cophenetic distance
  max_coph[seq_along(tree$tip.label)] <- 0
  max_depth[seq_along(tree$tip.label)] <- 0
  max_length[seq_along(tree$tip.label)] <- 0
  
  for (i in seq.int(length(max_coph), length(tree$tip.label) + 1, -1)) {
    edges <- which(tree$edge[,1] == i)
    if (length(edges) == 0) {
      max_coph[i] <- max_length[i] <- max_depth[i] <- 0
    } else {
      max_length[i] <- max(max_length[tree$edge[edges,2]] + tree$edge.length[edges])
      max_depth[i] <- max(max_depth[tree$edge[edges,2]] + 1)
      if (length(edges) == 1) {
        max_coph[i] <- max_coph[tree$edge[edges,2]]
      } else {
        max_coph[i] <- max(
          max_coph[tree$edge[edges,2]],
          sum(
            sort(
              max_length[tree$edge[edges,2]] + tree$edge.length[edges],
              decreasing = TRUE
            )[1:2]
          )
        )
      }
    }
  }
  tibble::tibble(
    max_depth,
    max_length,
    max_coph
  )
}

tip_depth <- function(tree) {
  d <- integer(ape::Ntip(tree) + ape::Nnode(tree))
  d[ape::Ntip(tree) + 1] <- 0
  for (i in seq.int(ape::Ntip(tree) + 1, length(d))) {
    edges <- which(tree$edge[,1] == i)
    for (j in edges) {
      d[tree$edge[j, 2]] <- d[i] + 1
    }
  }
  d[1:ape::Ntip(tree)]
}

depth_order <- function(tree, reverse = FALSE) {
  stopifnot(ape::is.binary(tree))
  tree <- reorder(tree, "postorder")
  mc <- max_cophenetic(tree)
  kids <- lapply(
    phangorn::Children(tree, seq_len(ape::Nnode(tree)) + ape::Ntip(tree)),
    sort
    )
  n_rotations <- 0
  on.exit(cat("\nRotated", n_rotations, "nodes.\n"))
  for (i in seq.int(ape::Nnode(tree) + ape::Ntip(tree), 1 + ape::Ntip(tree), -1)) {
    k_depth <- mc$max_depth[kids[[i - ape::Ntip(tree)]]]
    if (xor(is.unsorted(k_depth), reverse)) {
      ape::rotate(tree, i)
      n_rotations <- n_rotations + 1
    }
    if (i %% 100 == 0) cat(".")
  }
  tree
}

correlog <- function(physeq, metric, timelag,
                     break.pts = 0:13 - 0.5,
                     cutoff = FALSE) {
  dist_eco <- phyloseq::distance(physeq, method = metric)
  dist_sp <- phyloseq::sample_data(physeq) %>%
    with(x + 30000 * as.integer(site)) %>%
    dist()
  dist_t = phyloseq::sample_data(physeq) %>%
    with(as.integer(year)) %>%
    dist()
  dist_spt = dist_sp + 100000 * (dist_t - timelag)
  vegan::mantel.correlog(
    dist_eco,
    dist_spt,
    break.pts = break.pts,
    cutoff = cutoff
  )
}
