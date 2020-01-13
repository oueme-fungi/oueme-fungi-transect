epa_iterate <- function(subject, query, subject_tree, subject_model, iterations, threads = 1) {
  align <- c(subject, query)
  alnclass <- class(align)[1]
  unalign <- gsub("-", "", align)
  unalign <- methods::as(unalign, alnclass)
  for (i in seq_len(iterations)) {
    flog.info("Starting cycle %d of iterative EPA placement and alignment.", i)
    epa_result <- epa_ng(
      ref_msa = subject,
      query = query,
      tree = subject_tree,
      model = subject_model,
      threads = threads,
      redo = TRUE
    )
    graft_tree <- gappa_graft(
      jplace = epa_result,
      threads = threads,
      allow_file_overwriting = TRUE,
      fully_resolve = FALSE
    )
    # DECIPHER requires the guide tree to be a dendrogram object with
    # maximum height 0.5.
    # This also requires that the tree be ultrametric and dichotomous.
    # In addition, there can problems with excessive recursion which are
    # controlled by sorting the tree to have the smallest clades first
    graft_tree$root_edge <- 0
    graft_tree$edge.length <- pmax(graft_tree$edge.length, 1e-6)
    graft_tree <- ape::multi2di(graft_tree)
    graft_tree <- ape::chronoMPL(graft_tree, se = FALSE, test = FALSE)
    while (any(graft_tree$edge.length < 0)) {
      graft_tree$edge.length <- pmax(graft_tree$edge.length, 1e-6)
      graft_tree <- phytools::force.ultrametric(graft_tree, method = "extend")
    }
    #graft_tree <- depth_order(graft_tree)
    #graft_tree <- ape::rotateConstr(graft_tree, rev(graft_tree$tip.label))
    graft_stats <- max_cophenetic(graft_tree)
    graft_tree$edge.length <- graft_tree$edge.length / max(graft_stats$max_length) / 4
    graft_tree <- as.dendrogram(ape::as.hclust.phylo(graft_tree))
    realign <- DECIPHER::AlignSeqs(
      myXStringSet = unalign,
      guideTree = graft_tree,
      iterations = 0,
      refinements = 0,
      processors = threads
    )
    query <- realign[names(query)]
    subject <- realign[names(subject)]
    if (isTRUE(all.equal(realign, align))) {
      flog.info("Alignment converged after %d iterations.", i)
      break
    }
    align <- realign
  }
  epa_result <- epa_ng(
    ref_msa = subject,
    query = query,
    tree = subject_tree,
    model = subject_model,
    threads = threads,
    redo = TRUE
  )
  list(jplace = epa_result, query = query, subject = subject)
}