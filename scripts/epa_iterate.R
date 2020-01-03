epa_iterate <- function(subject, query, subject_tree, subject_model, iterations, threads = 1) {
  align <- c(subject, query)
  alnclass <- class(align)[1]
  unalign <- gsub("-", "", align)
  unalign <- methods::as(unalign, alnclass)
  for (i in seq_len(iterations)) {
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
      fully_resolve = TRUE
    )
    realign <- DECIPHER::AlignSeqs(
      myXStringSet = unalign,
      guideTree = graft_tree,
      iterations = 0,
      refinements = 0,
      processors = threads
    )
    query <- realign[names(query)]
    subject <- realign[names(subject)]
    if (isTRUE(all.equal(realign, align))) break
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