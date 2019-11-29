epa_ng <- function(ref_msa, tree, query, outdir = tempdir(), model, threads, exec = "epa-ng") {
  
  ref_msa_file <- tempfile("reference", fileext = ".fasta")
  if (methods::is(ref_msa, "XStringSet")) {
    assertthat::assert_that(length(unique(Biostrings::width(ref_msa))) == 1)
  } else if (methods::is(ref_msa, "MultipleAlignment")) {
    ref_msa <- methods::as(ref_msa, "XStringSet")
  } else if (is.character(ref_msa)) {
    if (length(ref_msa) == 1 && file.exists(ref_msa)) {
      ref_msa_file <- ref_msa
    }
    assertthat::assert_that(length(unique(nchar(ref_msa))) == 1)
    if (tzara::has_alphabet(ref_msa, Biostrings::DNA_ALPHABET)) {
      ref_msa <- Biostrings::DNAStringSet(ref_msa)
    } else if (tzara::has_alphabet(ref_msa, Biostrings::RNA_ALPHABET)) {
      ref_msa <- Biostrings::RNAStringSet(ref_msa)
    } else if (tzara::has_alphabet(ref_msa, Biostrings::AA_ALPHABET)) {
      ref_msa <- Biostrings::AAStringSet(ref_msa)
    } else {
      stop("Unknown alphabet in reference alignment.")
    }
  } else {
    stop("'ref_msa' should be an XStringSet, MultipleAlignment, character vector",
         "of aligned sequences, or filename.")
  }
  if (methods::is(ref_msa, "XStringSet")) {
    Biostrings::writeXStringSet(ref_msa, ref_msa_file)
    on.exit(file.remove(ref_msa_file))
  }
  
  query_file <- tempfile("query", fileext = ".fasta")
  if (methods::is(query, "XStringSet")) {
    assertthat::assert_that(length(unique(Biostrings::width(query))) == 1)
  } else if (methods::is(query, "MultipleAlignment")) {
    ref_msa <- methods::as(ref_msa, "XStringSet")
  } else if (is.character(query)) {
    if (length(query) == 1 && file.exists(query)) {
      query_file <- query
    }
    assertthat::assert_that(length(unique(nchar(query))) == 1)
    if (tzara::has_alphabet(query, Biostrings::DNA_ALPHABET)) {
      query <- Biostrings::DNAStringSet(query)
    } else if (tzara::has_alphabet(query, Biostrings::RNA_ALPHABET)) {
      query <- Biostrings::RNAStringSet(query)
    } else if (tzara::has_alphabet(query, Biostrings::AA_ALPHABET)) {
      query <- Biostrings::AAStringSet(query)
    } else {
      stop("Unknown alphabet in query alignment.")
    }
  } else {
    stop("'query' should be an XStringSet, MultipleAlignment, character vector",
         "of aligned sequences, or filename.")
  }
  if (methods::is(query, "XStringSet")) {
    Biostrings::writeXStringSet(query, query_file)
  }
  
  if (methods::is(tree, "phylo")) {
    tree_file <- tempfile("tree", tempdir(), ".tree")
    ape::write.tree(tree, tree_file)
    on.exit(file.remove(tree_file))
  } else if (is.character(tree) && file.exists(tree)) {
    tree_file <- tree
  } else {
    stop("'tree' should be a phylo object or a file name.")
  }
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  
  args <- c(
    "--ref-msa", ref_msa_file,
    "--query", query_file,
    "--tree", tree_file,
    "--outdir", outdir
  )
  
  for (m in model) {
    args <- c(
      args,
      "--model", model
    )
  }
  
  if (!missing(threads) && !is.null(threads)) {
    assertthat::assert_that(assertthat::is.count(threads))
    args <- c(
      args,
      "--threads", threads
    )
  }
  
  system2(exec, args = args)
}