# Run EPA-NG to place aligned short reads onto a tree based on long sequences
epa_ng <- function(ref_msa, tree, query, outdir = tempdir(), model, threads, exec = "epa-ng", redo = FALSE) {
  
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
    on.exit(unlink(outdir, recursive = TRUE))
  }
  
  if (length(model) == 1 && file.exists(model)) {
    model_file <- model
  } else {
    model_file <- tempfile("info")
    on.exit(unlink(model_file))
    writeLines(model, model_file)
  }
  
  args <- c(
    "--ref-msa", ref_msa_file,
    "--query", query_file,
    "--tree", tree_file,
    "--outdir", outdir,
    "--model", model_file
  )
  
  if (isTRUE(redo)) {
    args <- c(args, "--redo")
  }
  
  if (!missing(threads) && !is.null(threads)) {
    assertthat::assert_that(assertthat::is.count(threads))
    args <- c(
      args,
      "--threads", threads
    )
  }
  
  system2(exec, args = args)
  jsonlite::read_json(file.path(outdir, "epa_result.jplace"))
}

is_jplace <- function(jplace) {
  is.list(jplace) &&
    setequal(
      names(jplace),
      c("tree", "placements", "metadata", "version", "fields")
    )
}

# run "gappa examine graft" to add short reads to a tree based on placement
# results (e.g. from EPA)
gappa_graft <- function(jplace, outdir = tempdir(), threads = NULL,
                        allow_file_overwriting = FALSE, verbose = FALSE,
                        fully_resolve = FALSE) {
  if (is.character(jplace) && all(file.exists(jplace))) {
    jplace_file <- jplace
  } else if (is_jplace(jplace)) {
    jplace_file <- tempfile("result", fileext = ".jplace")
    jsonlite::write_json(jplace, jplace_file, auto_unbox = TRUE)
    on.exit(unlink(jplace_file))
  }
  
  out_file <- sub(".jplace", ".newick", jplace_file)
  
  args <- c(
    "examine", "graft",
    "--jplace-path", jplace_file,
    "--out-dir", outdir
  )
  
  if (isTRUE(fully_resolve)) {
    args <- c(args, "--fully-resolve")
  }
  
  if (!is.null(threads)) args <- c(args, "--threads", threads)
  if (isTRUE(allow_file_overwriting)) args <- c(args, "--allow-file-overwriting")
  if (isTRUE(verbose)) args <- c(args, "--verbose")
  
  if (missing(outdir)) on.exit(unlink(out_file))
  system2("gappa", args)
  ape::read.tree(out_file)
}

# Delete the branch leading to sequences placed on a tree using EPA-NG and
# GAPPA, so that the tree can be used as a guide tree for RAxML to determine
# relationships between the short reads.
grafts_to_polytomies <- function(graft_tree, base_tree) {
  allDesc <- phangorn::allDescendants(graft_tree)
  endtips <- allDesc %>%
    purrr::map_lgl(~all(. <= ape::Ntip(graft_tree)) &
                     length(.) > 1) %>%
    which()
  grafts <- purrr::keep(
    endtips,
    ~!any(graft_tree$tip.label[allDesc[[.]]] %in% base_tree$tip.label)
  )
  graft_tree$edge.length[graft_tree$edge[,2] %in% grafts] <- 0
  ape::di2multi(graft_tree)
}
