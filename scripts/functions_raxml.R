# function to interface with RAxML
# author Brendan Furneaux, modified from raxml function in ips package version
# 0.0.11 by C. Heibl

#' Call raxml from R
#'
#' @source modified from raxml function in ips package version 0.0.11 by C. Heibl
#' #'@title Maximum Likelihood Tree Estimation with RAxML
#'@description Provides an interface to the C program \bold{RAxML} (see
#'  Reference section) for maximum likelihood estimation of tree topology and/or
#'  branch lengths, rapid and conventional non-parametric bootstrapping, mapping
#'  splits onto individual topologies, and a lot more. See the RAxML manual for
#'  details, especially if you are a new user of RAxML.
#'@param RNAaln A matrix of RNA sequences of class
#'  \code{\link[Biostrings]{RNAMultipleStringSet}}.
#'@param m A vector of mode \code{"character"} defining a model of molecular
#'  evolution; currently only GTR model available.
#'@param f A vector of mode \code{"character"} selecting an RAxML algorithm
#'  analogous to the \code{-f} flag (see Detail section and RAxML manual).
#'@param N Either of mode \code{"integer"} or \code{"character"}. Integers give
#'  the number of independent searches on different starting tree or replicates
#'  in bootstrapping. Alternatively, one of four bootstopping criteria can be
#'  chosen: \code{"autoFC"}, \code{"autoMR"}, \code{"autoMRE"}, or
#'  \code{"autoMRE_IGN"}.
#'@param p Integer, setting a random seed for the parsimony starting trees.
#'@param b Integer, setting a random seed for bootstrapping.
#'@param x Integer, setting a random seed for rapid bootstrapping.
#'@param k Logical, if \code{TRUE}, the branch lengths of bootstrapped trees are
#'  recorded.
#'@param S A vector of mode \code{"character"} giving the secondary structure
#'  associated with the alignment.
#'@param A A vector of mode \code{"character"} defining a secondary structure
#'  substitution model.
#'@param weights A vector of mode \code{"numeric"} giving integers to assign
#'  individual weights to each column of the alignment. (-a)
#'@param partitions A data frame giving the partitions of the alignment.
#'@param outgroup A vector of mode \code{"character"} containing the names of
#'  the outgroup taxa.
#'@param backbone A \code{\link{phylo}} object representing a backbone tree.
#'@param file A vector of mode \code{"character"} giving a name for the output
#'  files.
#'@param dir Directory to generate output files in
#'@param exec A vector of mode \code{"character"} giving the path to the
#'  directory containing the RAxML executable. The default value will work on
#'  Mac OS X if the folder containing the executable is renamed to
#'  \code{"RAxML-8.0.3"}.
#'@param threads Integer, giving the number of parallel threads to use (PTHREADS
#'  only).
#'@details There are some limitations of this wrapper compared to RAxML run
#'  directly from the command line. \enumerate{ \item Only DNA is allowed as
#'  data type. \item Option \code{f} can only take a limited number of values
#'  (\code{d}, \code{a}). } % close enumerate
#'
#'  RAxML needs the specification of random seeds for parsimony estimation of
#'  starting trees and for bootstrap resampling. The corresponding argument
#'  names in \code{raxml} are identical to the flags used by RAxML (\code{-p},
#'  \code{-b}, and \code{-x}). If you choose not to give any values,
#'  \code{raxml} will generate a (different) value for each required random
#'  seed every time it is called. Be aware that \code{\link{set.seed}} will work
#'  only for \code{p}, but not for \code{b} or \code{x}.
#'@return A list with a variable number of elements, depending on the analysis
#'  chosen: \tabular{ll}{ \code{"info"} \tab RAxML log file as character
#'  string\cr \code{"bestTree"} \tab MLE of tree\cr \code{"bipartitions"} \tab
#'  MLE of tree annotated with bootstrap proportions\cr \code{"bootstrap"} \tab
#'  bootstrapped trees\cr }
#'@references (in chronolocigal order)
#'
#'  Stamatakis, A., T. Ludwig and H. Meier. 2004. RAxML-III: A fast program for
#'  maximum likelihood-based inference of large phylogenetic trees.
#'  \emph{Bioinformatics} \bold{1}: 1--8.
#'
#'  Stamatakis, A. 2006. RAxML-VI-HPC: Maximum likelihood-based phylogenetic
#'  analyses with thousands of taxa and mixed models. \emph{Bioinformatics}
#'  \bold{22}: 2688--2690.
#'
#'  Stamatakis, A., P. Hoover, and J. Rougemont. 2008. A rapid bootstrap
#'  algorithm for the RAxML web-servers. \emph{Syst. Biol.} \bold{75}: 758--771.
#'
#'  Pattengale, N. D., M. Alipour, O. R. P. Bininda-Emonds, B. M. E. Moret, and
#'  A. Stamatakis. 2010. How many bootstrap replicates are necessary?
#'  \emph{Journal of Computational Biology} \bold{17}: 337-354.
#'
#'  Stamatakis, A. 2014. RAxML Version 8: A tool for phylogenetic analysis and
#'  post-analysis of large phylogenies. \emph{Bioinformatics} Advance Access.
#' @note RAxML is a C program and the source code is not contained in this
#'   package. This means that in order to run this function you will need to
#'   install RAxML yourself. See
#'   \url{http://sco.h-its.org/exelixis/web/software/raxml/} for the most recent
#'   documentation and source code of RAxML. Depending on where you chose to
#'   install RAxML, you need to adjust the \code{exec} argument.
#'
#'  \code{raxml} was last tested and running fine on Mac OS X with RAxML 8.0.29.
#'  Please be aware that calling third-party software from within R is a
#'  platform-specific process and I cannot guarantee that \code{raxml} will
#'  behave properly on any system.
#'@seealso \code{\link{raxml.partitions}} to store partitioning information in a
#'  data frame suitable for input as \code{partitions} argument in \code{raxml}.
raxml_RNA <- function(RNAaln, m = "GTRCAT", f, N, p, b, x, k, S, A = "S16",
                      weights, partitions, outgroup, backbone = NULL,
                      file = paste0("fromR_", Sys.Date()), exec,
                      dir = tempdir(), threads){

  if (!inherits(RNAaln, "RNAMultipleAlignment")) stop("RNAaln is not of class
                                                      'RNAMultipleAlignment'")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  olddir <- setwd(dir)
  on.exit(setwd(olddir))

  # number of threads (PTHREADS only)
  # ---------------------------------
  if (!missing(threads))
    exec <- paste(exec, "-T", threads)

  ## input file names
  ## ----------------
  rin <- c(s = paste("-s ", file, ".phy", sep = ""),
           n = paste("-n", file),
           napt = paste("-n ", file, ".APT", sep = ""),
           weight = paste0(file, "-weights.txt"),
           partition = paste0(file, "-partitions.txt"),
           secstr = paste0(file, "-secstr.txt"),
           tree = paste0(file, "-backbone.tre"))

  # Clear previous runs with same tag
  # ---------------------------------
  unlink(rin)
  unlink(list.files(pattern = paste0("RAxML_[[:alpha:]]+[.]", file)))

  # Substitution model
  # ------------------
  m <- match.arg(m, c("GTRCAT", "GTRCATX",
                      "GTRCATI", "GTRCATIX",
                      "ASC_GTRCAT", "ASC_GTRCATX",
                      "GTRGAMMA", "GTRGAMMAX",
                      "GTRGAMMAI", "GTRGAMMAIX",
                      "ASC_GTRGAMMA", "ASC_GTRGAMMAX"))
  m <- paste("-m", m)

  # Secondary structure model
  # ------------------
  A <- match.arg(A, c("S6A", "S6B", "S6C", "S6D", "S6E",
                      "S7A", "S7B", "S7C", "S7D", "S7E", "S7F",
                      "S16", "S16A", "S16B"))
  A <- paste("-A", A)

  ## Number of searches/replicates
  ## -----------------------------
  if (is.character(N)){
    N <- match.arg(N, c("autoFC", "autoMR", "autoMRE", "autoMRE_IGN"))
  }
  N <- paste("-N", N)

  ## Random seeds
  ## ------------
  rs <- function(rseed, type = "p"){
    if (missing(rseed)) rseed <- sample(1:999999, 1)
    paste0("-", type, " ",  rseed)
  }
  p <- rs(p)
  if (!missing(x)) x <- rs(x, type = "x")

  ## Write sequences to input file
  ## -----------------------------
  Biostrings::write.phylip(RNAaln, paste(file, "phy", sep = "."))

  ## rout: raxml output file names
  ## -----------------------------
  output.types <- c("info", "bestTree", "bootstrap", "bipartitions")
  rout <- paste("RAxML_", output.types, ".", file, sep = "")
  names(rout) <- output.types

  ## Algorithms
  ## ----------
  if (missing(f)) f <- "d"
  f <- match.arg(f, c("d", "a"))
  alg <- paste("-f", f, p) # add parsimony seed to algorithm
  if (f == "a") alg <- paste(alg, x)
  if (!missing(b)) alg <- paste(alg, "-b", b)
  if (missing(N)) stop("the number of runs must be given (N)")

  ## Outgroup
  ## --------
  if (missing(outgroup)){
    o <- ""
  } else {
    if (length(grep(",", outgroup))) outgroup <- unlist(strsplit(outgroup, ","))
    o <- outgroup %in% rownames(RNAaln)
    if (!all(o)){
      o <- paste(paste("\n  -", outgroup[!o]), collapse = "")
      stop(paste("outgroup names not in 'DNAbin':",   o))
    }
    o <- paste(outgroup, collapse = ",")
    o <- paste("-o", o)
  }

  ## Columns weights
  ## ---------------
  if (!missing(weights)){
    if (length(weights) != ncol(RNAaln)) stop("number of 'weights' don't equal number of columns")
    write(weights, rin["weight"], length(weights))
    weights <- paste("-a", rin["weight"])
  } else {
    weights <- ""
  }

  # Write partition file
  ## ------------------
  if (!missing(partitions)) {
    if (is.character(partitions)){
      q <- partitions
    } else {
      q <- paste(partitions$type, ", ",
                 partitions$locus, " = ",
                 partitions$begin, "-",
                 partitions$end, sep = "")
    }
    write(q, rin["partition"])
    multipleModelFileName <- paste(" -q", rin["partition"], "")
  } else multipleModelFileName <- ""

  ## Backbone tree
  ## -------------
  if (!is.null(backbone)){
    ape::write.tree(backbone, rin["tree"])
    g <- paste(" -g", rin["tree"], "")
  } else {
    g <- " "
  }

  ## Secondary structure
  ## -------------
  if (!missing(S) && !is.null(S)) {
    writeLines(S, rin["secstr"])
    S <- paste("-S", rin["secstr"], A)
  } else {
    S <- ""
  }

  ## Save branch lengths of bootstrap replicates
  ## -------------------------------------------
  if (missing(k)) k <- FALSE
  k <- ifelse(k, "-k", "")

  ## Prepare and execute call
  ## ------------------------
  CALL <- paste(exec, alg, m, o, k,
                weights,
                multipleModelFileName, N, g,
                S,
                rin["s"], rin["n"])

  if (length(grep("MPI", exec))) system(paste("mpirun", CALL))
  print(CALL)
  system(CALL)

  res <- scan(rout["info"], quiet = TRUE, what = "char", sep = "\n")
  if (length(grep("exiting", res)))
    stop("\n", paste(res, collapse = "\n"))

  ## Read results
  ## ------------
  bestTree <- bipartitions <- bootstrap <- NULL
  info <- scan(rout["info"], what = "c", sep = "\n", quiet = TRUE)
  if (f %in% c("a", "d") & missing(b))
    bestTree <- ape::read.tree(rout["bestTree"])
  if (f %in% c("a"))
    bipartitions <- ape::read.tree(rout["bipartitions"])
  if (!missing(b) | !missing(x))
    bootstrap <- ape::read.tree(rout["bootstrap"])

  obj <- list(info = info,
              bestTree = bestTree,
              bipartitions = bipartitions,
              bootstrap = bootstrap)
  obj[sapply(obj, is.null)] <- NULL
  obj
}
