# This file generally does not need to be rerun, but I have needed to reset
# Packrat enough times that I'm sick of typing it all out.
packrat::init(infer.dependencies = FALSE)
setRepositories(graphics = FALSE)
install.packages(c("tidyverse",
                   "here",
                   "assertr",
                   "clustermq",
                   "remotes",
                   "tidyxl",
                   "glue",
                   "bookdown",
                   "tikzDevice",
                   "seqinr",
                   "visNetwork",
                   "networkD3",
                   "odseq",
                   "DECIPHER",
                   "tictoc",
                   "phyloseq",
                   "iNEXT",
                   "BiocGenerics",
                   "BiocParallel",
                   "Biostrings",
                   "Rsamtools",
                   "GenomicAlignments",
                   "ggtree"), Ncpus = 3)

remotes::install_github(c("brendanf/ShortRead",
                          "brendanf/rITSx",
                          "brendanf/FUNGuildR",
                          "benjjneb/dada2",
                          "ropensci/drake"), Ncpus = 3)
packrat::snapshot()
