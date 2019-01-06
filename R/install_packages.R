if (!exists("force.install")) force.install <- FALSE

maybe.install <- function(package, force.install = FALSE) {
  if (!requireNamespace(package, quietly = TRUE) || force.install)
    install.packages(package)
}

maybe.install.Bioc <- function(package, force.install = FALSE, version = "3.8") {
  if (!requireNamespace(package, quietly = TRUE) || force.install)
    BiocManager::install(package, version = version)
}

maybe.install.Git <- function(package, force.install = FALSE) {
  if (!requireNamespace(strsplit(package, "/")[[1]][2], quietly = TRUE) || force.install)
    devtools::install_github(package)
}


lapply(c("magrittr",
         "tidyverse",
         "glue",
         "BiocManager",
         "readxl",
         "plyr",
         "knitr",
         "devtools",
         "seqinr"),
       maybe.install,
       force.install = force.install)

lapply(c("ShortRead",
         "Biostrings"),
       maybe.install.Bioc,
       force.install = force.install)

lapply(c("benjjneb/dada2",
         # "hoesler/rwantshue",
         "mhahsler/rBLAST",
         "hadley/multidplyr"),
       maybe.install.Git,
       force.install = force.install)
