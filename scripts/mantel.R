
# rename the rows and columns of a sequence abundance table
# convert rownames from a fileame to {Seq.Run}{Plate}{Well}
# convert colnames from the full sequence to an 8-character hash.
# also pools forward and reverse reads from same well.

relabel_seqtable <- function(seqtable, hashlen = 8) {
  seqtable %>%
    # set the new column names
    magrittr::set_colnames(seqhash(colnames(.), len = hashlen)) %>%
    # convert to a tibble for easier name column operations
    tibble::as_tibble(rownames = "file") %>%
    # parse the filename
    tidyr::extract(
      col = "file",
      into = c("Seq.Run", "Plate", "Well", "Direction", "Region"),
      regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.qfilt\\.fastq\\.gz") %>%
    # create the ID
    tidyr::unite("ID", Seq.Run, Plate, Well, sep = "") %>%
    # remove unnecessary columns
    dplyr::select(-Direction, -Region) %>%
    # pool all reads from the same well
    dplyr::group_by(ID) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::ungroup() %>%
    # convert back to a matrix
    tibble::column_to_rownames("ID") %>%
    as.matrix
}

assemble_physeq <- function(platemap, datasets, seqtable, tree) {
  samp <- platemap %>%
    dplyr::mutate_at("Primer.Pair", tolower) %>%
    dplyr::left_join(
      datasets %>%
        dplyr::select(Dataset, Seq.Run, Tech, Forward, Reverse) %>%
        dplyr::mutate(
          Primer.Pair = paste(stringr::str_replace(Forward, "_tag.*$", ""),
                              stringr::str_replace(Reverse, "_tag.*$", ""),
                              sep = "_")),
      by = "Primer.Pair") %>%
    dplyr::mutate_at(c("Year", "Site", "Qual", "Plate", "Well", "Primer.Pair"),
                     factor) %>%
    tidyr::unite("ID", Seq.Run, Plate, Well, sep = "", remove = FALSE) %>%
    tibble::column_to_rownames("ID") %>%
    phyloseq::sample_data()
  
  asvtab <- seqtable %>%
    phyloseq::otu_table(taxa_are_rows = FALSE)
  
  phyloseq::phyloseq(samp, asvtab, tree)
}