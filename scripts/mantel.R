
# rename the rows and columns of a sequence abundance table
# convert rownames from a fileame to {Seq.Run}{Plate}{Well}
# convert colnames from the full sequence to an 8-character hash.
# also pools forward and reverse reads from same well.

relabel_seqtable <- function(seqtable) {
  seqtable %>%
    # set the new column names
    magrittr::set_colnames(tzara::seqhash(colnames(.))) %>%
    # convert to a tibble for easier name column operations
    tibble::as_tibble(rownames = "file") %>%
    # parse the filename
    tidyr::extract(
      col = "file",
      into = c("seq_run", "plate", "well", "direction", "region"),
      regex = "([:alpha:]+_\\d+)_(\\d+)-([A-H]1?\\d)([fr]?)-([:alnum:]+)\\.qfilt\\.fastq\\.gz") %>%
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

assemble_physeq <- function(platemap, datasets, seqtable, tree) {
  samp <- platemap %>%
    dplyr::mutate_at("primer_pair", tolower) %>%
    dplyr::left_join(
      datasets %>%
        dplyr::select(dataset, seq_run, tech, forward, reverse) %>%
        dplyr::mutate(
          primer_pair = paste(stringr::str_replace(forward, "_tag.*$", ""),
                              stringr::str_replace(reverse, "_tag.*$", ""),
                              sep = "_")),
      by = "primer_pair") %>%
    dplyr::mutate_at(c("year", "site", "qual", "plate", "well", "primer_pair"),
                     factor) %>%
    tidyr::unite("ID", seq_run, plate, well, sep = "", remove = FALSE) %>%
    tibble::column_to_rownames("ID") %>%
    phyloseq::sample_data()
  
  asvtab <- seqtable %>%
    phyloseq::otu_table(taxa_are_rows = FALSE)
  
  phyloseq::phyloseq(samp, asvtab)#, tree)
}