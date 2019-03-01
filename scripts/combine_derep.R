#' Combine DADA2 Derep Objects to Get a Master Mapping to Unique Sequences
#'
#' @param dereps A list of tibbles, with columns "file" and "derep",
#' representing the source file and the derep object, respectively.
#'
#' @return a list with two members:
#' \describe {
#' \item{map}{a tibble with columns "file", "idx", and "newmap", giving the mapping from the "idx"th sequence in "file" to a sequence in "fasta"}
#' \item{fasta}{a DNAStringSet giving all unique sequences; the name of the sequence is an integer which matches the values of "newmap" in "map"}
#' }
combine_derep <- function(...) {
  dereps <- dplyr::bind_rows(...)
  
  # get all the old mappings
  oldmap <- dereps %>%
    dplyr::mutate_at("derep", purrr::map, "map") %>%
    tidyr::unnest() %>%
    dplyr::rename(oldmap = derep) %>%
    dplyr::group_by(file) %>%
    dplyr::mutate(idx = seq_along(file)) %>%
    dplyr::ungroup()
  
  # get the old unique sequences
  olduniques <- dereps %>%
    dplyr::mutate_at("derep",
                     ~purrr::map(., .f = ~tibble::tibble(seq = names(.$uniques),
                                                         n = .$uniques))) %>%
    tidyr::unnest()
  
  # combine duplicate sequences among all files.
  newuniques <- olduniques %>%
    dplyr::group_by(seq) %>%
    dplyr::summarize(n = sum(n)) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::transmute(seq = seq,
              newmap = seq_along(seq))
  
  # create a mapping from the unique sequences list in each file
  # to the master unique sequence list
  newderep <- olduniques %>%
    {dplyr::left_join(
      dplyr::select(., file, seq) %>%
        dplyr::group_by(file) %>%
        dplyr::mutate(oldmap = seq_along(seq)) %>%
        dplyr::ungroup(),
      newuniques,
      .by = "seq")}
  
  out <- list()
  # map from the individual sequences in each file to the master unique
  # sequence list
  out$map <- dplyr::left_join(oldmap,
                       dplyr::select(newderep, "file", "oldmap", "newmap"),
                       by = c("file", "oldmap"))%>%
    dplyr::select(file, "idx", "newmap")
  #unique sequence list
  out$fasta <- Biostrings::DNAStringSet(x = set_names(newuniques$seq,
                                                      newuniques$newmap))
  out
}