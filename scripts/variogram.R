# Make a variogram from two distance matrices
variogram_dist <- function(eco_dist, sp_dist, breaks) {
  tibble::tibble(
    dist = unclass(sp_dist),
    gamma = unclass(eco_dist),
    bin = cut(dist, breaks = breaks)
  ) %>%
    dplyr::filter(!is.na(bin)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      np = as.numeric(dplyr::n()),
      dist = mean(dist),
      gamma = mean(gamma),
      dir.hor = 0,
      dir.ver = 0,
      id = "var1"
    ) %>%
    as.data.frame() %>%
    `class<-`(c("gstatVariogram", class(.)))
}

# Make a variogram from two distance matrices
variogramST_dist <- function(eco_dist, sp_dist, timelag, breaks) {
  tibble::tibble(
    dist = unclass(sp_dist),
    gamma = unclass(eco_dist),
    timelag = timelag,
    bin = cut(dist, breaks = breaks)
  ) %>%
    dplyr::filter(!is.na(bin)) %>%
    dplyr::group_by(bin, timelag) %>%
    dplyr::summarise(
      np = as.numeric(dplyr::n()),
      dist = mean(dist),
      gamma = mean(gamma),
      id = "var1"
    ) %>%
    dplyr::group_by(bin) %>%
    dplyr::mutate(avgDist = stats::weighted.mean(dist, np)) %>%
    dplyr::ungroup() %>%
    tidyr::extract(bin, c("min", "max"), regex = "\\(([-.0-9]+),([-.0-9]+)]") %>%
    dplyr::mutate(
      min = as.numeric(min),
      max = as.numeric(max),
      spacelag = (min + max)/2
    ) %>%
    as.data.frame() %>%
    `class<-`(c("StVariogram", class(.)))
}
