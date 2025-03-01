# functions for spatiotemporal variograms
# inspired by the gstat package, but ultimately totally reimplemented
# author Brendan Furneaux

variogram_table <- function(eco_dist, sp_dist, breaks) {
  tibble::tibble(
    dist = unclass(sp_dist),
    gamma = unclass(eco_dist),
    bin = cut(dist, breaks = breaks)
  )
}

as_variogram <- function(vario) {
  assertthat::assert_that(
    is.data.frame(vario),
    assertthat::has_name(vario, "dist"),
    assertthat::has_name(vario, "gamma"),
    assertthat::has_name(vario, "bin")
  )
  vario %>%
    dplyr::filter(!is.na(bin)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      np = as.numeric(dplyr::n()),
      dist = mean(dist),
      gamma = mean(gamma),
      # var = var(gamma),
      dir.hor = 0,
      dir.ver = 0,
      id = "var1"
    ) %>%
    as.data.frame()
}

# Make a variogram from two distance matrices
variogram_dist <- function(eco_dist, sp_dist, breaks) {
  variogram_table(eco_dist, sp_dist, breaks) %>%
    as_variogram()
}

variogramST_table <- function(eco_dist, sp_dist, timelag, breaks) {
  tibble::tibble(
    dist = unclass(sp_dist),
    gamma = unclass(eco_dist),
    timelag = unclass(timelag),
    bin = cut(dist, breaks = breaks)
  )
}

as_variogramST <- function(vario) {
  assertthat::assert_that(
    is.data.frame(vario),
    assertthat::has_name(vario, "dist"),
    assertthat::has_name(vario, "timelag"),
    assertthat::has_name(vario, "gamma"),
    assertthat::has_name(vario, "bin")
  )
  vario %>%
    dplyr::filter(!is.na(bin)) %>%
    dplyr::group_by(bin, timelag) %>%
    dplyr::summarise(
      np = as.numeric(dplyr::n()),
      dist = mean(dist),
      gamma = mean(gamma),
      # var = var(gamma),
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

# Make a variogram from two distance matrices
variogramST_dist <- function(eco_dist, sp_dist, timelag, breaks) {
  variogramST_table(eco_dist, sp_dist, timelag, breaks) %>%
    as_variogramST()
}

variog <- function(physeq, metric, breaks) {
  dist_eco <- phyloseq::distance(physeq, method = metric) %>%
    c()
  dist_sp <- phyloseq::sample_data(physeq) %>%
    with(x + 30000 * as.integer(site)) %>%
    dist() %>%
    c()
  dist_t = phyloseq::sample_data(physeq) %>%
    with(as.integer(year)) %>%
    dist() %>%
    c()
  dist_spt = dist_sp + 100000 * dist_t
  variogram_table(
    eco_dist = dist_eco,
    sp_dist = dist_spt,
    breaks = breaks
  )
}

variogST <- function(physeq, metric, breaks) {
  dist_eco <- phyloseq::distance(physeq, method = metric) %>%
    c()
  dist_sp <- phyloseq::sample_data(physeq) %>%
    with(x + 30000 * as.integer(site)) %>%
    dist() %>%
    c()
  dist_t = phyloseq::sample_data(physeq) %>%
    with(as.integer(year)) %>%
    dist() %>%
    c()
  variogramST_table(
    eco_dist = dist_eco,
    sp_dist = dist_sp,
    timelag = dist_t,
    breaks = breaks
  )
}

variog_predict = function(model, newdata) {
  mutate(newdata, gamma = predict(model, newdata = newdata))
}

variog_summary = function(model) {
  dplyr::mutate(
    broom::tidy(model),
    confint = purrr::map(
      term,
      ~ tryCatch(
        broom::confint_tidy(model, parm = .),
        error = function(e) tibble::tibble(conf.low = NA_real_, conf.high = NA_real_)
      )
    )
  ) %>%
    tidyr::unnest(confint)
}

fit_variogram <- function(variog) {
    as_variogram(variog) %>%
    stats::nls(
      # this it reparameterized to give better convergence
      # sill = 1 - (C0 + C1)
      # nugget = C1 / (C0 + C1)
      gamma ~ (1 - sill) - (1 - sill) * (1 - nugget) * exp(dist/range*log(0.05)),
      data = .,
      start = list(
        nugget = min(.$gamma),
        sill = 1 - max(.$gamma),
        range = 30
      ),
      upper = list(nugget = 0.99, sill = 0.99, range = Inf),
      lower = list(nugget = 0, sill = 0, range = 0.001),
      algorithm = "port",
      weights = pmin(.$np/.$dist, 1),
      control = list(warnOnly = TRUE, maxiter = 1000, minFactor = 1/4096)
    ) %>% {
      list(
        pars = as.list(.$m$getPars()),
        summary = variog_summary(.),
        predict = variog_predict(
          .,
          newdata = expand_grid(dist = seq(1, 25, 0.5), timelag = 0:1)
        )
      )
    }
}

fit_variogramST <- function(variogST, variofit2) {
  as_variogramST(variogST) %>%
    stats::nls(
      # this it reparameterized to give better convergence
      # sill = 1 - (C0 + C1)
      # nugget = C1 / (C0 + C1)
      gamma ~ (1 - sill) - (1 - sill) * (1 - nugget) * exp((dist/range + timelag/timerange)*log(0.05)),
      # add the parameters from the spatial-only fit.
      data = .,
      start = c(
        list(timerange = 2),
        variofit2$pars
      ),
      lower = c(
        list(timerange = 0.001),
        variofit2$pars
      ),
      algorithm = "port",
      weights = pmin(.$np/.$dist, 1),
      control = list(warnOnly = TRUE, maxiter = 1000, minFactor = 1/4096)
    ) %>% {
      list(
        pars = as.list(.$m$getPars()),
        summary = variog_summary(.),
        predict = variog_predict(
          .,
          newdata = expand_grid(dist = seq(1, 25, 0.5), timelag = 0:1)
        )
      )
    }
}
