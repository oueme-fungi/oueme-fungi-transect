# functions for mantel correlograms
# author Brendan Furneaux

correlog <- function(physeq, metric, timelag,
                     break.pts = 0:13 - 0.5,
                     cutoff = FALSE) {
  dist_eco <- phyloseq::distance(physeq, method = metric)
  dist_sp <- phyloseq::sample_data(physeq) %>%
    with(x + 30000 * as.integer(site)) %>%
    dist()
  dist_t = phyloseq::sample_data(physeq) %>%
    with(as.integer(year)) %>%
    dist()
  dist_spt = dist_sp + 100000 * (dist_t - timelag)
  vegan::mantel.correlog(
    dist_eco,
    dist_spt,
    break.pts = break.pts,
    cutoff = cutoff
  )
}
