library(drake)
library(futile.logger)
load("data/plan/drake.Rdata")

flog.info("Copying cache to temporary location.")
tmpcache <- tempdir()
file.copy(".drake", tmpcache, recursive = TRUE)
dc <- drake_cache(file.path(tmpcache, ".drake"))

cached <- dc$list()
flog.info("%d cached targets", length(cached))

cached <- cached[!startsWith(cached, "p-")]
cached <- cached[!startsWith(cached, "n-")]
flog.info("%d cached targets", length(cached))
cached <- setdiff(cached, plan$target)
for (t in plan$target[vapply(plan$dynamic, is.language, TRUE)]) {
  cached <- setdiff(cached, grepl(paste0(t, "_[0-9a-f]{8}"), cached))
} 
flog.info("Checking current disk space taken by cache...")
space <- system2("du", args = c("-sh", dc$path), stdout = TRUE)
flog.info(space)
if (length(cached) > 0) {
  flog.info("Removing %d cached targets no longer in plan", length(cached))
  clean(list = cached, cache = dc, garbage_collection = TRUE, verbose = 2L)
  flog.info("Checking disk space taken by cache after cleaning...")
  space <- system2("du", args = c("-sh", dc$path), stdout = TRUE)
  flog.info(space)
} else {
  flog.info("No unused targets to remove.")
}

flog.info("Copying back from temporary directory.")
unlink(".drake", recursive = TRUE)
file.copy(file.path(tmpcache, ".drake"), ".", recursive = TRUE)
