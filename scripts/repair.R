library(drake)
tmpcache <- tempdir()
file.copy(".drake", tmpcache, recursive = TRUE)
drake_cache(file.path(tmpcache, ".drake"))$repair(force = TRUE)
unlink(".drake", recursive = TRUE)
file.copy(file.path(tmpcache, ".drake"), ".", recursive = TRUE)

