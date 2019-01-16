options(packrat.dependency.discovery.disabled = TRUE)

local({
  r <- getOption("repos")             # hard code the US repo for CRAN
  r["CRAN"] <- "https://ftp.acc.umu.se/mirror/CRAN/"
  options(repos = r)
})
#### -- Packrat Autoloader (version 0.5.0) -- ####
source("packrat/init.R", echo = TRUE)
#### -- End Packrat Autoloader -- ####
local({
  s <- packrat::status()
  if (is.null(s) ||
      any(is.na(s$packrat.version)) ||
      any(is.na(s$library.version)) ||
      any(s$packrat.version != s$library.version)) {
    packrat::restore()
  }
})