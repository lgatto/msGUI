resetCache <- function() {  
  env$cache <- list()
  n <- length(spIndex())
  env$cache$spectra <- vector("list", length=n)
  env$cache$xic <- vector("list", length=n)
}

