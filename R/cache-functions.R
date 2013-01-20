resetCache <- function(env) {  
  if(class(env$cache)!="function") {
    files <- unlist(c(env$cache$spectra[!sapply(env$cache$spectra, is.null)], 
                      env$cache$xic[!sapply(env$cache$xic, is.null)]))
    if(length(files)>0) file.remove(files)
  }
  env$cache <- list()
  n <- length(env$spIndex())
  env$cache$spectra <- vector("list", length=n)
  env$cache$xic <- vector("list", length=n)
}

