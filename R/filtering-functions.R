saveFilterValues <- function(env) {
  filterValues <- list(rtfrom=svalue(filterInfo$rtfrom), 
                       rtto=svalue(filterInfo$rtto), 
                       pmzfrom=svalue(filterInfo$pmzfrom), 
                       pmzto=svalue(filterInfo$pmzto), 
                       spifrom=svalue(filterInfo$spifrom), 
                       spito=svalue(filterInfo$spito), 
                       indexfrom=svalue(filterInfo$indexfrom), 
                       indexto=svalue(filterInfo$indexto),
                       pcfrom=svalue(filterInfo$pcfrom), 
                       pcto=svalue(filterInfo$pcto),
                       massfrom=svalue(filterInfo$massfrom), 
                       massto=svalue(filterInfo$massto),
                       ms1=svalue(filterInfo$ms1), 
                       ms2=svalue(filterInfo$ms2))
  lapply(filterValues, as.numeric)
}

filterReset <- function(env) {
  mapply(blockHandler, filterInfo, filterSpectraHandlerIDs)
  svalue(filterInfo$rtactive) <- FALSE
  svalue(filterInfo$rtfrom) <- min(spRtime())
  svalue(filterInfo$rtto) <- max(spRtime())
  svalue(filterInfo$pmzactive) <- FALSE
  svalue(filterInfo$pmzfrom) <- min(spPrecMz())
  svalue(filterInfo$pmzto) <- max(spPrecMz())
  svalue(filterInfo$spiactive) <- FALSE
  svalue(filterInfo$spifrom) <- min(spPrecInt())
  svalue(filterInfo$spito) <- max(spPrecInt())
  svalue(filterInfo$indexactive) <- FALSE
  svalue(filterInfo$indexfrom) <- min(spIndex())
  svalue(filterInfo$indexto) <- max(spIndex())
  svalue(filterInfo$pcactive) <- FALSE
  svalue(filterInfo$pcfrom) <- min(spPrecCharge())
  svalue(filterInfo$pcto) <- max(spPrecCharge())
  svalue(filterInfo$massactive) <- FALSE
  svalue(filterInfo$massfrom) <- min(spPrecMz()*spPrecCharge())
  svalue(filterInfo$massto) <- max(spPrecMz()*spPrecCharge())
  svalue(filterInfo$ms1) <- TRUE
  svalue(filterInfo$ms2) <- TRUE
  filterValues <<- saveFilterValues(env)
  mapply(unblockHandler, filterInfo, filterSpectraHandlerIDs)
}

filterSwitch <- function(on) {
  lapply(filterInfo, function(x) enabled(x) <- on)
}


btwn <- function(x, from, to) {
  stopifnot(all(is.numeric(x), is.character(from), is.character(to)))
  x >= as.numeric(from) & x <= as.numeric(to)
}

validityCheck <- function(from, to, data, prevFrom, prevTo) {
  numFrom <- as.numeric(svalue(from))
  numTo <- as.numeric(svalue(to))
  prevValid <- prevFrom >= min(data) & prevFrom < prevTo & prevTo > min(data) & prevTo <= max(data)
  if(numFrom < min(data)) svalue(from) <- min(data)
  if(numFrom > max(data)) svalue(from) <- ifelse(prevValid, prevFrom, min(data))
  if(numTo < min(data)) svalue(to) <- ifelse(prevValid, prevTo, max(data))
  if(numTo > max(data)) svalue(to) <- max(data)
  if(numFrom > numTo) {
    if(numFrom!=prevFrom) svalue(from) <- prevFrom
    else svalue(to) <- prevTo
  }
}

filterSpectra <- function(h, ...) {
  
  keep <- TRUE
  
  mapply(blockHandler, filterInfo, filterSpectraHandlerIDs)
  
  validityCheck(filterInfo$rtfrom, filterInfo$rtto, 
                spRtime(), filterValues$rtfrom, filterValues$rtto)
  validityCheck(filterInfo$pmzfrom, filterInfo$pmzto, 
                spPrecMz(), filterValues$pmzfrom, filterValues$pmzto)
  validityCheck(filterInfo$spifrom, filterInfo$spito, 
                spPrecInt(), filterValues$spifrom, filterValues$spito)
  validityCheck(filterInfo$indexfrom, filterInfo$indexto, 
                spIndex(), filterValues$indexfrom, filterValues$indexto)
  validityCheck(filterInfo$pcfrom, filterInfo$pcto, 
                spPrecCharge(), filterValues$pcfrom, filterValues$pcto)
  validityCheck(filterInfo$massfrom, filterInfo$massto, 
                spPrecMz()*spPrecCharge(), filterValues$massfrom, filterValues$massto)
  
  if(svalue(filterInfo$rtactive)) 
    keep <- keep & btwn(spRtime(), 
                        svalue(filterInfo$rtfrom), svalue(filterInfo$rtto))    
  if(svalue(filterInfo$pmzactive)) 
    keep <- keep & (btwn(spPrecMz(), 
                        svalue(filterInfo$pmzfrom), svalue(filterInfo$pmzto)) | spMsLevel() == 1)
  if(svalue(filterInfo$spiactive)) 
    keep <- keep & btwn(spPrecInt(), 
                        svalue(filterInfo$spifrom), svalue(filterInfo$spito))
  if(svalue(filterInfo$indexactive)) 
    keep <- keep & btwn(spIndex(), 
                        svalue(filterInfo$indexfrom), svalue(filterInfo$indexto))
  if(svalue(filterInfo$pcactive)) 
    keep <- keep & btwn(spPrecCharge(), 
                        svalue(filterInfo$pcfrom), svalue(filterInfo$pcto))
  if(svalue(filterInfo$massactive)) 
    keep <- keep & btwn(spPrecCharge(), 
                        svalue(filterInfo$massfrom), svalue(filterInfo$massto))
  
  # If MS1 was deselected previously, and user deselects MS2, 
  # select MS1 in order to have something selected
  if(!svalue(filterInfo$ms1) & !filterValues$ms2) svalue(filterInfo$ms2) <- TRUE    
  # Likewise with MS2
  if(!svalue(filterInfo$ms2) & !filterValues$ms1) svalue(filterInfo$ms1) <- TRUE 
  
  # Disable Precursos MZ filter if only MS2 selected  
  enabled(filterInfo$rtactive) <- svalue(filterInfo$ms2)
  enabled(filterInfo$rtfrom) <- svalue(filterInfo$ms2)
  enabled(filterInfo$rtto) <- svalue(filterInfo$ms2)
  
  if(!svalue(filterInfo$ms1)) keep <- keep & spMsLevel() != 1
  if(!svalue(filterInfo$ms2)) keep <- keep & spMsLevel() != 2
  
  filterValues <<- saveFilterValues(env)
  
  mapply(unblockHandler, filterInfo, filterSpectraHandlerIDs)
  
  if(any(keep)) {
    newSequence <<- spIndex()[keep]
    # Compare sequences and update if necessary
    if(!identical(currSequence, newSequence)) {
      
      prevIndex <- currSequence[counter]
      currSequence <<- newSequence
      
      # Find the nearest spectrum to the one before
      counter <<- which.min(abs(newSequence - prevIndex))
      
      # also update graphs if index has changed
      if(prevIndex!=newSequence[counter]) updateSpectrum()
    } 
  } else cat("\nNot a single spectrum survived filtering!")
# 
#   if(!identical(c(filterValues$ms1, filterValues$ms2), 
#                 c(svalue(filterInfo$ms1), svalue(filterInfo$ms1)))) plotXIC()
}
