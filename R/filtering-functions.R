saveFilterValues <- function() {
  env$filterValues <- lapply(filterInfo, function(x) lapply(x, function(i) svalue(i)))
  env$filterValuesXIC <- lapply(filterInfoXIC, function(x) lapply(x, function(i) svalue(i)))
  env$filterValuesMS <- lapply(filterInfoMS, svalue)
}

blockFilters <- function(block=TRUE) {
  if(block) {
    mapply(function(obj, id) mapply(blockHandler, obj, id), filterInfo, filterSpectraHandlerIDs)
    mapply(function(obj, id) mapply(blockHandler, obj, id), filterInfoXIC, filterXICHandlerIDs)
    mapply(blockHandler, filterInfoMS, filterSpectraMSHandlerIDs)    
  } else {    
    mapply(function(obj, id) mapply(unblockHandler, obj, id), filterInfo, filterSpectraHandlerIDs)
    mapply(function(obj, id) mapply(unblockHandler, obj, id), filterInfoXIC, filterXICHandlerIDs)
    mapply(unblockHandler, filterInfoMS, filterSpectraMSHandlerIDs)   
  }
}

updateRanges <- function(obj, values) {
  svalue(obj$active) <- FALSE
  svalue(obj$from) <- min(values)
  svalue(obj$to) <- max(values)
}
  
filterReset <- function(env) {  
  blockFilters()
  
  # Filters for spectra
  mapply(updateRanges, filterInfo, filterData)
  
  # MS levels checkboxes
  lapply(filterInfoMS, function(i) svalue(i) <- TRUE)
  
  # Filter for XIC
  mapply(updateRanges, filterInfoXIC, list(spPrecMz()))
  
  saveFilterValues()  
  blockFilters(FALSE)
}

filterSwitch <- function(on) {
  lapply(filterInfo, function(x) lapply(x, function(i) enabled(i) <- on))
  lapply(filterInfoXIC, function(x) lapply(x, function(i) enabled(i) <- on))
  lapply(filterInfoMS, function(x) enabled(x) <- on)
}

btwn <- function(x, from, to) {
  stopifnot(all(is.numeric(x), is.character(from), is.character(to)))
  x >= as.numeric(from) & x <= as.numeric(to)
}


validityCheck <- function(object, data, pastValues) {
  
  from <- object$from
  to <- object$to
  
  numFrom <- suppressMessages(as.numeric(svalue(from)))
  numTo <- suppressMessages(as.numeric(svalue(to)))
  
  prevFrom <- as.numeric(pastValues$from)
  prevTo <- as.numeric(pastValues$to)
    
  prevValid <- prevFrom >= min(data) & prevFrom < prevTo & prevTo > min(data) & prevTo <= max(data)
  if(is.na(numFrom)) svalue(from) <- min(data) else {
    if(numFrom < min(data)) svalue(from) <- min(data)
    if(numFrom > max(data)) svalue(from) <- ifelse(prevValid, prevFrom, min(data))
  }
  if(is.na(numTo)) svalue(to) <- ifelse(prevValid, prevTo, max(data)) else {
    if(numTo < min(data)) svalue(to) <- ifelse(prevValid, prevTo, max(data))
    if(numTo > max(data)) svalue(to) <- max(data)
  }
  if(as.numeric(svalue(from)) > as.numeric(svalue(to))) {
    if(numFrom!=prevFrom) svalue(from) <- prevFrom
    else svalue(to) <- prevTo
  }
}

filterSpectra <- function(h, ...) {
  blockFilters()
  
  # Check entered values, fix if needed
  mapply(validityCheck, filterInfo, filterData, filterValues)
    
  # Make sure at least one MS level is selected
  if(!svalue(filterInfoMS$ms1) & !filterValuesMS$ms2) svalue(filterInfoMS$ms2) <- TRUE
  if(!svalue(filterInfoMS$ms2) & !filterValuesMS$ms1) svalue(filterInfoMS$ms1) <- TRUE 
  
  # Disable precursor-related filters when only MS1 are selected
  ms2 <- svalue(filterInfoMS$ms2)
  lapply(filterInfo[c("pmz", "spi", "pc", "mass")], 
         function(x) lapply(x, 
                            function(i) enabled(i) <- ms2))
  
  # Filter spectra
  keep <- mapply(function(data, obj) {
    if(svalue(obj$active)) {
      btwn(data, svalue(obj$from), svalue(obj$to))
    } else rep(TRUE, length(spRtime()))
  }, filterData, filterInfo)
  
  keep <- apply(keep, 1, sum) == length(filterInfo)
  
  if(!svalue(filterInfoMS$ms1)) keep <- keep & spMsLevel() != 1
  if(!svalue(filterInfoMS$ms2)) keep <- keep & spMsLevel() != 2
  
  saveFilterValues()
  
  blockFilters(FALSE)
  
  # Update sequence of spectra and update chart
  if(any(keep)) {
    newSequence <- spIndex()[keep]
    # Compare sequences and update if necessary
    if(!identical(currSequence, newSequence)) {
      
      prevIndex <- currSequence[counter]
      env$currSequence <- newSequence
      
      # Find the nearest spectrum to the one before
      env$counter <- which.min(abs(newSequence - prevIndex))
      
      # also update graphs if index has changed
      if(prevIndex!=newSequence[counter]) updateSpectrum()
      else updateSpectrumInfo()
    } 
  } else cat("\nNot a single spectrum survived filtering!")
}

filterXIC <- function(h, ...) {
  blockFilters()
  validityCheck(filterInfoXIC$XIC, filterData[[3]], filterValuesXIC)
  
  blockFilters(FALSE)
}