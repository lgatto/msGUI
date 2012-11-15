saveFilterValues <- function() {
  env$filterValues <- lapply(filterInfo, function(x) lapply(x, function(i) svalue(i)))
  env$filterValuesMS <- lapply(filterInfoMS, svalue)
}

blockFilters <- function(block=TRUE) {
  if(block) {
    mapply(function(obj, id) mapply(blockHandler, obj, id), filterInfo, filterSpectraHandlerIDs)
    mapply(blockHandler, filterInfoMS, filterSpectraMSHandlerIDs)    
  } else {    
    mapply(function(obj, id) mapply(unblockHandler, obj, id), filterInfo, filterSpectraHandlerIDs)
    mapply(unblockHandler, filterInfoMS, filterSpectraMSHandlerIDs)   
  }
}
  
filterReset <- function(env) {
  
  filterData <- list(spRtime(), spIndex(), spPrecMz(), spPrecInt(), spPrecCharge(), spPrecMz()*spPrecCharge())
    
  blockFilters()
  
  mapply(function(obj, values) {
    svalue(obj$active) <- FALSE
    svalue(obj$from) <- min(values)
    svalue(obj$to) <- max(values)
  }, filterInfo, filterData)
  
  lapply(filterInfoMS, function(i) svalue(i) <- TRUE)
  
  saveFilterValues()
  
  blockFilters(FALSE)
}

filterSwitch <- function(on) {
  lapply(filterInfo, function(x) lapply(x, function(i) enabled(i) <- on))
}

btwn <- function(x, from, to) {
  stopifnot(all(is.numeric(x), is.character(from), is.character(to)))
  x >= as.numeric(from) & x <= as.numeric(to)
}


validityCheck <- function(object, data, pastValues) {
  
  from <- object$from
  to <- object$to
  
  numFrom <- as.numeric(svalue(from))
  numTo <- as.numeric(svalue(to))
  
  prevFrom <- as.numeric(pastValues$from)
  prevTo <- as.numeric(pastValues$to)
  
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
  
  filterData <- list(spRtime(), spIndex(), spPrecMz(), spPrecInt(), spPrecCharge(), spPrecMz()*spPrecCharge())
  
  blockFilters()
  
  # Check current values, make fixes if necessary
  mapply(validityCheck, filterInfo, filterData, filterValues)
  
  # If MS1 was deselected previously, and user deselects MS2, 
  # select MS1 in order to have something selected
  if(!svalue(filterInfoMS$ms1) & !filterValuesMS$ms2) svalue(filterInfoMS$ms2) <- TRUE    
  # And vice versa
  if(!svalue(filterInfoMS$ms2) & !filterValuesMS$ms1) svalue(filterInfoMS$ms1) <- TRUE 
  
  # Disable Precursor MZ filter if only MS2 selected  
  lapply(filterInfo$rt, function(i) enabled(i) <- svalue(filterInfoMS$ms2))
  
  # Filter!
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
    } 
  } else cat("\nNot a single spectrum survived filtering!")
}
