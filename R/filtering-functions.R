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
  
  if(!env$anyMS1spectra) {
    svalue(filterInfoMS[[1]]) <- FALSE
    svalue(filterInfoMS[[2]]) <- TRUE
    lapply(filterInfoMS, function(x) enabled(x) <- FALSE)
    lapply(filterInfoXIC, function(x) lapply(x, function(i) enabled(i) <- FALSE))
  }
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
  mapply(validityCheck, filterInfoXIC, filterData[3], filterValuesXIC)
  
  # If XIC filter active, its values are forced upon the Precursor MZ filter
  if(env$anyMS1spectra & svalue(filterInfoXIC$XIC$active)) {
    svalue(filterInfo$pmz$active) <- TRUE
    svalue(filterInfo$pmz$from) <- svalue(filterInfoXIC$XIC$from)
    svalue(filterInfo$pmz$to) <- svalue(filterInfoXIC$XIC$to)
  }
  
  updateXIC <- FALSE
  if(svalue(filterInfoXIC$XIC$active)!=filterValuesXIC$XIC$active & anyMS1spectra)
    updateXIC <- TRUE
  
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
    } else rep(TRUE, nSpectra)
  }, filterData, filterInfo)
  
  colnames(keep) <- names(filterInfo)
  
  MS1filters <- c("rt", "index")
  
  keepMS1 <- apply(keep[, MS1filters], 1, sum) == length(MS1filters)
  keepMS2 <- apply(keep, 1, sum) == length(filterInfo)
  
  keep <- (keepMS1 | spMsLevel() != 1) & (keepMS2 | spMsLevel() != 2)
  
  if(!svalue(filterInfoMS$ms1)) keep <- keep & spMsLevel() != 1
  if(!svalue(filterInfoMS$ms2)) keep <- keep & spMsLevel() != 2
  
  saveFilterValues()
  
  blockFilters(FALSE)
  
  # XIC filtering
  if(env$anyMS1spectra) {
    
    if(svalue(filterInfoXIC$XIC$active)) {
      
      from <- svalue(filterInfoXIC$XIC$from)
      to <- svalue(filterInfoXIC$XIC$to)
      # This looks at each group of MS2 spectra and returns TRUE if any spectrum 
      # in that group survived filtering
      filtered <- apply(MS2indices, 1, function(x) any(btwn(filterData[[3]][x[1]:x[2]], 
                                                            from, to))) 
      # Now get the indices of MS1 spectra that are just before these MS2 spectra
      MS1indices <- MS2indices[filtered, 1] - 1
      # Recalculate these indices for data.frame with MS1 spectra only 
      
      MS1indicesL <- rep(FALSE, nSpectra)
      MS1indicesL[MS1indices] <- TRUE
      
      env$XICvalues <- MS1indicesL[spMsLevel()==1]
      updateXIC <- TRUE
    } else env$XICvalues <- TRUE
  }
  
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
      if(prevIndex!=newSequence[counter]) {
        updateXIC <- FALSE
        updateSpectrum()
      }
      else updateSpectrumInfo()
    } 
  } else cat("Not a single spectrum survived filtering!\n", 
             "Please select different values or turn off some filters.\n", sep="")
  if(updateXIC) {
    plotXIC(XICZoom, noCache=TRUE) 
    if(!zoomWindowClosed) {
      visible(env$plotZoom) <- TRUE      
      plotSpectrumZoom(spectrumZoom)
    }
  }
}

filterXIC <- function(h, ...) {
  blockFilters()
  validityCheck(filterInfoXIC$XIC, filterData[[3]], filterValuesXIC)
  
  blockFilters(FALSE)
}