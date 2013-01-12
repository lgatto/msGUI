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

updateRanges <- function(obj, values, transform) {
  svalue(obj$active) <- FALSE
  svalue(obj$from) <- transform(round(min(values), digits=settings$digits))
  svalue(obj$to) <- transform(round(max(values), digits=settings$digits))
}

filterStats <- function(names, data, mslevels, sequence) {
  cat("Filtering summary\n")
  cat(format("Number of spectra", width=18), length(currSequence), "\n", sep="")
  cat(format("MS levels", width=18), paste(unique(spMsLevel(currSequence)), collapse=", "), "\n", sep="")
  mapply(function(name, value) cat(format(name, width=18), paste(range(value[currSequence]), collapse=" - "), "\n", sep=""), 
         filterNames, filterData)
}
  
filterReset <- function(env) {  
  blockFilters()
  
  # Filters for spectra
  mapply(updateRanges, filterInfo, filterData, filterTransform)
  
  # MS levels checkboxes
  lapply(filterInfoMS, function(i) svalue(i) <- TRUE)
  
  # Filter for XIC
  mapply(updateRanges, filterInfoXIC, list(spPrecMz()), list(identity))
    
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

btwn <- function(x, from, to) x >= from & x <= to

validityCheck <- function(object, data, pastValues, detransform, transform) {

  range <- range(data)
  
  numFrom <- detransform(svalue(object$from))
  numTo <- detransform(svalue(object$to))
  
  prevFrom <- detransform(pastValues$from)
  prevTo <- detransform(pastValues$to)
    
  prevValid <- prevFrom >= range[1] & prevFrom < prevTo & prevTo > range[1] & prevTo <= range[2]
  if(is.na(numFrom)) numFrom <- range[1] else {
    if(numFrom < range[1]) numFrom <- range[1]
    if(numFrom > range[2]) numFrom <- ifelse(prevValid, prevFrom, range[1])
  }
  if(is.na(numTo)) numTo <- ifelse(prevValid, prevTo, range[2]) else {
    if(numTo < range[1]) numTo <- ifelse(prevValid, prevTo, range[2])
    if(numTo > range[2]) numTo <- range[2]
  }
  if(numFrom > numTo) {
    if(numFrom!=prevFrom) numFrom <- prevFrom
    else numTo <- prevTo
  }
  svalue(object$from) <- transform(numFrom)
  svalue(object$to) <- transform(numTo)
}

filterSpectra <- function(h, ...) {
  blockFilters()
  
  # Check entered values, fix if needed
  mapply(validityCheck, filterInfo, filterData, filterValues, 
         filterDetransform, filterTransform)
  mapply(validityCheck, filterInfoXIC, filterData[3], filterValuesXIC, 
         list(as.numeric), list(identity))
  
  # If XIC filter active, its values are forced upon the Precursor MZ filter
  if(env$anyMS1spectra & svalue(filterInfoXIC$XIC$active)) {
    svalue(filterInfo$pmz$active) <- TRUE
    svalue(filterInfo$pmz$from) <- svalue(filterInfoXIC$XIC$from)
    svalue(filterInfo$pmz$to) <- svalue(filterInfoXIC$XIC$to)
  }
  
  updateXIC <- FALSE
  if(svalue(filterInfoXIC$XIC$active)!=filterValuesXIC$XIC$active & anyMS1spectra)
    updateXIC <- TRUE
  
  # XIC filtering
  if(env$anyMS1spectra) {
    
    if(svalue(filterInfoXIC$XIC$active)) {
      
      from <- as.numeric(svalue(filterInfoXIC$XIC$from))
      to <- as.numeric(svalue(filterInfoXIC$XIC$to))
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
  
  # Make sure at least one MS level is selected
  if(!svalue(filterInfoMS$ms1) & !filterValuesMS$ms2) svalue(filterInfoMS$ms2) <- TRUE
  if(!svalue(filterInfoMS$ms2) & !filterValuesMS$ms1) svalue(filterInfoMS$ms1) <- TRUE 
  
  # Disable precursor-related filters when only MS1 are selected
  ms2 <- svalue(filterInfoMS$ms2)
  lapply(filterInfo[c("pmz", "spi", "pc", "mass")], 
         function(x) lapply(x, 
                            function(i) enabled(i) <- ms2))
  
  # Filter spectra
  keep <- mapply(function(data, obj, detransform) {
    if(svalue(obj$active)) {
      btwn(data, detransform(svalue(obj$from)), detransform(svalue(obj$to)))
    } else rep(TRUE, nSpectra)
  }, filterData, filterInfo, filterDetransform)
  
  colnames(keep) <- names(filterInfo)
  
  MS1filters <- c("rt", "index")
  
  saveFilterValues()
  
  keepMS1 <- apply(keep[, MS1filters], 1, sum) == length(MS1filters)
  keepMS2 <- apply(keep, 1, sum) == length(filterInfo)
  
  keep <- (keepMS1 | spMsLevel() != 1) & (keepMS2 | spMsLevel() != 2)
  
  if(!svalue(filterInfoMS$ms1)) keep <- keep & spMsLevel() != 1
  if(!svalue(filterInfoMS$ms2)) keep <- keep & spMsLevel() != 2
  
  blockFilters(FALSE)
  
  # Update sequence of spectra and update chart
  if(any(keep)) {
    newSequence <- spIndex()[keep]
    # Compare sequences and update if necessary
    if(!identical(currSequence, newSequence) | forceRedraw) {

      prevIndex <- currSequence[counter]
      env$currSequence <- newSequence
      
      # Find the nearest spectrum to the one before
      env$counter <- which.min(abs(newSequence - prevIndex))
      
      # also update graphs if index has changed
      if(prevIndex!=newSequence[counter] | forceRedraw) {
        updateXIC <- FALSE
        updateSpectrum()
        env$forceRedraw <- FALSE
      }
      else updateSpectrumInfo()
    } 
    
    if (verbose) filterStats()
    
  } else {
    
    visible(plotTop) <- TRUE
    plotMsg()
    visible(plotBottom) <- TRUE
    plotMsg()
    
    if(!zoomWindowClosed) {
      visible(plotZoom) <- TRUE      
      plotMsg()
    }
    
    if(!XICWindowClosed) {
      visible(plotXICw) <- TRUE      
      plotMsg()
    }
    
    visible(msGUIWindow) <- TRUE
    
    buttonSwitch(FALSE)
    clickSwitch(FALSE)
    
    env$forceRedraw <- TRUE
    updateSpectrumInfo(blank=TRUE)
    
    gmessage(message="No spectra survived filtering!\nPlease select different filter values.", 
             title="")
  }
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