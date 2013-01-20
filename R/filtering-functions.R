saveFilterValues <- function(env) {
  env$filterValues <- lapply(env$filterInfo, function(x) lapply(x, function(i) svalue(i)))
  env$filterValuesXIC <- lapply(env$filterInfoXIC, function(x) lapply(x, function(i) svalue(i)))
  env$filterValuesMS <- lapply(env$filterInfoMS, svalue)
}

blockFilters <- function(block=TRUE, env) {
  if(block) {
    mapply(function(obj, id) mapply(blockHandler, obj, id), env$filterInfo, env$filterSpectraHandlerIDs)
    mapply(function(obj, id) mapply(blockHandler, obj, id), env$filterInfoXIC, env$filterXICHandlerIDs)
    mapply(blockHandler, env$filterInfoMS, env$filterSpectraMSHandlerIDs)    
  } else {    
    mapply(function(obj, id) mapply(unblockHandler, obj, id), env$filterInfo, env$filterSpectraHandlerIDs)
    mapply(function(obj, id) mapply(unblockHandler, obj, id), env$filterInfoXIC, env$filterXICHandlerIDs)
    mapply(unblockHandler, env$filterInfoMS, env$filterSpectraMSHandlerIDs)   
  }
}

updateRanges <- function(obj, values, transform, digits, format) {
  svalue(obj$active) <- FALSE
  svalue(obj$from) <- transform(min(values), format, digits=digits)
  svalue(obj$to) <- transform(max(values), format, digits=digits)
}

filterStats <- function(env) {
  cat("Filtering summary\n")
  cat(format("Number of spectra", width=18), length(env$currSequence), "\n", sep="")
  cat(format("MS levels", width=18), paste(unique(env$spMsLevel(env$currSequence)), 
                                           collapse=", "), "\n", sep="")
  mapply(function(name, value) cat(format(name, width=18), 
                                   paste(range(value[env$currSequence]), 
                                         collapse=" - "), "\n", sep=""), 
         env$filterNames, env$filterData)
}
  
filterReset <- function(env) {  
  blockFilters(TRUE, env)
  
  # Filters for spectra
  mapply(updateRanges, env$filterInfo, env$filterData, env$filterTransform, 
         env$settings$digits, env$settings$RtFormat)
  
  # MS levels checkboxes
  lapply(env$filterInfoMS, function(i) svalue(i) <- TRUE)
  
  # Filter for XIC
  svalue(env$filterInfoXIC$XIC$x) <- round(mean(range(env$filterData[[3]]), digits=0))
    
  saveFilterValues(env)  
  blockFilters(FALSE, env)
  
}

filterSwitch <- function(on, env) {
  lapply(env$filterInfo, function(x) lapply(x, function(i) enabled(i) <- on))
  lapply(env$filterInfoXIC, function(x) lapply(x, function(i) enabled(i) <- on))
  lapply(env$filterInfoMS, function(x) enabled(x) <- on)
  
  if(!env$anyMS1spectra) {
    svalue(env$filterInfoMS[[1]]) <- FALSE
    svalue(env$filterInfoMS[[2]]) <- TRUE
    lapply(env$filterInfoMS, function(x) enabled(x) <- FALSE)
    lapply(env$filterInfoXIC, function(x) lapply(x, function(i) enabled(i) <- FALSE))
  }
}

btwn <- function(x, from, to) x >= from & x <= to

validityCheck <- function(object, data, pastValues, detransform, transform, format, digits) {

  range <- range(data)
  
  numFrom <- detransform(svalue(object$from), format)
  numTo <- detransform(svalue(object$to), format)
  
  prevFrom <- detransform(pastValues$from, format)
  prevTo <- detransform(pastValues$to, format)
    
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
  svalue(object$from) <- transform(numFrom, format, digits)
  svalue(object$to) <- transform(numTo, format, digits)
}

filterSpectra <- function(h, ...) {
  env <- h$action
  blockFilters(TRUE, env)
  
  # Check entered values, fix if needed
  mapply(validityCheck, env$filterInfo, env$filterData, env$filterValues, 
         env$filterDetransform, env$filterTransform, env$settings$RtFormat, env$settings$digits)
  
  pmzRange <- range(env$filterData[[3]])
  if (!btwn(svalue(env$filterInfoXIC$XIC$x), pmzRange[1], pmzRange[2]))
    svalue(env$filterInfoXIC$XIC$x) <- round(mean(pmzRange), digits=0)
  
  # If XIC filter active, its values are forced upon the Precursor MZ filter
  if(env$anyMS1spectra & svalue(env$filterInfoXIC$XIC$active)) {
    svalue(env$filterInfo$pmz$active) <- TRUE
    svalue(env$filterInfo$pmz$from) <- max(pmzRange[1], svalue(env$filterInfoXIC$XIC$x) - env$settings$Da/2)
    svalue(env$filterInfo$pmz$to) <- min(pmzRange[2], svalue(env$filterInfoXIC$XIC$x) + env$settings$Da/2)
  }
  
  updateXIC <- FALSE
  if(svalue(env$filterInfoXIC$XIC$active)!=env$filterValuesXIC$XIC$active & env$anyMS1spectra)
    updateXIC <- TRUE
  
  # XIC filtering
  if(env$anyMS1spectra) {
    
    if(svalue(env$filterInfoXIC$XIC$active)) {
      
      from <- max(pmzRange[1], svalue(env$filterInfoXIC$XIC$x) - env$settings$Da)
      to <- min(pmzRange[2], svalue(env$filterInfoXIC$XIC$x) + env$settings$Da)
      # This looks at each group of MS2 spectra and returns TRUE if any spectrum 
      # in that group survived filtering
      filtered <- apply(env$MS2indices, 1, function(x) any(btwn(env$filterData[[3]][x[1]:x[2]], 
                                                            from, to))) 
      # Now get the indices of MS1 spectra that are just before these MS2 spectra
      MS1indices <- env$MS2indices[filtered, 1] - 1
      # Recalculate these indices for data.frame with MS1 spectra only 
      
      MS1indicesL <- rep(FALSE, env$nSpectra)
      MS1indicesL[MS1indices] <- TRUE
      
      env$XICvalues <- MS1indicesL[env$spMsLevel()==1]
      updateXIC <- TRUE
    } else env$XICvalues <- TRUE
  }
  
  # Make sure at least one MS level is selected
  if(!svalue(env$filterInfoMS$ms1) & !env$filterValuesMS$ms2) svalue(env$filterInfoMS$ms2) <- TRUE
  if(!svalue(env$filterInfoMS$ms2) & !env$filterValuesMS$ms1) svalue(env$filterInfoMS$ms1) <- TRUE 
  
  # Disable precursor-related filters when only MS1 are selected
  ms2 <- svalue(env$filterInfoMS$ms2)
  lapply(env$filterInfo[c("pmz", "spi", "pc", "mass")], 
         function(x) lapply(x, 
                            function(i) enabled(i) <- ms2))
  
  # Filter spectra
  keep <- mapply(function(data, obj, detransform) {
    if(svalue(obj$active)) {
      btwn(data, detransform(svalue(obj$from)), detransform(svalue(obj$to)))
    } else rep(TRUE, env$nSpectra)
  }, env$filterData, env$filterInfo, env$filterDetransform)
  
  colnames(keep) <- names(env$filterInfo)
  
  MS1filters <- c("rt", "index")
  
  saveFilterValues(env)
  
  keepMS1 <- apply(keep[, MS1filters], 1, sum) == length(MS1filters)
  keepMS2 <- apply(keep, 1, sum) == length(env$filterInfo)
  
  keep <- (keepMS1 | env$spMsLevel() != 1) & (keepMS2 | env$spMsLevel() != 2)
  
  if(!svalue(env$filterInfoMS$ms1)) keep <- keep & env$spMsLevel() != 1
  if(!svalue(env$filterInfoMS$ms2)) keep <- keep & env$spMsLevel() != 2
  
  blockFilters(FALSE, env)
  
  # Update sequence of spectra and update chart
  if(any(keep)) {
    newSequence <- env$spIndex()[keep]
    # Compare sequences and update if necessary
    if(!identical(env$currSequence, newSequence) | env$forceRedraw) {

      prevIndex <- env$currSequence[env$counter]
      env$currSequence <- newSequence
      
      # Find the nearest spectrum to the one before
      env$counter <- which.min(abs(newSequence - prevIndex))
      
      # also update graphs if index has changed
      if(prevIndex!=newSequence[env$counter] | env$forceRedraw) {
        updateXIC <- FALSE
        updateSpectrum(h=list(action=list(0, env)))
        env$forceRedraw <- FALSE
      }
      else updateSpectrumInfo(FALSE, env)
      buttonSwitch(TRUE, env)
    } 
    
    if (env$verbose) filterStats(env)
    
  } else {
    
    visible(env$plotTop) <- TRUE
    plotMsg()
    visible(env$plotBottom) <- TRUE
    plotMsg()
    
    if(!env$zoomWindowClosed) {
      visible(env$plotZoom) <- TRUE      
      plotMsg()
    }
    
    if(!env$XICWindowClosed) {
      visible(env$plotXICw) <- TRUE      
      plotMsg()
    }
    
    visible(env$msGUIWindow) <- TRUE
    
    buttonSwitch(FALSE, env)
    clickSwitch(FALSE, env)
    
    env$forceRedraw <- TRUE
    updateSpectrumInfo(TRUE, env)
    
    gmessage(message="No spectra survived filtering!\nPlease select different filter values.", 
             title="")
  }
  if(updateXIC) {
    plotXIC(env$XICZoom, noCache=TRUE, env) 
    if(!env$zoomWindowClosed) {
      visible(env$plotZoom) <- TRUE      
      plotSpectrumZoom(env$spectrumZoom, env)
    }
    if(!env$XICWindowClosed) {
      visible(env$plotXICw) <- TRUE      
      plotChromaZoom(env)
    }
  }
}