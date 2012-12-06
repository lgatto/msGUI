plotMsg <- function(message) {
  par(mar=rep(0, 4))
  plot(x=c(0, 1), y=c(0, 1), col="white", , bg="grey10", axes=FALSE)
  text(x=.5, y=.5, labels=message)
}

plotXIC <- function() {
  
  if(device=="png") {
    if(is.null(cache$xic[[currSequence[counter]]])) {
      
      filename <- tempfile(fileext=".png")    
      png(filename, width=500, height=250, restoreConsole=FALSE)    
      plotChromatogram()    
      dev.off()    
      env$cache$xic[[currSequence[counter]]] <- filename
      
    } else filename <- cache$xic[[currSequence[counter]]]
    
    visible(plotBottom) <- TRUE 
    plotMsg("Refreshing...")  
    lim <- par()$usr
    rasterImage(readPNG(filename), lim[1], lim[3], lim[2], lim[4], interpolate=FALSE)
    
  } else if(device=="cairo") {
    
    visible(plotBottom) <- TRUE 
    plotChromatogram()    
    
  } else if(device=="tkrplot") {
    
    if(plotBottomDrawn) {
      
      tkrreplot(plotBottom)
      
    } else {
      
      plotBottom <<- tkrplot(getToolkitWidget(groupPlots), 
                          fun = plotChromatogram)
      add(groupPlots, plotBottom)
      
      plotBottomDrawn <<- TRUE
      
    }
    
  } else if(device=="gimage") {    
    
    filename <- tempfile(fileext=".png")
    png(filename, width=500, height=250)
    
    plotChromatogram()
    
    dev.off()
    
    svalue(plotBottom) <- filename 
    
  }
  
}

plotSpectrum <- function(zoom=NULL) {
  
  if(device=="png") {
    if(is.null(cache$spectra[[currSequence[counter]]]) | !is.null(zoom)) {
      
      if(is.null(zoom) & verbose) cat("\ncaching", currSequence[counter])
      
      filename <- tempfile(fileext=".png")
      
      png(filename, width=500, height=250, restoreConsole=FALSE)
      
      plotSpectrumGraph(zoom=zoom)
      
      dev.off()
      
      if(is.null(zoom)) env$cache$spectra[[currSequence[counter]]] <- filename
      
    } else {
      if(verbose) cat("\nloading", currSequence[counter], "from cache")
      filename <- cache$spectra[[currSequence[counter]]]
    }
        
    visible(plotTop) <- TRUE 
#     grid::grid.raster(readPNG(filename), width=1, height=1) # Not specifying 
    # whidth and height results in handlerChanged retrieving NaN and Inf 
    # coordinates. 
    
    # Correction: specifying these does not help after all. The reason for 
    # occasional bad coordinates is not clear. 
    
    # Another method without package:grob. Also works but graphs 
    # flicker when reloading which is not super amazing.
    
    plotMsg("Refreshing...")
    
    lim <- par()$usr
    rasterImage(readPNG(filename), lim[1], lim[3], lim[2], lim[4], interpolate=FALSE)
        
  } else if(device=="cairo") {
    
    visible(plotTop) <- TRUE 
    plotSpectrumGraph(zoom=zoom)    
    
  } else if(device=="tkrplot") {
    
    if(plotTopDrawn) {
      
      tkrreplot(plotTop)
      
    } else {
      
      plotTop <<- tkrplot(getToolkitWidget(groupPlots), 
                     fun = plotSpectrumGraph)
      add(groupPlots, plotTop)
      
      plotTopDrawn <<- TRUE
      
    }
    
  } else if(device=="gimage") {    
    
    filename <- tempfile(fileext=".png")
    png(filename, width=500, height=250)
    
    plotSpectrumGraph(zoom=zoom)
    
    dev.off()
    
    svalue(plotTop) <- filename 
    
  }
  
}

plotChromatogram <- function() {
  
  dt <- xic(n=1, FALSE)
  dtmax <- max(dt[, 2])
  
  dt[, 2] <- dt[, 2]/dtmax
  
  par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.75)
  
  time <- proc.time()
  .GlobalEnv$.msgui <- FALSE
  
  plot(dt, type = "l", ylim=c(0, 1.075), xlab="Retention time", ylab="Total ion count")
#   abline(v=spRtime(index), col="red", lty=3)
  # abline went too high. 
  lines(x=rep(spRtime(currSequence[counter]), 2), y=c(0, 1), col="red", lty=3)
  
  par(cex=.75, adj=0)
  text(x=min(dt[, 1]), y=1.075, 
       labels=paste("Max total ions ", dtmax))
  
  par(parSave)
  time <- proc.time() - time
  
  .GlobalEnv$.msgui <- TRUE
  
  if(verbose) cat("\nxic:", dim(dt)[1], "data points plotted in", 
      time[3]*1000, "miliseconds")
}

plotSpectrumGraph <- function(zoom=NULL) {
  
  pks <- peaks(currSequence[counter])
  pksmax <- max(pks[, 2])
  
  mx <- pks[order(pks[, 2], decreasing=TRUE), ][1:5, ]
  pks[, 2] <- pks[, 2]/pksmax
  
  time <- proc.time()
  .GlobalEnv$.msgui <- FALSE
  
#   browser()
  
  env$parSave <- par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.75)
  plot(pks, xlab="Mass to charge ratio (M/Z)", ylab="Intensity", 
       type = ifelse(spMsLevel(index)==1, "l", "h"), 
       xlim=c(min(spLowMZ()[spMsLevel()==spMsLevel(index)]), 
              max(spHighMZ()[spMsLevel()==spMsLevel(index)])), 
       ylim=c(0, 1.075)) 
  text(x=mx[, 1], y=mx[, 2]/pksmax, labels=round(mx[, 1], digits=3), 
       col="grey50", adj=c(0, 0))
  par(cex=.75, adj=0)
  text(x=min(spLowMZ()[spMsLevel()==spMsLevel(index)]), y=1.075, 
       labels=paste("Base peak intensity ", pksmax))
  abline(h=0, col="grey50")
  if(!is.null(zoom)) {
    lines(zoom$x[c(1, 1, 2, 2, 1)], zoom$y[c(1, 2, 2, 1, 1)], 
          col="grey25", lty=3)
  }
  
  par(parSave)
  time <- proc.time() - time
  
  .GlobalEnv$.msgui <- TRUE
  
  if(verbose) cat("\nspectrum:", dim(pks)[1], "data points plotted in", time[3]*1000, "miliseconds")
  
}

plotSpectrumZoom <- function(limits=NULL) {
  if(verbose) cat("\nlimits:", limits$x, limits$y)
  
  .GlobalEnv$.msgui <- FALSE
  pks <- peaks(index)
  pksmax <- max(pks[, 2])
  
# attempt to make nice labels in the visible area of zoom chart
#   mx <- pks[pks[, 1] > limits$x[1] | pks[, 1] < limits$x[2] | pks[, 2] > limits$y[1] | pks[, 2] < limits$y[2], ]
#   mx <- mx[order(mx[, 2], decreasing=TRUE), ][1:5, ]
  pks[, 2] <- pks[, 2]/pksmax
  
#   par(mar=rep(.1, 4), mgp=c(3,1,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
#       adj=.5, las=1, cex=0.5)
#   par(parSave)
  par(mar=c(3, 3, 0, 0), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.65)
  plot(pks, xlab="Mass to charge ratio (M/Z)", ylab="Intensity", #zero.line=TRUE, 
       type = ifelse(spMsLevel(index)==1, "l", "h"), 
       xlim=limits$x, ylim=limits$y) 
#   text(x=mx[, 1], y=mx[, 2]/pksmax, labels=round(mx[, 1], digits=3), 
#        col="grey50", adj=c(0, 0))
  
  .GlobalEnv$.msgui <- TRUE
  
}
