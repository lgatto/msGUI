plotXIC <- function() {
  
  if(device=="png") {
    
    filename <- tempfile(fileext=".png")
    png(filename, width=500, height=250)
    
    plotChromatogram()
    
    dev.off()
    
    visible(plotBottom) <- TRUE 
    grid::grid.raster(readPNG(filename))
    
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
    
    svalue(plotTop) <- filename 
    
  }
  
}

plotSpectrum <- function() {
  
  if(device=="png") {
    
    filename <- tempfile(fileext=".png")
    png(filename, width=500, height=250)
    
    plotSpectrumGraph()
    
    dev.off()
    
    visible(plotTop) <- TRUE 
    grid::grid.raster(readPNG(filename))
    
  } else if(device=="cairo") {
    
    visible(plotTop) <- TRUE 
    plotSpectrumGraph()    
    
  } else if(device=="tkrplot") {
    
    if(plotTopDrawn) {
      
      tkrreplot(plotTop)
      
    } else {
      
      plotBottom <<- tkrplot(getToolkitWidget(groupPlots), 
                     fun = plotSpectrumGraph)
      add(groupPlots, plotBottom)
      
      plotTopDrawn <<- TRUE
      
    }
    
  } else if(device=="gimage") {    
    
    filename <- tempfile(fileext=".png")
    png(filename, width=500, height=250)
    
    plotSpectrumGraph()
    
    dev.off()
    
    svalue(plotTop) <- filename 
    
  }
  
}

plotChromatogram <- function() {
  
  dt <- xic(n=c(1, 2)[c(svalue(filterInfo$ms1), svalue(filterInfo$ms2))], 
            FALSE)
  par(mar=c(3,3,.5,1), mgp=c(2,.7,0), tck=-.01, bty="n", lab=c(5, 3, 7), las=0)
  
  time <- system.time({plot(dt, type = "l"); 
                       abline(v=spRtime(index), col="red", lty=3)})
  
  cat("\nxic:", dim(dt)[1], "data points plotted in", 
      time[3]*1000, "miliseconds")
}

plotSpectrumGraph <- function() {
  
  pks <- peaks(index)
  pksmax <- max(pks[, 2])
  
  mx <- pks[order(pks[, 2], decreasing=TRUE), ][1:5, ]
  pks[, 2] <- pks[, 2]/pksmax
  
  time <- proc.time()
  par(mar=c(3,3,1,1), mgp=c(2,.7,0), tck=-.01, bty="n", lab=c(5, 3, 7), las=1)
  plot(pks, xlab="Intensity", ylab="M/Z", 
       type = ifelse(spMsLevel(index)==1, "l", "h"), 
       xlim=c(min(spLowMZ()[spMsLevel()==spMsLevel(index)]), 
              max(spHighMZ()[spMsLevel()==spMsLevel(index)]))) 
  text(x=mx[, 1], y=mx[, 2]/pksmax, labels=round(mx[, 1], digits=3), col="grey50", adj=c(0, 0))
  text(x=0, y=1.05, labels=mx[1, 2])
  abline(h=0, col="grey50")
  time <- proc.time() - time
  
  cat("\nspectrum:", dim(pks)[1], "data points plotted in", time[3]*1000, "miliseconds")
  
}
