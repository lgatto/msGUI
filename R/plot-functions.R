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
    
    svalue(plotBottom) <- filename 
    
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
      
      plotTop <<- tkrplot(getToolkitWidget(groupPlots), 
                     fun = plotSpectrumGraph)
      add(groupPlots, plotTop)
      
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
  dtmax <- max(dt[, 2])
  
  dt[, 2] <- dt[, 2]/dtmax
  
  par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.75)
  
  time <- proc.time()
  
  plot(dt, type = "l", ylim=c(0, 1.075), xlab="Retention time", ylab="Total ion count")
#   abline(v=spRtime(index), col="red", lty=3)
  lines(x=rep(spRtime(index), 2), y=c(0, 1), col="red", lty=3)
  
  par(cex=.75, adj=0)
  text(x=min(dt[, 1]), y=1.075, 
       labels=paste("Max total ions ", pksmax))
  
  time <- proc.time() - time
  
  cat("\nxic:", dim(dt)[1], "data points plotted in", 
      time[3]*1000, "miliseconds")
}

plotSpectrumGraph <- function() {
  
  pks <- peaks(index)
  pksmax <- max(pks[, 2])
  
  mx <- pks[order(pks[, 2], decreasing=TRUE), ][1:5, ]
  pks[, 2] <- pks[, 2]/pksmax
  
  time <- proc.time()
  
  par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
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
  
  time <- proc.time() - time
  
  cat("\nspectrum:", dim(pks)[1], "data points plotted in", time[3]*1000, "miliseconds")
  
}
