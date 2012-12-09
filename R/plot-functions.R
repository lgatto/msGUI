plotMsg <- function(message="Nothing to display", col="grey50") {
  par(mar=rep(0, 4))
  plot(0:1, 0:1, col="transparent", axes=FALSE)
  text(.5, .5, labels=message, cex=.65, col=col)
}

plotXIC <- function(zoom=NULL) 
  plotGeneric(zoom=zoom, cacheName="xic", 
              height=settings$chromaHeight, width=settings$width, 
              fun=plotChromatogram, plotObject=plotBottom)

plotSpectrum <- function(zoom=NULL) 
  plotGeneric(zoom=zoom)

plotGeneric <- function(zoom=NULL, cacheName="spectra", 
                        height=settings$spectrumHeight, width=settings$width, 
                        fun=plotSpectrumGraph, 
                        plotObject=plotTop) {
  
  if(device=="png") {
    
    if(is.null(cache[[cacheName]][[index]]) | !is.null(zoom)) {
      if(is.null(zoom) & verbose) cat("\ncaching", index)
      
      filename <- tempfile(fileext=".png") 
      if(.Platform$OS.type=="windows") 
        png(filename, width=settings$width, height=settings$spectrumHeight, restoreConsole=FALSE)
      else 
        png(filename, width=500, height=250)  
      fun(zoom=zoom)    
      dev.off()    
      if(is.null(zoom)) env$cache[[cacheName]][[index]] <- filename
    } else filename <- cache[[cacheName]][[index]]
    
    visible(plotObject) <- TRUE
    lim <- par()$usr
    rasterImage(readPNG(filename), lim[1], lim[3], lim[2], lim[4], interpolate=FALSE)
    
  } else if(device=="cairo") {
    
    visible(plotObject) <- TRUE 
    fun(zoom=zoom)    
    
  } 
  
}

plotChromatogram <- function(zoom=NULL) {
  
  xic <- xic(n=1, FALSE)
  maxTotalIons <- max(xic[, 2])  
  xic[, 2] <- xic[, 2]/maxTotalIons
  par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.75)
  
  time <- proc.time()
  
  plot(xic, type = "l", ylim=c(0, 1.075), xlab="Retention time", ylab="Total ion count")
  lines(x=rep(spRtime(index), 2), y=c(0, 1), col="red", lty=3)
  text(x=min(xic[, 1]), y=1.075, adj=0, 
       labels=paste("Max total ions ", maxTotalIons))
  
  if(!is.null(zoom)) {
    lines(zoom$x[c(1, 1, 2, 2, 1)], zoom$y[c(1, 2, 2, 1, 1)], 
          col="grey25", lty=3)
  }
  
  if(verbose) {
    time <- proc.time() - time
    cat("\nxic:", dim(xic)[1], "data points plotted in", 
      time[3]*1000, "miliseconds")
  }
}

plotSpectrumGraph <- function(zoom=NULL) {

  MsLevel <- spMsLevel(index)
  peaks <- peaks(index)
  basePeak <- max(peaks[, 2])
  peaks[, 2] <- peaks[, 2]/basePeak  
  labels <- peaks[order(peaks[, 2], decreasing=TRUE), ]  
  xCoords <- labels[, 1]  
  current <- 1  
  out <- 1
  
  while(length(out) < settings$labelNumber) {
    s <- 1
    while(min(abs(xCoords[out] - xCoords[current+s])) <= settings$labelThreshold) s <- s + 1
    out <- c(out, current <- current + s)
  }
  
  labels <- labels[out, ]
  
  time <- proc.time()
  
  par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.75)
  plot(peaks, xlab="Mass to charge ratio (M/Z)", ylab="Intensity", 
       type = ifelse(MsLevel==1, settings$MS1PlotType, settings$MS2PlotType), 
       xlim=xLimits[, MsLevel], ylim=c(0, 1.075)) 
  text(labels, labels=round(labels[, 1], settings$digits), 
       col="grey50", adj=c(0, 0))
  par(adj=0)   
  text(x=min(xLimits[, MsLevel]), y=1.075, 
       labels=paste("Base peak intensity ", round(basePeak, settings$digits)))
  
  abline(h=0, col="grey50")
  
  if(!is.null(zoom)) {
    lines(zoom$x[c(1, 1, 2, 2, 1)], zoom$y[c(1, 2, 2, 1, 1)], 
          col="grey25", lty=3)
  }
  
  time <- proc.time() - time
  
  if(verbose) cat("\nspectrum:", dim(peaks)[1], "data points plotted in", time[3]*1000, "miliseconds")
  
}

plotSpectrumZoom <- function(limits=NULL) {
  
  if(verbose) cat("\nlimits:", limits$x, limits$y)
  
  peaks <- peaks(index)
  peaks[, 2] <- peaks[, 2]/max(peaks[, 2])  
  
  inFocus <- peaks[, 1] > limits$x[1] & peaks[, 1] < limits$x[2]
  inFocus <- inFocus & peaks[, 2] > limits$y[1] & peaks[, 2] < limits$y[2]
  numberInFocus <- sum(inFocus)
  
  par(mar=c(3, 3, 0, 0), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.65)
  plot(peaks, xlab="Mass to charge ratio (M/Z)", ylab="Intensity", 
       type = ifelse(spMsLevel(index)==1, settings$MS1PlotType, settings$MS2PlotType), 
       xlim=limits$x, ylim=limits$y) 
  if(numberInFocus>0) {    
    labels <- peaks[inFocus, , drop=FALSE]
    labels <- labels[order(labels[, 2], decreasing=TRUE), , drop=FALSE]
    labels <- labels[1:min(numberInFocus, settings$labelNumber), , drop=FALSE]
    text(labels, labels=round(labels[, 1], digits=settings$digits), 
         col="grey50", adj=c(0, 0))    
  }
  par(adj=1)
  text(limits$x[2], limits$y[2], labels=paste("Acquisition number:", spIndex(index)))
}

plotChromaZoom <- function() {
  
  xic <- xic(n=1, FALSE)
  maxTotalIons <- max(xic[, 2])
  xic[, 2] <- xic[, 2]/maxTotalIons
  
  par(mar=c(3, 3, 0, 0), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.65)
  plot(xic, type = "l", xlim=env$XICZoom$x, ylim=env$XICZoom$y, 
       xlab="Retention time", ylab="Total ion count")
  
}