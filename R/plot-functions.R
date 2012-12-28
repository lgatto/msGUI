plotMsg <- function(message="Nothing to display", col="grey50") {
  par(mar=rep(0, 4))
  plot(0:1, 0:1, col="transparent", axes=FALSE)
  text(.5, .5, labels=message, cex=.65, col=col)
}

plotXIC <- function(zoom=NULL, noCache=FALSE)
  plotGeneric(zoom=zoom, cacheName="xic", 
              height=settings$chromaHeight, width=settings$width, 
              fun=plotChromatogram, plotObject=plotBottom, noCache=noCache, 
              name="chromatogram")

plotSpectrum <- function(zoom=NULL, noCache=FALSE) 
  plotGeneric(zoom=zoom, noCache=noCache)

plotGeneric <- function(zoom=NULL, cacheName="spectra", 
                        height=settings$spectrumHeight, width=settings$width, 
                        fun=plotSpectrumGraph, 
                        plotObject=plotTop, noCache, name="spectrum") {
  
  if(device=="png") {
    
    if(is.null(cache[[cacheName]][[index]]) | !is.null(zoom) | noCache) {
      if(is.null(zoom) & verbose) cat("caching", name, index, "\n")
      
      time <- proc.time()
      
      filename <- tempfile(fileext=".png") 
      if(.Platform$OS.type=="windows") 
        png(filename, width=settings$width, height=settings$spectrumHeight, restoreConsole=FALSE)
      else 
        png(filename, width=settings$width, height=settings$spectrumHeight)  
      fun(zoom=zoom)    
      dev.off()    
      if(is.null(zoom) & !noCache) env$cache[[cacheName]][[index]] <- filename
      
    } else filename <- cache[[cacheName]][[index]]
    
    time <- proc.time()
    
    visible(plotObject) <- TRUE
    lim <- par()$usr
    rasterImage(readPNG(filename), lim[1], lim[3], lim[2], lim[4], interpolate=FALSE)
    if(verbose) cat(name, index, "loaded from cache in", (proc.time()-time)[3]*1000, "miliseconds\n")
    
  } else if(device=="cairo") {
    
    visible(plotObject) <- TRUE 
    fun(zoom=zoom)    
    
  } 
  
}

plotChromatogram <- function(zoom=NULL) {
  time <- proc.time()
  
  xic <- xic(n=1, settings$chromaMode)
  maxTotalIons <- max(xic[, 2])  
  
  xic <- xic[XICvalues, ]
  xic[, 2] <- xic[, 2]/maxTotalIons
  par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.75)
  plot(xic, type = "h", ylim=c(0, 1.075), xlab="Retention time", xlim=xLimitsXIC, 
       ylab=ifelse(settings$chromaMode, "Base peak intensity", "Total ion count"))
  lines(x=rep(spRtime(index), 2), y=c(0, 1), col="red", lty=3)
  text(x=xLimitsXIC[1], y=1.075, adj=0, 
       labels=paste("Max total ions ", maxTotalIons))
  abline(h=0, col="grey50")
  
  if(!is.null(zoom)) {
    lines(zoom$x[c(1, 1, 2, 2, 1)], zoom$y[c(1, 2, 2, 1, 1)], 
          col="grey25", lty=3)
  }
  if(verbose) 
    cat("    ", dim(xic)[1], "data points plotted in", (proc.time() - time)[3]*1000, "miliseconds\n")
}

plotSpectrumGraph <- function(zoom=NULL) {
  time <- proc.time()
  
  MsLevel <- spMsLevel(index)
  peaks <- peaks(index)
  basePeak <- max(peaks[, 2])
  peaks[, 2] <- peaks[, 2]/basePeak  
  
  # Make labels
  labels <- peaks[order(peaks[, 2], decreasing=TRUE), ]  
  xCoords <- labels[, 1]  
  n <- length(xCoords)
  current <- 1  
  out <- 1
  
  while(length(out) < settings$labelNumber & current < n) {
    s <- 1
    while(min(abs(xCoords[out] - xCoords[current+s])) <= settings$labelThreshold & current + s <= n) s <- s + 1
    out <- c(out, current <- current + s)
  }
  
  labels <- labels[out[out<=n], ]
    
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
  
  if(verbose) 
    cat("    ", dim(peaks)[1], "data points plotted in", (proc.time() - time)[3]*1000, "miliseconds\n")
}

plotSpectrumZoom <- function(limits=NULL) {
  
  if(verbose) cat("spectrum zoom limits:", limits$x, limits$y, "\n")
  
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
  abline(h=0, col="grey50")
  text(limits$x[2], limits$y[2], labels=paste("Acquisition number:", spIndex(index)))
}

plotChromaZoom <- function() {
  
  xic <- xic(n=1, settings$chromaMode)
  maxTotalIons <- max(xic[, 2]) 
  xic <- xic[XICvalues, ]
  xic[, 2] <- xic[, 2]/maxTotalIons
  
  par(mar=c(3, 3, 0, 0), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.65)
  plot(xic, type = "h", xlim=env$XICZoom$x, ylim=env$XICZoom$y, 
       xlab="Retention time", 
       ylab=ifelse(settings$chromaMode, "Base peak intensity", "Total ion count"))
  abline(h=0, col="grey50")
}
