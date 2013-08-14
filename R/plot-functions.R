plotMsg <- function(message="Nothing to display", col="grey50") {
  par(mar=rep(0, 4))
  plot(0:1, 0:1, col="transparent", axes=FALSE)
  text(.5, .5, labels=message, cex=.65, col=col)
}

plotXIC <- function(zoom=NULL, int=NULL, noCache=FALSE, env)
  plotGeneric(zoom=zoom, int=int, cacheName="xic", 
              height=env$settings$chromaHeight, width=env$settings$width, 
              fun=plotChromatogram, plotObject=env$plotBottom, noCache=noCache, 
              name="chromatogram", env=env)

plotSpectrum <- function(zoom=NULL, int=NULL, noCache=FALSE, env) 
  plotGeneric(zoom=zoom, int=int, cacheName="spectra", 
              height=env$settings$spectrumHeight, width=env$settings$width, 
              fun=plotSpectrumGraph, plotObject=env$plotTop, noCache=noCache, 
              name="spectrum", env=env)

plotGeneric <- function(zoom=NULL, int=NULL, cache, cacheName, height, width, fun, 
                        plotObject, noCache=FALSE, name, env) {
  
  if(env$device=="png") {

    if(is.null(env$cache[[cacheName]][[env$index]]) | !is.null(zoom) | !is.null(int) | noCache) {
      if(is.null(zoom) & is.null(int) & env$verbose) cat("caching", name, env$index, "\n")
      
      time <- proc.time()
      
      filename <- tempfile(fileext=".png") 
      if(.Platform$OS.type=="windows") 
        png(filename, width=env$settings$width, height=env$settings$spectrumHeight, restoreConsole=FALSE)
      else 
        png(filename, width=env$settings$width, height=env$settings$spectrumHeight)  
      fun(zoom=zoom, int=int, env=env)    
      dev.off()    
      if(is.null(zoom) & is.null(int) & !noCache) env$cache[[cacheName]][[env$index]] <- filename
      
    } else filename <- env$cache[[cacheName]][[env$index]]
    
    time <- proc.time()
    
    visible(plotObject) <- TRUE
    lim <- par()$usr
    rasterImage(readPNG(filename), lim[1], lim[3], lim[2], lim[4], interpolate=FALSE)
    if(env$verbose) cat(name, env$index, "loaded from cache in", (proc.time()-time)[3]*1000, "miliseconds\n")
    
  } else if(env$device=="cairo") {
    
    visible(env$plotObject) <- TRUE 
    fun(zoom=zoom, int=int)    
    
  } 
  
}

plotChromatogram <- function(zoom=NULL, int=NULL, env) {
  time <- proc.time()
  
  xic <- env$xic(n=1, env$settings$chromaMode)
  maxTotalIons <- max(xic[, 2])  
  
  xic <- xic[env$XICvalues, ]
  xic[, 2] <- xic[, 2]/maxTotalIons
  par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.75)
  plot(xic, type = "h", ylim=c(0, 1.075), xlab="Retention time", xlim=env$xLimitsXIC, 
       ylab=ifelse(env$settings$chromaMode, "Base peak intensity", "Total ion count"), 
       xaxt = "n")
  axis(1, at=axTicks(1), labels=formatRt2(axTicks(1), env$settings$RtFormat, env$settings$digits))
  lines(x=rep(env$spRtime(env$index), 2), y=c(0, 1), col="red", lty=3)
  text(x=env$xLimitsXIC[1], y=1.075, adj=0, 
       labels=paste("Max total ions ", maxTotalIons))
  abline(h=0, col="grey50")
  
  if(!is.null(zoom)) {
    lines(zoom$x[c(1, 1, 2, 2, 1)], zoom$y[c(1, 2, 2, 1, 1)], 
          col="grey25", lty=3)
  }
  if(env$verbose) 
    cat("    ", dim(xic)[1], "data points plotted in", (proc.time() - time)[3]*1000, "miliseconds\n")
}

plotSpectrumGraph <- function(zoom=NULL, int=NULL, env) {
  time <- proc.time()
  MsLevel <- env$spMsLevel(env$index)
  peaks <- env$peaks(env$index)
  basePeak <- max(peaks[, 2])
  peaks[, 2] <- peaks[, 2]/basePeak  
  
  # Make labels
  if (env$settings$labelNumber > 1) {
    labels <- peaks[order(peaks[, 2], decreasing=TRUE), ]  
    xCoords <- labels[, 1]  
    n <- length(xCoords)
    current <- 1  
    out <- 1
    while(length(out) < env$settings$labelNumber & current < n) {
      s <- 1
      while(min(abs(xCoords[out] - xCoords[current+s])) <= env$settings$labelDist & current + s <= n) s <- s + 1
      out <- c(out, current <- current + s)
    }
    labels <- labels[out[out<=n], ]
  } else if (env$settings$labelNumber == 1) {
    labels <- peaks[which.max(peaks[, 2]), , drop=FALSE] 
  }
  
  par(mar=c(3,3,0,1), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.75)
  plot(x=env$xLimits[, MsLevel][1], y=0, xlab="Mass to charge ratio (M/Z)", ylab="Intensity", 
       xlim=env$xLimits[, MsLevel], ylim=c(0, 1.075)) 
  if(!is.null(int)) {
    polygon(int$x[c(1, 1, 2, 2, 1)], c(0, 1, 1, 0, 0), 
          col="lightpink1", border="transparent")
    text(env$xLimits[, MsLevel][2], 1.075, pos=2, 
         labels=paste("Integrand:", sum(peaks[peaks[, 1]>=int$x[1] & peaks[, 1]<=int$x[2], 2])))
  }
  lines(peaks, 
        type = ifelse(MsLevel==1, env$settings$MS1PlotType, env$settings$MS2PlotType))
  if (env$settings$labelNumber > 0)
    text(labels, labels=round(labels[, 1], env$settings$digits), 
         col="grey50", adj=c(0, 0))  
  par(adj=0)   
  text(x=min(env$xLimits[, MsLevel]), y=1.075, 
       labels=paste("Base peak intensity ", round(basePeak, env$settings$digits)))
  
  abline(h=0, col="grey50")
  
  if(!is.null(zoom)) {
    lines(zoom$x[c(1, 1, 2, 2, 1)], zoom$y[c(1, 2, 2, 1, 1)], 
          col="grey25", lty=3)
  }
  
  if(env$verbose) 
    cat("    ", dim(peaks)[1], "data points plotted in", (proc.time() - time)[3]*1000, "miliseconds\n")
}

plotSpectrumZoom <- function(limits=NULL, env) {
  
  if(env$verbose) cat("spectrum zoom limits:", limits$x, limits$y, "\n")
  
  peaks <- env$peaks(env$index)
  peaks[, 2] <- peaks[, 2]/max(peaks[, 2])  
  
  inFocus <- peaks[, 1] > limits$x[1] & peaks[, 1] < limits$x[2]
  inFocus <- inFocus & peaks[, 2] > limits$y[1] & peaks[, 2] < limits$y[2]
  numberInFocus <- sum(inFocus)
  
  par(mar=c(3, 3, 0, 0), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.65)
  plot(peaks, xlab="Mass to charge ratio (M/Z)", ylab="Intensity", 
       type = ifelse(env$spMsLevel(env$index)==1, env$settings$MS1PlotType, 
                     env$settings$MS2PlotType), xlim=limits$x, ylim=limits$y) 
  if(numberInFocus>0) {    
    labels <- peaks[inFocus, , drop=FALSE]
    labels <- labels[order(labels[, 2], decreasing=TRUE), , drop=FALSE]
    labels <- labels[1:min(numberInFocus, env$settings$labelNumber), , drop=FALSE]
    text(labels, labels=round(labels[, 1], digits=env$settings$digits), 
         col="grey50", adj=c(0, 0))    
  }
  par(adj=1)
  abline(h=0, col="grey50")
  text(limits$x[2], limits$y[2], labels=paste("Acquisition number:", env$spIndex(env$index)))
}

plotChromaZoom <- function(env) {
  
  xic <- env$xic(n=1, env$settings$chromaMode)
  maxTotalIons <- max(xic[, 2]) 
  xic <- xic[env$XICvalues, ]
  xic[, 2] <- xic[, 2]/maxTotalIons
  
  par(mar=c(3, 3, 0, 0), mgp=c(2,0.45,0), tck=-.01, bty="n", lab=c(5, 3, 7), 
      adj=.5, las=1, cex=0.65)
  plot(xic, type = "h", xlim=env$XICZoom$x, ylim=env$XICZoom$y, xaxt = 'n',
       xlab="Retention time", 
       ylab=ifelse(env$settings$chromaMode, "Base peak intensity", "Total ion count"))
  axis(1, at=axTicks(1), labels=formatRt2(axTicks(1), env$settings$RtFormat, env$settings$digits))
  abline(h=0, col="grey50")
}
