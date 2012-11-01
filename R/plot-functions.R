plotXIC <- function() {
  
  # Make plot
  filename <- tempfile(fileext=".png")
  png(filename, width=500, height=250)
  
  dt <- xic(n=c(1, 2)[c(svalue(filterInfo$ms1), svalue(filterInfo$ms2))], 
            FALSE)
  par(mar=c(3,3,.5,1), mgp=c(2,.7,0), tck=-.01, bty="n", lab=c(5, 3, 7))
  
  time <- system.time({plot(dt, type = "l"); 
                      abline(v=spRtime(index), col="red", lty=3)})
  
  cat("\nxic:", dim(dt)[1], "data points plotted in", 
      time[3]*1000, "miliseconds")
  
  dev.off()
  
  visible(plotBottom) <- TRUE
  
#   svalue(plotBottom) <- filename
  
  grid::grid.raster(readPNG(filename))
}

plotSpectrum <- function() {
  filename <- tempfile(fileext=".png")
  png(filename, width=500, height=250)
  
  pks <- peaks(index)
  pksmax <- max(pks[, 2])
  
  mx <- pks[order(pks[, 2], decreasing=TRUE), ][1:5, ]
  pks[, 2] <- pks[, 2]/pksmax
    
  time <- proc.time()
  par(mar=c(3,3,.5,1), mgp=c(2,.7,0), tck=-.01, bty="n", lab=c(5, 3, 7))
  plot(pks, xlab="Intensity", ylab="M/Z", 
       type = ifelse(spMsLevel(index)==1, "l", "h"), 
                            xlim=c(min(spLowMZ()[spMsLevel()==spMsLevel(index)]), 
                                   max(spHighMZ()[spMsLevel()==spMsLevel(index)]))) 
  text(x=mx[, 1], y=mx[, 2]/pksmax, labels=round(mx[, 2], digits=3), col="grey50", adj=c(0, 0))
  time <- proc.time() - time
  
  cat("\nspectrum:", dim(pks)[1], "data points plotted in", time[3]*1000, "miliseconds")
  dev.off()
    
  visible(plotTop) <- TRUE
#   svalue(plotTop) <- filename
  
  grid::grid.raster(readPNG(filename))
  
}

