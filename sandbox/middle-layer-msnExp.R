makeDataEnvMSnExp <- function(object) {
  ## Package environment to store private data
  ## Arguments
  ##  filename: name of a file to be opened by
  ##            mzR::openMSfile
  ## Value:
  ##  An environment storing the filehandle,
  ##  the data header and runInfo
  e <- new.env(parent = emptyenv(), hash = TRUE)   
  assign("fh", assayData(object), env = e)
  hd <- header(object)
  assign("hd", hd[order$header(object), ], env = e)
  assign("runInfo", list(scanCount=length(object), 
                         lowMZ=min(e$hd$precursor.mz), 
                         highMZ=max(e$hd$precursor.mz), 
                         dStartTime=min(e$hd$retention.time), 
                         dEndTime=max(e$hd$retention.time), 
                         msLevels=unique(e$hd$ms.level)), env = e)
  lockEnvironment(e, bindings = TRUE)
  return(e)
}
makeMSnExpAccessors <- function(object, assignEnv) {
  ## Creates the appropriate raw data file accessors
  ## Arguments:
  ##  filename: a character (of length 1), containing
  ##            the name of the raw file to be accessed
  dataEnv <- makeDataEnvMSnExp(object)
  if (missing(assignEnv)) 
    assignEnv <- sys.frame(-1)
  makeAccessors(assignEnv)
  ## Experiment info
  assign("expRtRange",
         function() c(dataEnv$runInfo$dStartTime, dataEnv$runInfo$dEndTime),
         pos = assignEnv)
  assign("expPrecMzRange",
         function() {
           pmz <- dataEnv$hd$precursor.mz
           range(pmz[pmz != 0])
         }, pos = assignEnv)
  assign("nMS1",
         function(n = 1) sum(dataEnv$hd$ms.level == n),
         pos = assignEnv)
  assign("nMS2",
         function(n = 2) sum(dataEnv$hd$ms.level == n),
         pos = assignEnv)
  ## Spectrum info
  assign("spRtime",
         function(i) dataEnv$hd[i, "retention.time"]  ,
         pos = assignEnv)
  assign("spIndex",
         function(i) dataEnv$hd[i, "acquisition.number"],
         pos = assignEnv)
  assign("spPrecMz", 
         function(i) dataEnv$hd[i, "precursor.mz"],
         pos = assignEnv)
  assign("spPrecInt",
         function(i) dataEnv$hd[i, "precursor.intensity"] ,
         pos = assignEnv)
  assign("spMsLevel",
         function(i) dataEnv$hd[i, "ms.level"] ,
         pos = assignEnv)
  assign("spPrecCharge",
         function(i) 0, #rep(0, dataEnv$runInfo$scanCount), #dataEnv$hd[i, "precursorCharge"],
         pos = assignEnv)
#   assign("spPrecScanNum",
#          function(i) dataEnv$hd[i, "precursorScanNum"],
#          pos = assignEnv)
#   assign("spBasePeakMz",
#          function(i) dataEnv$hd[i, "basePeakMZ"],
#          pos = assignEnv)
  assign("spLowMZ",
         function(i) min(mz(get(rownames(dataEnv$hd)[i], envir=dataEnv$fh))),
         pos = assignEnv)
  assign("spHighMZ",
         function(i) max(mz(get(rownames(dataEnv$hd)[i], envir=dataEnv$fh))),
         pos = assignEnv)
#   assign("spBasePeakInt", 
#          function(i) dataEnv$hd[i, "basePeakIntensity"],
#          pos = assignEnv)
  assign("spPeaksCount",
         function(i) peaksCount(get(rownames(dataEnv$hd)[i], envir=dataEnv$fh)),
         pos = assignEnv)
  ## Plotting
  assign("peaks",
         function(i) cbind(mz(get(rownames(dataEnv$hd)[i], envir=dataEnv$fh)), 
                           intensity(get(rownames(dataEnv$hd)[i], envir=dataEnv$fh))),
         pos = assignEnv)
#   assign("xic",
#          function(n = 1, basePeaks=FALSE) {    
#            if (is.null(n)) {
#              lev <- TRUE
#            } else {
#              lev <- dataEnv$hd$msLevel %in% n
#            }
#            dataEnv$hd[lev, c("retentionTime",
#                              ifelse(basePeaks, "basePeakIntensity", 
#                                     "totIonCurrent"))]
#          }, pos = assignEnv)
#   assign("close",
#          function() mzR::close(dataEnv$fh),
#          pos = assignEnv)
  invisible(TRUE)
}
