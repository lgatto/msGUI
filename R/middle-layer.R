## We want to be able to access data from the following backends
## - raw data files via mzRramp instances - mzR package
## - raw experiments via MSnExp instances - MSnbase package
backends <- c("mzRramp", "MSnExp")

## List of data accessors that we need
## (to be completed with Andrius' list)
accessors <- c(## Experiment info
               "expRtRange",
               "nMS1", "nMS2", ## possibly other
               "expPrecMzRange",
               ## Spectrum info
               "spRtime",
               "spIndex",
               "spPrecMz",
               "spPrecInt",
               "spMsLevel",
               "spPrecCharge",
               "spPrecScanNum",
               "spBasePeakMz",
               "spBasePeakInt",
               "spLowMZ", 
               "spHighMZ",
               "spPeaksCount",
               ## Plotting
               "peaks",
               "xic",
               ## Housekeeing
               "close")

makeDataEnv <- function(filename) {
  ## Package environment to store private data
  ## Arguments
  ##  filename: name of a file to be opened by
  ##            mzR::openMSfile
  ## Value:
  ##  An environment storing the filehandle,
  ##  the data header and runInfo
  e <- new.env(parent = emptyenv(), hash = TRUE)   
  assign("fh", openMSfile(filename), envir = e)
  assign("hd", header(e$fh), envir = e)
  assign("runInfo", runInfo(e$fh), envir = e)
  lockEnvironment(e, bindings = TRUE)
  return(e)
}


makeAccessors <- function(assignEnv) {
  ## Initialise accessors to 'not yet implemented'
  ## Arguments:
  ##  assignEnv: an environment, where the empty  
  ##             default accessors will be assigned
  ##             to. If missing, defaults to calling
  ##             environment.
  ## Value:
  ##  invisibly returns TRUE
  if (missing(assignEnv))
    assignEnv <- sys.frame(-1)
  tmp <- sapply(accessors,
                assign, value = function(x) stop("Not yet implemented"),
                pos = assignEnv)
  invisible(TRUE)
}


makeMzRrampAccessors <- function(filename, assignEnv) {
  ## Creates the appropriate raw data file accessors
  ## Arguments:
  ##  filename: a character (of length 1), containing
  ##            the name of the raw file to be accessed
  dataEnv <- makeDataEnv(filename)
  if (missing(assignEnv)) 
    assignEnv <- sys.frame(-1)
  makeAccessors(assignEnv)
  ## Experiment info
  assign("expRtRange",
         function() c(dataEnv$runInfo$dStartTime, dataEnv$runInfo$dEndTime),
         pos = assignEnv)
  assign("expPrecMzRange",
         function() {
           pmz <- dataEnv$hd$precursorMZ
           range(pmz[pmz != 0])
         }, pos = assignEnv)
  assign("nMS1",
         function(n = 1) sum(dataEnv$hd$msLevel == n),
         pos = assignEnv)
  assign("nMS2",
         function(n = 2) sum(dataEnv$hd$msLevel == n),
         pos = assignEnv)
  ## Spectrum info
  assign("spRtime",
         function(i) dataEnv$hd[i, "retentionTime"]  ,
         pos = assignEnv)
  assign("spIndex",
         function(i) dataEnv$hd[i, "acquisitionNum"],
         pos = assignEnv)
  assign("spPrecMz", 
         function(i) dataEnv$hd[i, "precursorMZ"],
         pos = assignEnv)
  assign("spPrecInt",
         function(i) dataEnv$hd[i, "precursorIntensity"] ,
         pos = assignEnv)
  assign("spMsLevel",
         function(i) dataEnv$hd[i, "msLevel"] ,
         pos = assignEnv)
  assign("spPrecCharge",
         function(i) dataEnv$hd[i, "precursorCharge"],
         pos = assignEnv)
  assign("spPrecScanNum",
         function(i) dataEnv$hd[i, "precursorScanNum"],
         pos = assignEnv)
  assign("spBasePeakMz",
         function(i) dataEnv$hd[i, "basePeakMZ"],
         pos = assignEnv)
  assign("spLowMZ",
         function(i) dataEnv$hd[i, "lowMZ"],
         pos = assignEnv)
  assign("spHighMZ",
         function(i) dataEnv$hd[i, "highMZ"],
         pos = assignEnv)
  assign("spBasePeakInt", 
         function(i) dataEnv$hd[i, "basePeakIntensity"],
         pos = assignEnv)
  assign("spPeaksCount",
         function(i) dataEnv$hd[i, "peaksCount"],
         pos = assignEnv)
  ## Plotting
  assign("peaks",
         function(i) mzR::peaks(dataEnv$fh, i),
         pos = assignEnv)
  assign("xic",
         function(n = 1, basePeaks=FALSE) {    
           if (is.null(n)) {
             lev <- TRUE
           } else {
             lev <- dataEnv$hd$msLevel %in% n
           }
           dataEnv$hd[lev, c("retentionTime",
                               ifelse(basePeaks, "basePeakIntensity", 
                                                 "totIonCurrent"))]
         }, pos = assignEnv)
  assign("close",
         function() mzR::close(dataEnv$fh),
         pos = assignEnv)
  invisible(TRUE)
}


makeMSnExpAccessors <- function() {
  ## Creates the MSnExp accessors - currently default ones
  assignEnv <- sys.frame(-1)
  makeAccessors(assignEnv)
  invisible(TRUE)
}


makeDataEnvMSnExp <- function(object) {
  ## Package environment to store private data
  ## Arguments
  ## object: name of object of class MSnExp
  ## Value:
  ##  An environment storing the filehandle,
  ##  the data header and runInfo
  e <- new.env(parent = emptyenv(), hash = TRUE)   
  assign("fh", assayData(object), envir = e)
  hd <- header(object)
  assign("hd", hd[order(hd$acquisition.number), ], envir = e)
  assign("runInfo", list(scanCount=length(object), 
                         lowMZ=min(e$hd$precursor.mz), 
                         highMZ=max(e$hd$precursor.mz), 
                         dStartTime=min(e$hd$retention.time), 
                         dEndTime=max(e$hd$retention.time), 
                         msLevels=unique(e$hd$ms.level)), envir = e)
  e$hd[, c("lowMZ", "highMZ")] <- t(sapply(rownames(hd), 
                                           function(x) 
                                             range(mz(get(x, envir=e$fh)))))
  lockEnvironment(e, bindings = TRUE)
  return(e)
}


makeMSnExpAccessors <- function(object, assignEnv) {
  ## Creates the appropriate raw data file accessors
  ## Arguments:
  ## object: name of object of class MSnExp
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
         function(i) dataEnv$hd[i, "retention.time"],
         pos = assignEnv)
  assign("spIndex",
         function(i) dataEnv$hd[i, "acquisition.number"],
         pos = assignEnv)
  assign("spPrecMz", 
         function(i) dataEnv$hd[i, "precursor.mz"],
         pos = assignEnv)
  assign("spPrecInt",
         function(i) dataEnv$hd[i, "precursor.intensity"],
         pos = assignEnv)
  assign("spMsLevel",
         function(i) dataEnv$hd[i, "ms.level"],
         pos = assignEnv)
  assign("spPrecCharge",
         function(i) dataEnv$hd[i, "charge"],
         pos = assignEnv)
  #   assign("spPrecScanNum",
  #          function(i) dataEnv$hd[i, "precursorScanNum"],
  #          pos = assignEnv)
  #   assign("spBasePeakMz",
  #          function(i) dataEnv$hd[i, "basePeakMZ"],
  #          pos = assignEnv)
  assign("spLowMZ",
         function(i) dataEnv$hd[i, "lowMZ"],
         pos = assignEnv)
  assign("spHighMZ",
         function(i) dataEnv$hd[i, "highMZ"], 
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
