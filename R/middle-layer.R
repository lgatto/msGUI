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
               "spPeaksCount",
               ## Plotting
               "peaks",
               "xic")

## Package environment to store private data
.msGUIenv <- new.env(parent = emptyenv(), hash = TRUE) 
populateMsGuiEnv <- function(filename) {
  assign("fh", openMSfile(filename), env = .msGUIenv)
  assign("hd", header(.msGUIenv$fh), env = .msGUIenv)
  assign("runInfo", runInfo(.msGUIenv$fh), env = .msGUIenv)
  lockEnvironment(.msGUIenv, bindings = TRUE)
  invisible(TRUE)
}
  
## Initialise
emptyAccessors <- function() {
  acc <- lapply(1:length(accessors),
                function(i) return(function() stop("Not yet implemented")))
  names(acc) <- accessors
  return(acc)
}


makeMSnExpAccessors <- function() {
  emptyAccessors()
}

makeMzRrampAccessors <- function() {
  mzRrampAccessors <- emptyAccessors()  
  ## Experiment info
  mzRrampAccessors[["expRtRange"]] <-
    function() c(.msGUIenv$runInfo$dStartTime, .msGUIenv$runInfo$dEndTime)
  mzRrampAccessors[["expPrecMzRange"]] <-
    function() {
      pmz <- .msGUIenv$hd$precursorMZ
      range(pmz[pmz != 0])
    }
  mzRrampAccessors[["nMS1"]] <-
    function(n = 1) sum(.msGUIenv$hd$msLevel == n)
  mzRrampAccessors[["nMS2"]] <-
    function(n = 2) sum(.msGUIenv$hd$msLevel == n)  
  ## Spectrum info
  mzRrampAccessors[["spRtime"]] <-
    function(i) .msGUI$hd[i, "retentionTime"]  
  mzRrampAccessors[["spIndex"]] <-
    function(i) .msGUI$hd[i, "acquisitionNum"]  
  mzRrampAccessors[["spPrecMz"]] <- 
    function(i) .msGUI$hd[i, "precursorMZ"]  
  mzRrampAccessors[["spPrecInt"]] <-
    function(i) .msGUI$hd[i, "precursorIntensity"]  
  mzRrampAccessors[["spMsLevel"]] <-
    function(i) .msGUI$hd[i, "msLevel"]  
  mzRrampAccessors[["spPrecCharge"]] <-
    function(i) .msGUI$hd[i, "precursorCharge"]
  mzRrampAccessors[["spPrecScanNum"]] <-
    function(i) .msGUI$hd[i, "precursorScanNum"]
  mzRrampAccessors[["spBasePeakMz"]] <-
    function(i) .msGUI$hd[i, "basePeakMZ"]
  mzRrampAccessors[["spBasePeakInt"]] <- 
    function(i) .msGUI$hd[i, "basePeakIntensity"]
  mzRrampAccessors[["spPeaksCount"]] <-
    function(i) .msGUI$hd[i, "peaksCount"]  
  ## Plotting
  mzRrampAccessors[["peaks"]] <-
    function(i) peaks(.msGUIenv$fh, i)
  mzRrampAccessors[["xic"]] <-
    function(n = 1) {    
      if (is.null(n)) {
        lev <- TRUE
      } else {
        lev <- .msGUIenv$hd$msLevel %in% n
      }
      .msGUIenv$hd[lev, c("retentionTime",
                          "totIonCurrent")]
    }
  return(mzRrampAccessors)
}



