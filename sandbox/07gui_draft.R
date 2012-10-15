library(gWidgets)
library(mzR)
# source("../R/middle-layer.R")

# filename <- "d:/Dropbox/Documents/cambridge/r project/Data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"

options(guiToolkit = "RGtk2")
# options(guiToolkit = "tcltk")
 
.msGUIenv <- new.env(parent = emptyenv(), hash = TRUE)
assign("file.types", 
       list("Mass spectrometry data files" = list(patterns = c("*.netCDF", "*.mzXML", "*.mzData", "*.mzML")), 
            "R data files" = list(patterns = c("*.RData","*.rda")), 
            "All files" = list(patterns = c("*"))), 
       pos=.msGUIenv)
assign("digits", 4, pos=.msGUIenv)
assign("textHead", list(weight="bold", family="sans", 
                                size=12, color="grey40"), pos=.msGUIenv)
assign("textReg", list(weight="normal", family="sans", 
                                size=8, color="grey10"), pos=.msGUIenv)
assign("textRegGrey", list(weight="normal", family="sans", 
                       size=8, color="grey50"), pos=.msGUIenv)
assign("textDeemp", list(weight="normal", family="sans", 
                               size=8, color="grey50"), pos=.msGUIenv)

### Settings
# msGUIsettings <- list()
# 
# msGUIsettings <- list(
#   list("digits", 2, "Number of digits displayed"), 
#   list("textHead", list(weight="bold", family="sans", 
#                         size=12, color="grey04"), 
#        "Font settings", "Headings")
#   list("textReg", list(weight="normal", family="sans", 
#                        size=8, color="grey10"))


drawMain <- function(updateOnLoad=FALSE) {
#   )

  
  mget <- function(what) get(what, pos=.msGUIenv)
  mset <- function(what, value) assign(what, value, .msGUIenv)
    
  updateExperimentInfo <- function() {
    rtRange <- round(.msGUIenv$expRtRange(), digits=.msGUIenv$digits)
    precMzRange <- round(.msGUIenv$expPrecMzRange(), digits=.msGUIenv$digits)
    svalue(expInfo$rtfrom) <- rtRange[1]
    svalue(expInfo$rtto) <- rtRange[2]
    svalue(expInfo$pmzfrom) <- precMzRange[1]
    svalue(expInfo$pmzto) <- precMzRange[2]
    svalue(expInfo$ms1) <- .msGUIenv$nMS1()
    svalue(expInfo$ms2) <- .msGUIenv$nMS2()
  }
  
  counterReset <- function() {    
    .msGUIenv$currSequence <- .msGUIenv$spIndex()
    .msGUIenv$counter <- 1
  }
  
  updateExperiment <- function() {
    
    updateExperimentInfo()
    counterReset()    
    filterReset()    
    filterSwitch(TRUE)
    
    plotXIC()
    updateSpectrum()
  }
  
  plotXIC <- function(ms=NULL) {    
    par(mar=rep(0, 4))
#     if(svalue(.msGUIenv$))
    time <- system.time(plot(.msGUIenv$xic(ms), type = "l", frame.plot=FALSE, axes=TRUE, xlab="", ylab=""))
    cat("\ndatapoints:", dim(.msGUIenv$xic(ms))[1], "xic plot drawn in:", time[3])
  }
  
  updateSpectrum <- function(h=list(action=0), ...) {
    .msGUIenv$counter <- .msGUIenv$counter + h$action
    index <- .msGUIenv$currSequence[.msGUIenv$counter]
    svalue(specInfo$rt) <- .msGUIenv$spRtime(index)
    svalue(specInfo$ind) <- paste(.msGUIenv$spIndex(index), "of", length(.msGUIenv$currSequence))
    svalue(specInfo$pmz) <- round(.msGUIenv$spPrecMz(index), digits=.msGUIenv$digits)
    svalue(specInfo$int) <- round(.msGUIenv$spPrecInt(index), digits=.msGUIenv$digits)
    
    if(length(.msGUIenv$currSequence)==1) {
      enabled(buttonLeft) <- FALSE
      enabled(buttonRight) <- FALSE      
    } else 
    if(.msGUIenv$counter==1) {
      enabled(buttonLeft) <- FALSE
      enabled(buttonRight) <- TRUE
    } else 
    if(.msGUIenv$counter==length(.msGUIenv$currSequence)) {
      enabled(buttonLeft) <- TRUE
      enabled(buttonRight) <- FALSE
    } else {
      enabled(buttonLeft) <- TRUE
      enabled(buttonRight) <- TRUE
    }
    
    visible(plotTop) <- TRUE
    par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01) # some random settings for now
    time <- system.time(plot(.msGUIenv$peaks(index), type = "l", frame.plot=FALSE, axes=TRUE, xlab="", ylab=""))
    cat("\ndatapoints:", dim(.msGUIenv$peaks(index))[1], "spectrum drawn in:", time[3])
  }
    
  btwn <- function(x, from, to) {
    stopifnot(all(is.numeric(x), is.character(from), is.character(to)))
    x >= as.numeric(from) & x <= as.numeric(to)
  }
  
  validityCheck <- function(from, to, data, prevFrom, prevTo) {
    numFrom <- as.numeric(svalue(from))
    numTo <- as.numeric(svalue(to))
    prevValid <- prevFrom >= min(data) & prevFrom < prevTo & prevTo > min(data) & prevTo <= max(data)
    if(numFrom < min(data)) svalue(from) <- min(data)
    if(numFrom > max(data)) svalue(from) <- ifelse(prevValid, prevFrom, min(data))
    if(numTo < min(data)) svalue(to) <- ifelse(prevValid, prevTo, max(data))
    if(numTo > max(data)) svalue(to) <- max(data)
    if(numFrom > numTo) {
      if(numFrom!=prevFrom) svalue(from) <- prevFrom
      else svalue(to) <- prevTo
    }
  }
  
  filterSpectra <- function(h, ...) {    
    keep <- TRUE
        
    mapply(blockHandler, filterInfo, filterSpectraHandlerIDs)
    
    validityCheck(filterInfo$rtfrom, filterInfo$rtto, 
                  .msGUIenv$spRtime(), .msGUIenv$filterValues$rtfrom, .msGUIenv$filterValues$rtto)
    validityCheck(filterInfo$pmzfrom, filterInfo$pmzto, 
                  .msGUIenv$spPrecMz(), .msGUIenv$filterValues$pmzfrom, .msGUIenv$filterValues$pmzto)
    validityCheck(filterInfo$spifrom, filterInfo$spito, 
                  .msGUIenv$spPrecInt(), .msGUIenv$filterValues$spifrom, .msGUIenv$filterValues$spito)
    
    if(svalue(filterInfo$rtactive)) 
      keep <- keep & btwn(.msGUIenv$spRtime(), 
                          svalue(filterInfo$rtfrom), svalue(filterInfo$rtto))    
    if(svalue(filterInfo$pmzactive)) 
      keep <- keep & btwn(.msGUIenv$spPrecMz(), 
                          svalue(filterInfo$pmzfrom), svalue(filterInfo$pmzto))
    if(svalue(filterInfo$spiactive)) 
      keep <- keep & btwn(.msGUIenv$spPrecInt(), 
                          svalue(filterInfo$spifrom), svalue(filterInfo$spito))
    if(!svalue(filterInfo$ms1)) keep <- keep & .msGUIenv$spMsLevel() != 1
    if(!svalue(filterInfo$ms2)) keep <- keep & .msGUIenv$spMsLevel() != 2
    
    saveFilterValues()
    
    mapply(unblockHandler, filterInfo, filterSpectraHandlerIDs)
    
    newSequence <- .msGUIenv$spIndex()[keep]
    if(!identical(.msGUIenv$currSequence, newSequence)) {
      .msGUIenv$currSequence <- newSequence
      .msGUIenv$counter <- 1
      updateSpectrum()
    } 
    
#     if(!identical(c(.msGUIenv$filterValues$ms1, .msGUIenv$filterValues$ms2), 
#                   c(svalue(filterInfo$ms1), svalue(filterInfo$ms1)))) plotXIC()
    # to be implemented: when MSn selection changes, XIC should be redrawn accordingly
  }
  
  saveFilterValues <- function() {
    .msGUIenv$filterValues <- list(rtfrom=svalue(filterInfo$rtfrom), 
                                   rtto=svalue(filterInfo$rtto), 
                                   pmzfrom=svalue(filterInfo$pmzfrom), 
                                   pmzto=svalue(filterInfo$pmzto), 
                                   spifrom=svalue(filterInfo$spifrom), 
                                   spito=svalue(filterInfo$spito), 
                                   ms1=svalue(filterInfo$ms1), 
                                   ms2=svalue(filterInfo$ms2))
    .msGUIenv$filterValues <- lapply(.msGUIenv$filterValues, as.numeric)
    }
    
  filterReset <- function() {
    mapply(blockHandler, filterInfo, filterSpectraHandlerIDs)
    svalue(filterInfo$rtactive) <- FALSE
    svalue(filterInfo$rtfrom) <- min(.msGUIenv$spRtime())
    svalue(filterInfo$rtto) <- max(.msGUIenv$spRtime())
    svalue(filterInfo$pmzactive) <- FALSE
    svalue(filterInfo$pmzfrom) <- min(.msGUIenv$spPrecMz())
    svalue(filterInfo$pmzto) <- max(.msGUIenv$spPrecMz())
    svalue(filterInfo$spiactive) <- FALSE
    svalue(filterInfo$spifrom) <- min(.msGUIenv$spPrecInt())
    svalue(filterInfo$spito) <- max(.msGUIenv$spPrecInt())
    svalue(filterInfo$ms1) <- TRUE
    svalue(filterInfo$ms2) <- TRUE
    saveFilterValues()
    mapply(unblockHandler, filterInfo, filterSpectraHandlerIDs)
  }
  
  filterSwitch <- function(on) {
    lapply(filterInfo, function(x) enabled(x) <- on)
  }
    
  openFileHandler <- function(h, ...) {
    filename <- gfile(text = "Open file", type = "open", 
                      filter = .msGUIenv$file.types, handler=function(h, ...) NULL)
    if(!is.na(filename)) {
      makeMzRrampAccessors(filename)
      updateExperiment()      
    }
  }
  
  # Window and structure
  msGUIWindow <- gwindow("msGUI", visible=FALSE)
  groupMain <- ggroup(container=msGUIWindow, horizontal=FALSE, expand=TRUE)
  groupUpper <- ggroup(container=groupMain)
  groupMiddle <- ggroup(container=groupMain)
  groupMiddleLeft <- ggroup(container=groupMiddle, horizontal=FALSE, spacing=5)
  groupMiddleRight <- gframe(container=groupMiddle, horizontal=FALSE, expand=TRUE)
  groupPlots <- gpanedgroup(container=groupMiddleRight, horizontal=FALSE, expand=TRUE)
  
  ## Top group
  buttonOpenFile <- gbutton("Open file", container=groupUpper, 
                            handler=openFileHandler)
  buttonOpenObject <- gbutton("Open R object", container=groupUpper)
  addSpring(groupUpper)
  buttonSettings <- gbutton("Settings", container=groupUpper)
  buttonHelp <- gbutton("Help", container=groupUpper)
    
  # Left pane content
  
  # Lists for formatting and retrieval
  headings <- list()
  regular <- list()
  deemp <- list()
  separator <- list()
  expInfo <- list()
  specInfo <- list()
  filterInfo <- list()
  
  le <- glayout(container=groupMiddleLeft, spacing=2)
  
  #Experiment info
  
  le[1, 1:3] <- (separator$t1 <- glabel("", container=le))
  
  le[2, 1:3] <- (headings$t1 <- glabel("Experiment information", container=le))
  le[3, 2, anchor=c(-1,-1)] <- (deemp$t1 <- glabel("from", container=le))
  le[3, 3, anchor=c(-1,-1)] <- (deemp$t2 <- glabel("to", container=le))
  le[4, 1] <- (regular$t1 <- glabel("Retention time", container=le))
  le[4, 2] <- (expInfo$rtfrom <- glabel(container=le))
  le[4, 3] <- (expInfo$rtto <- glabel(container=le))
  le[5, 1] <- (regular$t2 <- glabel("Precursor M/Z", container=le))
  le[5, 2] <- (expInfo$pmzfrom <- glabel(container=le))
  le[5, 3] <- (expInfo$pmzto <- glabel(container=le))
  le[6, 1] <- (regular$t3 <- glabel("Spectra count", container=le))
  le[6, 2] <- (regular$t4 <- glabel("MS1", container=le))
  le[6, 3] <- (expInfo$ms1 <- glabel(container=le))
  le[7, 2] <- (regular$t15 <- glabel("MS2", container=le))
  le[7, 3] <- (expInfo$ms2 <- glabel(container=le))
  
  le[8, 1:3] <- (separator$t2 <- glabel("", container=le))
  
  # Spectrum Info
  
  le[9, 1:3] <- (headings$t2 <- glabel("Spectrum information", container=le))
  le[10, 1] <- (regular$t5 <- glabel("Retention time", container=le))
  le[10, 2:3] <- (specInfo$rt <- glabel("", container=le))
  le[11, 1] <- (regular$t6 <- glabel("Index", container=le))
  le[11, 2:3] <- (specInfo$ind <- glabel("", container=le, editable=TRUE))
  le[12, 1] <- (regular$t7 <- glabel("Precursor M/Z", container=le))
  le[12, 2:3] <- (specInfo$pmz <- glabel("", container=le))
  le[13, 1] <- (regular$t8 <- glabel("Precursor intensity", container=le))
  le[13, 2:3] <- (specInfo$int <- glabel("", container=le))

  le[14, 1:5] <- (separator$t3 <- glabel("", container=le))
  
  # Filter
  
  le[15, 1] <- (headings$t3 <- glabel("Filter", container=le))
  le[15, 2, anchor=c(-1,-1)] <- (deemp$t5 <- glabel("from", container=le))
  le[15, 3, anchor=c(-1,-1)] <- (deemp$t6 <- glabel("to", container=le))
  
  le[16, 1] <- (filterInfo$rtactive <- gcheckbox("Retention time", container=le, checked=TRUE))
  le[16, 2] <- (filterInfo$rtfrom <- gedit("", container=le, width=5))
  le[16, 3] <- (filterInfo$rtto <- gedit("", container=le, width=5))
  le[17, 1] <- (filterInfo$pmzactive <- gcheckbox("Precursor M/Z", container=le, checked=TRUE))
  le[17, 2] <- (filterInfo$pmzfrom <- gedit("", container=le, width=5))
  le[17, 3] <- (filterInfo$pmzto <- gedit("", container=le, width=5))
  le[18, 1] <- (filterInfo$spiactive <- gcheckbox("Precursor intensity", container=le, checked=TRUE))
  le[18, 2] <- (filterInfo$spifrom <- gedit("", container=le, width=5))
  le[18, 3] <- (filterInfo$spito <- gedit("", container=le, width=5))
  
  le[19, 1:5] <- (separator$t3 <- glabel("", container=le))
  
  le[20, 1, anchor=c(-1,0)] <- (regular$t12 <- glabel("Display MS levels", container=le))
  le[20, 2] <- (filterInfo$ms1 <- gcheckbox("MS1", checked=TRUE, container=le))
  le[20, 3] <- (filterInfo$ms2 <- gcheckbox("MS2", checked=TRUE, container=le))
  
  groupLeftButtons <- ggroup(container=groupMiddleLeft, expand=TRUE)
  buttonLeft <- gbutton(text=gettext("Previous"), handler=updateSpectrum, 
                         action=-1, cont=groupLeftButtons)
  
  buttonRight <- gbutton(text=gettext("Next"), handler=updateSpectrum, 
                             action=1, cont=groupLeftButtons)
  
  plotTop <- ggraphics(container=groupPlots, width=400, height=250, dpi=72)
  plotBottom <- ggraphics(container=groupPlots, width=400, height=250, dpi=72)

  lapply(headings, function(x) font(x) <- .msGUIenv$textHead)
  lapply(regular, function(x) font(x) <- .msGUIenv$textReg)
  lapply(deemp, function(x) font(x) <- .msGUIenv$textDeemp)
  lapply(expInfo, function(x) font(x) <- .msGUIenv$textReg)
  lapply(specInfo, function(x) font(x) <- .msGUIenv$textReg)
  lapply(filterInfo, function(x) font(x) <- .msGUIenv$textReg)
  lapply(separator, function(x) font(x) <- list(size=4))
  
  filterSpectraHandlerIDs <- lapply(filterInfo, addHandlerChanged, handler=filterSpectra)
  
  visible(msGUIWindow) <- TRUE
  
  enabled(buttonLeft) <- FALSE
  enabled(buttonRight) <- FALSE
  filterSwitch(FALSE)
  
  if(updateOnLoad) updateExperiment()
}
