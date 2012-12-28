wrapper <- function(filename=NULL, object=NULL, device="png", verbose=FALSE) {  
  
  env <- environment()
  
  settings <- defaultSettings()
  
  options(guiToolkit = "RGtk2")
  
  zoomWindowClosed <- TRUE
  XICWindowClosed <- TRUE
  optionsWindowClosed <- TRUE
  spectrumZoom <- NULL
  XICZoom <- NULL
  anyMS1spectra <- TRUE
  XICvalues <- TRUE
  experimentLoaded <- FALSE
  closingMsGUI <- FALSE
  
  environment(openFileHandler) <- env  
  environment(drawMain) <- env  
  environment(updateExperiment) <- env
  environment(updateExperimentInfo) <- env
  environment(updateSpectrum) <- env  
  environment(updateSpectrumInfo) <- env
  environment(counterReset) <- env
  environment(saveFilterValues) <- env
  environment(validityCheck) <- env
  environment(filterReset) <- env
  environment(filterSpectra) <- env
  environment(filterXIC) <- env
  environment(filterSwitch) <- env
  environment(clickSwitch) <- env
  environment(initialiseGUI) <- env
  environment(plotXIC) <- env
  environment(plotSpectrum) <- env  
  environment(plotGeneric) <- env
  environment(plotSpectrumGraph) <- env 
  environment(plotChromatogram) <- env 
  environment(buttonSwitch) <- env
  environment(plotSpectrumZoom) <- env
  environment(blockFilters) <- env
  environment(getObjects) <- env
  environment(openObject) <- env
  environment(drawVarBrowser) <- env
  environment(plotChromaZoom) <- env
  environment(resetCache) <- env
  environment(drawOptions) <- env
  environment(optsHandlerDefaults) <- env
  
  drawMain(env)
  initialiseGUI()
  
  if(!is.null(filename)) {
    makeMzRrampAccessors(filename, env)
    updateExperiment(env)
    experimentLoaded <- TRUE
  }   
  
  if(!is.null(object)) {
    makeMSnExpAccessors(object, env)
    updateExperiment(env)
    experimentLoaded <- TRUE
  }
}

updateExperimentInfo <- function() {
  precMzRange <- round(expPrecMzRange(), digits=settings$digits)
  svalue(expInfo$rtfrom) <- formatRt(expRtRange()[1])
  svalue(expInfo$rtto) <- formatRt(expRtRange()[2])
  svalue(expInfo$pmzfrom) <- precMzRange[1]
  svalue(expInfo$pmzto) <- precMzRange[2]
  svalue(expInfo$ms1) <- nMS1()
  svalue(expInfo$ms2) <- nMS2()
}

counterReset <- function(env) {    
  env$currSequence <- spIndex()
  env$counter <- 1
}

updateExperiment <- function(env) {  
  updateExperimentInfo()
  
  env$filterData <- list(spRtime(), spIndex(), spPrecMz(), 
                         spPrecInt(), spPrecCharge(), spPrecMz()*spPrecCharge())
  # filterData stores data for filters for fast access. 
  
  env$nSpectra <- length(filterData[[1]])
  env$xLimits <- sapply(1:2, function(i) if(any(spMsLevel()==i)) 
    c(min(spLowMZ()[spMsLevel()==i]), 
      max(spHighMZ()[spMsLevel()==i])) 
                        else rep(NA, 2))
  env$anyMS1spectra <- any(spMsLevel()==1)
  if(anyMS1spectra) {
    dt <- cbind(spMsLevel(), spPrecMz())
    where <- which(dt[, 1]==1)
    frame <- cbind(where[-length(where)] + 1, where[-1] - 1)
    env$MS2indices <- frame[frame[, 1] < frame[, 2], ]
    
    env$xLimitsXIC <- range(xic()[, 1])
  } else {
    visible(plotBottom) <- TRUE
    plotMsg("Experiment contains no MS1 spectra.")
  }
  
  resetCache()
  counterReset(env)    
  filterReset(env)    
  filterSwitch(TRUE)
  clickSwitch(TRUE)
  updateSpectrum()  
}

updateSpectrumInfo <- function() {
  svalue(specInfo$rt) <- formatRt(spRtime(index))
  svalue(specInfo$ind) <- paste(counter, " of ", length(currSequence))
  svalue(specInfo$acno) <- spIndex(index)
  svalue(specInfo$mslvl) <- spMsLevel(index)  
  svalue(specPrecInfo$pmz) <- round(spPrecMz(index), digits=settings$digits)
  svalue(specPrecInfo$int) <- round(spPrecInt(index), digits=settings$digits)
  svalue(specPrecInfo$charge) <- round(spPrecCharge(index), digits=settings$digits) 
}

updateSpectrum <- function(h=list(action=0), ...) {
  # If called by buttons Previous or Next, then h$action will have value -1 or 1. 
  env$counter <- counter + h$action
  env$index <- currSequence[counter]
  
  updateSpectrumInfo()
  
  # kill precursor info for MS1
  dispPrec <- spMsLevel(index)==2  
  lapply(specPrecInfo, function(i) enabled(i) <- dispPrec) 
  lapply(regularPrec, function(i) enabled(i) <- dispPrec)
  
  # Turn off buttons, click handlers and filters while drawing
  buttonSwitch(FALSE)
  clickSwitch(FALSE)  
  
  plotSpectrum(spectrumZoom)
  
  if(!zoomWindowClosed) {
    visible(env$plotZoom) <- TRUE      
    plotSpectrumZoom(spectrumZoom)
  }
  
  if(env$anyMS1spectra) {
    plotXIC(XICZoom) 
    if(!XICWindowClosed) {
      visible(env$plotXICw) <- TRUE      
      plotChromaZoom()
    }
  }
  
  clickSwitch(TRUE)
  buttonSwitch(TRUE)
}

clickSwitch <- function(on) {
  if(on) {
    mapply(unblockHandler, list(plotTop, plotBottom), clickHandlerIDs) 
#     if(verbose) cat("Plot interactivity activated\n")
  } else {
    mapply(blockHandler, list(plotTop, plotBottom), clickHandlerIDs) 
#     if(verbose) cat("Plot interactivity disabled\n")   
  }
  if(!env$anyMS1spectra) {
    blockHandler(plotBottom, clickHandlerIDs[[2]])
  }
}

buttonSwitch <- function(action=TRUE) {
  # (Selectively) turn on Next/Previous buttons
  if(action) {    
    # If current spectra selection has a single spectrum, both buttons disabled
    if(length(currSequence)==1) {
      enabled(buttonLeft) <- FALSE
      enabled(buttonRight) <- FALSE      
    } else 
      # If we have multiple spectra and the first is current
      if(counter==1) {
        enabled(buttonLeft) <- FALSE
        enabled(buttonRight) <- TRUE
      } else 
        # ... the last one is current...
        if(counter==length(currSequence)) {
          enabled(buttonLeft) <- TRUE
          enabled(buttonRight) <- FALSE
          # ...or we're in the middle
        } else {
          enabled(buttonLeft) <- TRUE
          enabled(buttonRight) <- TRUE
        }
  } else {
    
    enabled(buttonLeft) <- FALSE
    enabled(buttonRight) <- FALSE
    
  }
}

openFileHandler <- function(h, ...) {
  filename <- gfile(text = "Open file", type = "open", 
                    filter = settings$fileTypes, 
                    handler=function(h, ...) NULL)
  if(!is.na(filename)) {
    if(experimentLoaded) {      
      visible(plotTop) <- TRUE
      plotMsg()
      visible(plotBottom) <- TRUE
      plotMsg()
    }
    makeMzRrampAccessors(filename, env) #envir)
    updateExperiment(env)    
    experimentLoaded <- TRUE
  }
}

drawMain <- function(env) {
  
  # Window and structure
  env$msGUIWindow <- gwindow("msGUI", visible=FALSE, handler=function(h, ...) {
    env$closingMsGUI <- TRUE
    if(!zoomWindowClosed) dispose(zoomWindow)
    if(!XICWindowClosed) dispose(XICWindow)
    # Clean-up cache
    if(class(cache)!="function") {
      files <- unlist(cache$spectra[!sapply(cache$spectra, is.null)])
      if(length(files)>0) file.remove(files)
      files <- unlist(cache$xic[!sapply(cache$xic, is.null)])
      if(length(files)>0) file.remove(files)      
    }
    # Clean-up graphics devices    
    sapply(dev.list()[grepl("^png", names(dev.list()))], dev.off)
  })
  env$groupMain <- ggroup(container=msGUIWindow, horizontal=FALSE) 
  env$groupUpper <- ggroup(container=groupMain)
  env$groupMiddle <- ggroup(container=groupMain) 
  env$groupMiddleLeft <- ggroup(container=groupMiddle, horizontal=FALSE)
  env$groupMiddleRight <- ggroup(container=groupMiddle, horizontal=FALSE)
  env$groupPlots <- ggroup(container=groupMiddleRight, horizontal=FALSE)
  
  ## Top group
  env$buttonOpenFile <- gbutton("Open file", container=groupUpper, 
                                handler=openFileHandler)
  env$buttonOpenObject <- gbutton("Open R object", container=groupUpper, 
                                  handler=drawVarBrowser)
  addSpring(groupUpper)
  env$buttonSettings <- gbutton("Settings", container=groupUpper, 
                                handler=drawOptions)
  env$buttonHelp <- gbutton("Help", container=groupUpper)
  
  # Lists for formatting and retrieval
  env$headings <- list()
  env$regular <- list()
  env$regularPrec <- list() # to distinguish precursor-related text fields
  env$deemp <- list()
  env$separator <- list()
  env$expInfo <- list()
  env$specInfo <- list()
  env$specPrecInfo <- list() # to distinguish precursor-related edit fields
  env$filterInfo <- list()
  env$filterInfoMS <- list()
  env$filterInfoXIC <- list()
  
  # Left pane content
  
  env$le <- glayout(container=groupMiddleLeft, spacing=1)
  
  # Experiment info
  
  i <- 1
  
  le[i, 1:3] <- (env$separator$t1 <- glabel("", container=le))
  
  le[i + 1, 1:3] <- (env$headings$t1 <- glabel("Experiment information", container=le))
  le[i + 2, 2, anchor=c(-1,-1)] <- (env$deemp$t1 <- glabel("from", container=le))
  le[i + 2, 3, anchor=c(-1,-1)] <- (env$deemp$t2 <- glabel("to", container=le))
  le[i + 3, 1] <- (env$regular$t1 <- glabel("Retention time", container=le))
  le[i + 3, 2] <- (env$expInfo$rtfrom <- glabel(container=le))
  le[i + 3, 3] <- (env$expInfo$rtto <- glabel(container=le))
  le[i + 4, 1] <- (env$regular$t2 <- glabel("Precursor M/Z", container=le))
  le[i + 4, 2] <- (env$expInfo$pmzfrom <- glabel(container=le))
  le[i + 4, 3] <- (env$expInfo$pmzto <- glabel(container=le))
  le[i + 5, 1] <- (env$regular$t3 <- glabel("Spectra count", container=le))
  le[i + 5, 2] <- (env$regular$t4 <- glabel("MS1", container=le))
  le[i + 5, 3] <- (env$expInfo$ms1 <- glabel(container=le))
  le[i + 6, 2] <- (env$regular$t15 <- glabel("MS2", container=le))
  le[i + 6, 3] <- (env$expInfo$ms2 <- glabel(container=le))
  
  le[i + 7, 1:3] <- (env$separator$t2 <- glabel("", container=le))
  
  # Spectrum Info
  
  i <- 9
  
  le[i, 1:3] <- (env$headings$t2 <- glabel("Spectrum information", container=le))
  le[i + 1, 1] <- (env$regular$t5 <- glabel("Retention time", container=le))
  le[i + 1, 2:3] <- (env$specInfo$rt <- glabel("", container=le))
  le[i + 2, 1] <- (env$regular$t6 <- glabel("Index", container=le))
  le[i + 2, 2:3] <- (env$specInfo$ind <- glabel("", container=le))
  le[i + 3, 1] <- (env$regular$t16 <- glabel("Acquisition number", container=le))
  le[i + 3, 2:3] <- (env$specInfo$acno <- glabel("", container=le))
  le[i + 4, 1] <- (env$regular$t17 <- glabel("MS level", container=le))
  le[i + 4, 2:3] <- (env$specInfo$mslvl <- glabel("", container=le))
  le[i + 5, 1] <- (env$regularPrec$t1 <- glabel("Precursor M/Z", container=le))
  le[i + 5, 2:3] <- (env$specPrecInfo$pmz <- glabel("", container=le))
  le[i + 6, 1] <- (env$regularPrec$t2 <- glabel("Precursor intensity", container=le))
  le[i + 6, 2:3] <- (env$specPrecInfo$int <- glabel("", container=le))
  le[i + 7, 1] <- (env$regularPrec$t3 <- glabel("Precursor charge", container=le))
  le[i + 7, 2:3] <- (env$specPrecInfo$charge <- glabel("", container=le))
  
  le[i + 8, 1:5] <- (env$separator$t3 <- glabel("", container=le))
  
  # Filters
  
  i <- 18
  
  le[i, 1] <- (env$headings$t3 <- glabel("Filter", container=le))
  le[i, 2, anchor=c(-1,-1)] <- (env$deemp$t5 <- glabel("from", container=le))
  le[i, 3, anchor=c(-1,-1)] <- (env$deemp$t6 <- glabel("to", container=le))
  
  le[i + 1, 1] <- (env$filterInfo$rt$active <- gcheckbox("Retention time", container=le))
  le[i + 1, 2] <- (env$filterInfo$rt$from <- gedit("", container=le, width=5))
  le[i + 1, 3] <- (env$filterInfo$rt$to <- gedit("", container=le, width=5))
  le[i + 2, 1] <- (env$filterInfo$index$active <- gcheckbox("Acquisition number", container=le))
  le[i + 2, 2] <- (env$filterInfo$index$from <- gedit("", container=le, width=5))
  le[i + 2, 3] <- (env$filterInfo$index$to <- gedit("", container=le, width=5))
  le[i + 3, 1] <- (env$filterInfo$pmz$active <- gcheckbox("Precursor M/Z", container=le))
  le[i + 3, 2] <- (env$filterInfo$pmz$from <- gedit("", container=le, width=5))
  le[i + 3, 3] <- (env$filterInfo$pmz$to <- gedit("", container=le, width=5))
  le[i + 4, 1] <- (env$filterInfo$spi$active <- gcheckbox("Precursor intensity", container=le))
  le[i + 4, 2] <- (env$filterInfo$spi$from <- gedit("", container=le, width=5))
  le[i + 4, 3] <- (env$filterInfo$spi$to <- gedit("", container=le, width=5))
  le[i + 5, 1] <- (env$filterInfo$pc$active <- gcheckbox("Precursor charge", container=le))
  le[i + 5, 2] <- (env$filterInfo$pc$from <- gedit("", container=le, width=5))
  le[i + 5, 3] <- (env$filterInfo$pc$to <- gedit("", container=le, width=5))
  le[i + 6, 1] <- (env$filterInfo$mass$active <- gcheckbox("Precursor mass", container=le))
  le[i + 6, 2] <- (env$filterInfo$mass$from <- gedit("", container=le, width=5))
  le[i + 6, 3] <- (env$filterInfo$mass$to <- gedit("", container=le, width=5))
  
  le[i + 7, 1:5] <- (env$separator$t3 <- glabel("", container=le))
  
  le[i + 8, 1, anchor=c(-1,0)] <- (env$regular$t12 <- glabel("Display MS levels", container=le))
  le[i + 8, 2] <- (env$filterInfoMS$ms1 <- gcheckbox("MS1", checked=TRUE, container=le))
  le[i + 8, 3] <- (env$filterInfoMS$ms2 <- gcheckbox("MS2", checked=TRUE, container=le))
  
  le[i + 9, 1:5] <- (env$separator$t9 <- glabel("", container=le))
  
  le[i + 10, 1] <- (env$filterInfoXIC$XIC$active <- gcheckbox("Prec M/Z for XIC", container=le))
  le[i + 10, 2] <- (env$filterInfoXIC$XIC$from <- gedit("", container=le, width=5))
  le[i + 10, 3] <- (env$filterInfoXIC$XIC$to <- gedit("", container=le, width=5))
  
  le[i + 11, 1:5] <- (env$separator$t19 <- glabel("", container=le))
  
  # Buttons
  
  addSpring(groupMiddleLeft)
  
  env$groupLeftButtons <- ggroup(container=groupMiddleLeft)
  env$buttonLeft <- gbutton(text=gettext("Previous"), handler=updateSpectrum, 
                            action=-1, cont=groupLeftButtons)
  
  env$buttonRight <- gbutton(text=gettext("Next"), handler=updateSpectrum, 
                             action=1, cont=env$groupLeftButtons)
  
  env$plotTop <- ggraphics(container=groupPlots, width=settings$width,
                           height=settings$spectrumHeight, ps=12, dpi=75)
  env$plotBottom <- ggraphics(container=groupPlots, width=settings$width, 
                              height=settings$chromaHeight, ps=12, dpi=75)
  
  # Styling
  
  lapply(headings, setFont, settings$fontHead)
  lapply(regular, setFont, settings$fontReg)
  lapply(regularPrec, setFont, settings$fontReg)
  lapply(deemp, setFont, settings$fontDeemp)
  lapply(expInfo, setFont, settings$fontReg)
  lapply(specInfo, setFont, settings$fontReg)
  lapply(separator, setFont, list(size=4))
  lapply(specPrecInfo, setFont, settings$fontReg)
  lapply(filterInfo, function(x) lapply(x, setFont, settings$fontReg)) 
  lapply(filterInfoXIC, function(x) lapply(x, setFont, settings$fontReg)) 
  
  # Zoom handlers and GUI functions
  
  drawZoom <- function() {    
    env$zoomWindowClosed <- FALSE
    env$zoomWindow <- gwindow("Spectrum fragment", visible=FALSE, 
                              handler=function(h, ...) {
                                env$zoomWindowClosed <- TRUE
                                env$spectrumZoom <- NULL
                                if(!closingMsGUI) plotSpectrum()
                              })
    env$groupZoomMain <- ggroup(container=zoomWindow, horizontal=FALSE, expand=TRUE)
    env$plotZoom <- ggraphics(width=250, height=250, container=groupZoomMain, dpi=96, ps=12)
    visible(zoomWindow) <- TRUE
    env$clickHandlerZoomIDs <- addHandlerChanged(env$plotZoom, handler=handlerClickZoom)
  }   
  
  drawZoomXIC <- function() {    
    env$XICWindowClosed <- FALSE
    env$XICWindow <- gwindow("XIC", visible=FALSE, 
                             handler=function(h, ...) {
                               env$XICWindowClosed <- TRUE
                               env$XICZoom <- NULL
                               if(!closingMsGUI) plotXIC()
                             })
    env$groupXICMain <- ggroup(container=XICWindow, horizontal=FALSE, expand=TRUE)
    env$plotXICw <- ggraphics(width=250, height=250, container=groupXICMain, dpi=96, ps=12)
    visible(XICWindow) <- TRUE
    env$clickHandlerZoomXICIDs <- addHandlerChanged(env$plotXICw, handler=handlerClickZoomXIC)
  } 
  
  fixX <- function(x, lower, upper) {
    if(device=="png") {
      x <- ((x + .04) / 1.08 - 50/settings$width) * settings$width / (settings$width - 50 - 25) * (upper - lower) + lower
      sapply(sapply(x, max, lower), min, upper) # so that coordinates don't exceed limits
    } else x
  }
  
  fixY <- function(x, lower, upper, height) {
    if(device=="png") {
      x <- ((x + .04) / 1.08 - 40/height) * height / (height - 40 - 22) * (upper - lower) + lower
      sapply(sapply(x, max, lower), min, upper)
    } else x
  }  
  
  handlerClickSpectrum = function(h,...) {
    env$spectrumZoom <- list(x=fixX(h$x, 
                                    xLimits[, spMsLevel(index)][1], 
                                    xLimits[, spMsLevel(index)][2]),
                             y=fixY(h$y, 0, 1.05, settings$spectrumHeight))
                             
    if(verbose) cat("coords: x", h$x, "y", h$y, 
                    "  recalculated coords: x", spectrumZoom$x, 
                    "y", spectrumZoom$y, "\n") 
    
    plotSpectrum(zoom=spectrumZoom)
    if(zoomWindowClosed) drawZoom()
    visible(plotZoom) <- TRUE
    plotSpectrumZoom(spectrumZoom)
  }   
  
  handlerClickZoom <- function(h,...) {
    if(verbose) cat("coords: x", h$x, "y", h$y, "\n")
    env$spectrumZoom <- h
    plotSpectrum(zoom=h)
    if(zoomWindowClosed) drawZoom()
    visible(env$plotZoom) <- TRUE          
    plotSpectrumZoom(h)
  } 
  
  handlerClickZoomXIC <- function(h, ...) {
    if(verbose) cat("coords: x", h$x, "y", h$y, "\n")
    env$XICZoom <- h    
    plotXIC(zoom=h)    
    if(XICWindowClosed) drawZoomXIC(env)
    visible(env$plotXICw) <- TRUE          
    plotChromaZoom()
  }
  
  handlerClickXIC <- function(h,...) {
    clickSwitch(FALSE)
    xicRangeX <- range(xic(n=1, FALSE)[, 1])
    xCoord <- fixX(h$x, xicRangeX[1], xicRangeX[2])[2]
    
    if(verbose) cat("coords: x", h$x, "y", h$y, "recalculated coords", xCoord, "\n")
    
    if(same(h)) {
      prevCounter <- counter    
      env$counter <- which.min(abs(spRtime(currSequence)-xCoord))
      
      # update graphs if index has changed
      if(prevCounter!=counter) updateSpectrum() 
      
    } else {
      env$XICZoom$x <- fixX(h$x, xicRangeX[1], xicRangeX[2])
      env$XICZoom$y <- fixY(h$y, 0, 1.05, settings$chromaHeight) 
      
      if(XICWindowClosed) drawZoomXIC()
      visible(plotXICw) <- TRUE
      plotChromaZoom()
      plotXIC(XICZoom)
    }
    
    
    clickSwitch(TRUE)   
  }
  
  env$clickHandlerIDs <- list()
  env$clickHandlerIDs[[1]] <- addHandlerChanged(env$plotTop, handler=handlerClickSpectrum)
  env$clickHandlerIDs[[2]] <- addHandlerChanged(env$plotBottom, handler=handlerClickXIC)
  
  env$filterSpectraHandlerIDs <- lapply(env$filterInfo, 
                                        function(i) lapply(i, addHandlerChanged, 
                                                           handler=filterSpectra))
  
  env$filterSpectraMSHandlerIDs <- lapply(env$filterInfoMS, addHandlerChanged, 
                                          handler=filterSpectra)
  env$filterXICHandlerIDs <- lapply(env$filterInfoXIC, 
                                    function(i) lapply(i, addHandlerChanged, 
                                                       handler=filterSpectra))
}

getObjects <- function(classes="All classes") {
  objects <- ls(envir=globalenv())
  x <- data.frame(Object=objects, 
                  Class=sapply(objects, function(x) class(get(x))), 
                  stringsAsFactors=FALSE)
  if(classes=="All classes") 
    rows <- x$Class!="function" else rows <- x$Class==classes
  if(sum(rows)==0) return(data.frame(Object="No objects found", 
                                     Comment="", 
                                     stringsAsFactors=FALSE))
  x <- x[rows, -2, drop=FALSE]
  x$Comment <- paste(sapply(x$Object, function(i) length(get(i))), "unique MZs")
  return(x)
}

openObject <- function(object) {
  if(verbose) cat("pening...", object)
  if(experimentLoaded) {      
    visible(plotTop) <- TRUE
    plotMsg()
    visible(plotBottom) <- TRUE
    plotMsg()
  }
  makeMSnExpAccessors(get(object), env)
  if(verbose) cat("     done!\n")
  updateExperiment(env)
}

drawVarBrowser <- function(h, ...) {
  windowVB <- gwindow(title="Browse R objects", visible=FALSE, 
                      width=400, height=300, parent=msGUIWindow)
  panelVB <- ggroup(container=windowVB, horizontal=FALSE, expand=TRUE)
  panelVBtop <- ggroup(container=panelVB)
  panelVBmiddle <- ggroup(container=panelVB, expand=TRUE)
  panelVBbottom <- ggroup(container=panelVB) 
  text <- glabel("Filter objects by class: ", container=panelVBtop)
  filterVB <- gcombobox(container=panelVBtop, 
                        handler=function(h, ...) 
                          tableVB[] <- getObjects(svalue(filterVB)), 
                        items=c("MSnExp"))
  
  tableVB <- gtable(container=panelVBmiddle, items=getObjects(svalue(filterVB)), expand=TRUE)
  
  openObjectHandler <- function(h, ...) {
    if(length(svalue(tableVB))==0) 
      gmessage(title=" ", message="Please select an object!")
    else if(svalue(tableVB)=="No objects found") return(FALSE)
    else {
      openObject(svalue(tableVB))
      dispose(windowVB)
    }
  }
  
  addHandlerDoubleclick(tableVB, handler=openObjectHandler)
  
  buttonRefresh <- gbutton(text="Refresh", container=panelVBbottom, 
                           handler=function(h, ...) tableVB[] <- getObjects(svalue(filterVB)))
  addSpring(panelVBbottom)
  buttonOpen <- gbutton(text="Open", container=panelVBbottom, 
                        handler=openObjectHandler)
  visible(windowVB) <- TRUE
}

optsHandlerDefaults <- function(h, ...) {
  defaults <- defaultSettings()
  svalue(opts$spectrumHeight) <- defaults$spectrumHeight
  svalue(opts$chromaHeight) <- defaults$chromaHeight
  svalue(opts$width) <- defaults$width
  svalue(opts$labelNumber) <- defaults$labelNumber
  svalue(opts$MS1PlotType) <- ifelse(defaults$MS1PlotType=="h", "", " ")
  svalue(opts$MS2PlotType) <- ifelse(defaults$MS2PlotType=="h", "", " ")
  svalue(opts$chromaMode) <- ifelse(defaults$chromaMode, "Total ion count", 
                                    "Base peak intensity")
}

drawOptions <- function (h, ...) {
  
  if(!optionsWindowClosed) return(NULL)
  
  env$optionsWindowClosed <- FALSE
  
  env$optsWindow <- gwindow("Options", visible=FALSE, height=50, width=50, 
                            parent=msGUIWindow, handler=function(h, ...) 
                              env$optionsWindowClosed <- TRUE)
  env$optsGroup <- ggroup(container=optsWindow, horizontal=FALSE)
  env$l <- glayout(container=optsGroup, homogeneous=TRUE, spacing=30)
  env$l[1, 1] <- (env$l1 <- glayout(container=l, spacing=2))
  env$l[1, 2] <- (env$l2 <- glayout(container=l, spacing=2))
  
  env$l1[1, 1:3] <- (env$opts$headings$t1 <- glabel("Graph sizes", container=l1))
  env$l1[2, 2, anchor=c(0, 0)] <- (env$opts$text$t1 <- glabel("height", container=l1))  
  env$l1[2, 3, anchor=c(0, 0)] <- (env$opts$text$t2 <- glabel("width", container=l1))
  env$l1[3, 2] <- (env$opts$spectrumHeight <- gspinbutton(from=250, to=500, by=10, 
                                                          value=settings$spectrumHeight, 
                                                          digits=0, container=l1))
  env$l1[3, 3] <- (env$opts$width <- gspinbutton(from=500, to=800, by=10, 
                                                 value=settings$width, digits=0, 
                                                 container=l1))
  env$l1[3, 1, anchor=c(-1, 0)] <- (env$opts$text$t1 <- glabel("Spectrum graph", 
                                                                 container=l1))
  env$l1[4, 1, anchor=c(-1, 0)] <- (env$opts$text$t2 <- glabel("Chromatogram", 
                                                                 container=l1))
  env$l1[4, 2] <- (env$opts$chromaHeight <- gspinbutton(from=250, to=500, by=10, 
                                                        value=settings$chromaHeight, 
                                                        digits=0, container=l1))  
  i <- 5
  
  env$l1[i, 1:3] <- (env$opts$separator$t1 <- glabel("", container=l1))
  env$l1[i + 1, 1:3] <- (env$opts$headings$t3 <- glabel("Labels", container=l1))
  env$l1[i + 2, 1:2] <- (env$opts$text$t3 <- glabel("Number of peaks to label", 
                                                    container=l1))
  env$l1[i + 2, 3] <- (env$opts$labelNumber <- gedit(settings$labelNumber, 
                                                     coerce.with=as.numeric, 
                                                     width=2, container=l1))
  
  env$l2[1, 1:3] <- (env$opts$headings$t2 <- glabel("MS mode", container=l2))
  env$l2[2, 2, anchor=c(1, 0)] <- (env$opts$text$t4 <- glabel("MS1     ", container=l2))
  env$l2[2, 3, anchor=c(1, 0)] <- (env$opts$text$t5 <- glabel("MS2     ", container=l2))
  env$l2[3:4, 2] <- (env$opts$MS1PlotType <- gradio(c("", " "), 
                                                    ifelse(settings$MS1PlotType=="h", 1, 2), 
                                                    container=l2))
  env$l2[3:4, 3] <- (env$opts$MS2PlotType <- gradio(c("", " "), 
                                                    ifelse(settings$MS2PlotType=="h", 1, 2), 
                                                    container=l2))
  env$l2[3, 1, anchor=c(-1, 0)] <- (env$opts$text$t6 <- glabel("Centroided             ", container=l2))
  env$l2[4, 1, anchor=c(-1, 0)] <- (env$opts$text$t7 <- glabel("Profile", container=l2))
  
  i <- 5
  
  env$l2[i, 1:3] <- (env$opts$separator$t2 <- glabel("", container=l2))
  env$l2[i + 1, 1:3] <- (env$opts$headings$t4 <- glabel("Chromatogram mode", 
                                                        container=l2))
  env$l2[i + 2, 1:3] <- (env$opts$chromaMode <- gradio(c("Total ion count", 
                                                               "Base peak intensity"), 
                                                       ifelse(settings$chromaMode, 2, 1), 
                                                       container=l2))
  env$l2[i + 3, 1:3] <- (env$opts$separator$t3 <- glabel("", container=l2))
  
  env$optsButtons <- ggroup(container=optsGroup, horizontal=TRUE)
  env$opts$defaults <- gbutton("Restore defaults", container=optsButtons, 
                               width=50, handler=optsHandlerDefaults)
  
  addSpring(env$optsButtons)
  
  env$opts$ok <- gbutton("OK", container=optsButtons, width=settings$minButtonWidth, handler=function(h, ...) {
    if(any(c(settings$spectrumHeight!=svalue(opts$spectrumHeight), 
             settings$chromaHeight!=svalue(opts$chromaHeight), 
             settings$width!=svalue(opts$width), 
             settings$labelNumber!=svalue(opts$labelNumber), 
             settings$MS1PlotType!=ifelse(svalue(opts$MS1PlotType)=="", "h", "l"), 
             settings$MS2PlotType!=ifelse(svalue(opts$MS2PlotType)=="", "h", "l"), 
             settings$chromaMode!=(svalue(opts$chromaMode)=="Base peak intensity")
             ))) {
      if(verbose) cat("Applying changes... ")
      settings$spectrumHeight <- svalue(opts$spectrumHeight)
      settings$chromaHeight <- svalue(opts$chromaHeight)
      settings$width <- svalue(opts$width)
      settings$labelNumber <- svalue(opts$labelNumber)
      settings$MS1PlotType <- ifelse(svalue(opts$MS1PlotType)=="", "h", "l")
      settings$MS2PlotType <- ifelse(svalue(opts$MS2PlotType)=="", "h", "l")
      settings$chromaMode <- (svalue(opts$chromaMode)=="Base peak intensity")
      
      size(plotTop) <- c(settings$width, settings$spectrumHeight)
      size(plotBottom) <- c(settings$width, settings$chromaHeight)
      if(!zoomWindowClosed) {
        visible(env$plotZoom) <- TRUE      
        plotSpectrumZoom(spectrumZoom)        
      }
      if(!XICWindowClosed) {
        visible(env$plotXICw) <- TRUE      
        plotChromaZoom()        
      }        
      resetCache()
      updateSpectrum()
      if(verbose) cat("     done!\n")
    }
    dispose(env$optsWindow)
  })
  env$opts$cancel <- gbutton("Cancel", container=optsButtons, width=50, 
                             handler=function(h, ...) dispose(env$optsWindow))
  addSpace(env$optsButtons, 18)
  
  lapply(opts$headings, setFont, settings$fontHead)
  
  visible(env$optsWindow) <- TRUE
}

setFont <- function(x, style) font(x) <- style

same <- function(h) abs(h$x[1]-h$x[2]) < 1/200 & abs(h$y[1]-h$y[2]) < 1/200
