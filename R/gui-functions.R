wrapper <- function(filename=NULL, device="png", verbose=FALSE) {  
  
  env <- environment()
  
  settings <- defaultSettings()
  
  zoomWindowClosed <- TRUE
  XICWindowClosed <- TRUE
  spectrumZoom <- NULL
  
  environment(openFileHandler) <- env  
  environment(drawMain) <- env  
  environment(updateExperiment) <- env
  environment(updateExperimentInfo) <- env
  environment(updateSpectrum) <- env  
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
  
  drawMain(env)
  initialiseGUI()
  
  if(!is.null(filename)) {
    makeMzRrampAccessors(filename, env)
    updateExperiment(env)
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
  env$xLimits <- sapply(unique(spMsLevel()), 
                        function(i) c(min(spLowMZ()[spMsLevel()==i]), 
                                      max(spHighMZ()[spMsLevel()==i])))
  resetCache()
  counterReset(env)    
  filterReset(env)    
  filterSwitch(TRUE)
  updateSpectrum()  
}

updateSpectrum <- function(h=list(action=0), ...) {
  # If called by buttons Previous or Next, then h$action will have value -1 or 1. 
  env$counter <- counter + h$action
  env$index <- currSequence[counter]
  svalue(specInfo$rt) <- formatRt(spRtime(index))
  svalue(specInfo$ind) <- paste(counter, " of ", length(currSequence))
  svalue(specInfo$acno) <- spIndex(index)
  svalue(specInfo$mslvl) <- spMsLevel(index)  
  svalue(specPrecInfo$pmz) <- round(spPrecMz(index), digits=settings$digits)
  svalue(specPrecInfo$int) <- round(spPrecInt(index), digits=settings$digits)
  svalue(specPrecInfo$charge) <- round(spPrecCharge(index), digits=settings$digits) 
  
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
  plotXIC() 
  
  clickSwitch(TRUE)
  buttonSwitch(TRUE)
}

clickSwitch <- function(on) {
  if(on) {
    mapply(unblockHandler, list(plotBottom, plotTop), clickHandlerIDs) 
    if(verbose) cat("\nPlots interactive")
  } else {
    mapply(blockHandler, list(plotBottom, plotTop), clickHandlerIDs) 
    if(verbose) cat("\nPlot interactivity disabled")   
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
#     envir <- parent.env(environment())
    makeMzRrampAccessors(filename, env) #envir)
    updateExperiment(env) #ir)      
  }
}

drawMain <- function(env) {
  
  # Window and structure
  env$msGUIWindow <- gwindow("msGUI", visible=FALSE, handler=function(h, ...) {
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
  env$buttonSettings <- gbutton("Settings", container=groupUpper)
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
  
  le[i + 10, 1] <- (env$filterInfoXIC$XIC$active <- gcheckbox("XIC", container=le))
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
  
  env$plotTop <- ggraphics(container=env$groupPlots, width=500, height=250, ps=12, dpi=75)
  env$plotBottom <- ggraphics(container=env$groupPlots, width=500, height=250, ps=12, dpi=75)
  
  # Styling
  
  setFont <- function(x, style) font(x) <- style
  
  lapply(headings, setFont, settings$fontHead)
  lapply(regular, setFont, settings$fontReg)
  lapply(regularPrec, setFont, settings$fontReg)
  lapply(deemp, setFont, settings$fontDeemp)
  lapply(expInfo, setFont, settings$fontReg)
  lapply(specInfo, setFont, settings$fontReg)
  lapply(separator, setFont, list(size=4))
  lapply(filterInfoMS, setFont, settings$fontReg)
  lapply(filterInfo, function(x) lapply(x, setFont, settings$fontReg)) 
  lapply(filterInfoXIC, function(x) lapply(x, setFont, settings$fontReg)) 
  
  # Zoom handlers and GUI functions
  
  drawZoom <- function(env) {    
    env$zoomWindowClosed <- FALSE
    env$zoomWindow <- gwindow("Spectrum fragment", visible=FALSE, 
                              handler=function(h, ...) {
                                env$zoomWindowClosed <- TRUE
                                env$spectrumZoom <- NULL
                                plotSpectrum()
                              })
    env$groupZoomMain <- ggroup(container=env$zoomWindow, horizontal=FALSE, expand=TRUE)
    env$plotZoom <- ggraphics(width=250, height=250, container=env$groupZoomMain, dpi=96, ps=12)
    visible(zoomWindow) <- TRUE
    env$clickHandlerZoomIDs <- addHandlerChanged(env$plotZoom, handler=handlerClickZoom)
  }   
  
  drawZoomXIC <- function() {    
    env$XICWindowClosed <- FALSE
    env$XICWindow <- gwindow("XIC", visible=FALSE, 
                             handler=function(h, ...) {
                               env$XICWindowClosed <- TRUE
                               env$XICZoom <- NULL
                               plotXIC()
                             })
    env$groupXICMain <- ggroup(container=env$XICWindow, horizontal=FALSE, expand=TRUE)
    env$plotXICw <- ggraphics(width=250, height=250, container=env$groupXICMain, dpi=96, ps=12)
    visible(XICWindow) <- TRUE
#     env$clickHandlerZoomXICIDs <- addHandlerChanged(env$plotZoomXIC, handler=handlerClickZoomXIC)
  } 
  
  fixX <- function(x, lower, upper) {
    if(device=="png") {
      x <- ((x + .04) / 1.08 - 50/500) * 500 / (500 - 50 - 25) * (upper - lower) + lower
      #       x <- ((x + .04) / 1.08 - 40/500) * 500 / (500 - 40 - 20) * (upper - lower) + lower
      sapply(sapply(x, max, lower), min, upper) # so that coordinates don't exceed limits
    } else x
  }
  
  fixY <- function(x, lower, upper) {
    if(device=="png") {
      x <- ((x + .04) / 1.08 - 40/250) * 250 / (250 - 40 - 22) * (upper - lower) + lower
      sapply(sapply(x, max, lower), min, upper)
    } else x
  }  
  
  handlerClickSpectrum = function(h,...) {
    env$spectrumZoom <- list(x=fixX(h$x, 
                                    min(spLowMZ()[spMsLevel()==spMsLevel(index)]), 
                                    max(spHighMZ()[spMsLevel()==spMsLevel(index)])), 
                             y=fixY(h$y, 0, 1))
    
    if(verbose) cat("\nx:", h$x, "y: ", h$y, 
                    "  fixed: x:", spectrumZoom$x, 
                    "y:", spectrumZoom$y) 
    
    plotSpectrum(zoom=spectrumZoom)
    if(zoomWindowClosed) drawZoom(env)
    visible(plotZoom) <- TRUE
    plotSpectrumZoom(spectrumZoom)
  }   
  
  handlerClickZoom = function(h,...) {
    if(verbose) cat("\nx:", h$x, "y: ", h$y)
    env$spectrumZoom <- h
    plotSpectrum(zoom=h)
    if(zoomWindowClosed) drawZoom(env)
    visible(env$plotZoom) <- TRUE          
    plotSpectrumZoom(h)
  } 
  
  handlerClickXIC = function(h,...) {
    clickSwitch(FALSE)
    x <- xic(n=1, FALSE)
    xicRangeX <- range(x[, 1])
#     xicRangeY <- range(x[, 2])
    env$XICZoom$x <- fixX(h$x, xicRangeX[1], xicRangeX[2])
    env$XICZoom$y <- fixY(h$y, 0, 1) #xicRangeY[1], xicRangeY[2])
    
    if(verbose) cat("\nXIC clicked on:", c(h$x, h$y), " fixed x:", XICZoom$x) 
    
    if(h$x[1]==h$x[2] & h$y[1]==h$y[2]) {
      prevCounter <- counter    
      env$counter <- which.min(abs(spRtime(currSequence)-XICZoom$x[1]))
      
      # update graphs if index has changed
      if(prevCounter!=counter) updateSpectrum() 
      
      if(verbose) cat("displaying spectrum", currSequence[counter])
    } else {
      if(verbose) cat("displaying xic zoom")
      if(XICWindowClosed) drawZoomXIC()
      visible(plotXICw) <- TRUE
      plotChromaZoom()
      plotXIC(XICZoom)
    }
    
    
    clickSwitch(TRUE)   
  }
  
  env$clickHandlerIDs <- list()
  env$clickHandlerIDs[[1]] <- addHandlerChanged(env$plotBottom, handler=handlerClickXIC)
  env$clickHandlerIDs[[2]] <- addHandlerChanged(env$plotTop, handler=handlerClickSpectrum)
  
  env$filterSpectraHandlerIDs <- lapply(env$filterInfo, 
                                        function(i) lapply(i, addHandlerChanged, 
                                                           handler=filterSpectra))
  #   env$filterSpectraHandlerIDs <- c(filterSpectraHandlerIDs, 
  #                                    lapply(env$filterInfo, 
  #                                         function(i) lapply(i, addHandlerBlur, 
  #                                                            handler=filterSpectra)))
  #   # error:
  #   (rsession.exe:7992): Gtk-WARNING **: GtkEntry - did not receive focus-out-event. If you
  #   connect a handler to this signal, it must return
  #   FALSE so the entry gets the event as well
  
  env$filterSpectraMSHandlerIDs <- lapply(env$filterInfoMS, addHandlerChanged, 
                                          handler=filterSpectra)
  env$filterXICHandlerIDs <- lapply(env$filterInfoXIC, 
                                    function(i) lapply(i, addHandlerChanged, 
                                                       handler=filterXIC))
}

getObjects <- function(classes="All classes") {
  objects <- ls(envir=globalenv())
  x <- cbind(Object=objects, 
             Class=sapply(objects, function(x) class(get(x))), 
             Comment=sapply(objects, function(x) ifelse(class(get(x))=="MSnExp", "Some attribute", "")))
  rows <- x[, 2]!="function"
  if(classes!="All classes") rows <- x[, 2]==classes
  if(sum(rows)==0) return(data.frame(Object="No objects found", 
                                     Class="", 
                                     Comment="", 
                                     stringsAsFactors=FALSE))
  x[rows, , drop=FALSE]
}

openObject <- function(object) {
  cat("\nYay, will open", object)
  makeMSnExpAccessors(get(object), env)
  updateExperiment(env)
}

drawVarBrowser <- function(h, ...) {
  windowVB <- gwindow(title="Browse R objects", visible=FALSE, 
                      width=400, height=300, parent=msGUIWindow)
  panelVB <- ggroup(container=windowVB, horizontal=FALSE, expand=TRUE)
  panelVBtop <- ggroup(container=panelVB)
  panelVBmiddle <- ggroup(container=panelVB, expand=TRUE)
  panelVBbottom <- ggroup(container=panelVB) #, horizontal=TRUE
  text <- glabel("Filter objects by class: ", container=panelVBtop)
  filterVB <- gcombobox(container=panelVBtop, 
                        handler=function(h, ...) 
                          tableVB[] <- getObjects(svalue(filterVB)), 
                        items=c("MSnExp", "All classes"))
  #   
  #   text <- glabel("Filter objects by class: ", container=panelVBtop)
  #   filterVB <- gcombobox(container=panelVBtop, 
  #                         handler=function(h, ...) 
  #                           tableVB[] <- getObjects(svalue(filterVB)), 
  #                         items=c("MSnExp", "All classes"))
  #   layoutVB <- glayout(spacing=2, container=panelVBtop)
  #   layoutVB[1, 1] <- glabel("Filter objects by class: ", container=layoutVB)
  #   layoutVB[1, 2] <- (filterVB <- gcombobox(container=layoutVB, 
  #                         handler=function(h, ...) 
  #                           tableVB[] <- getObjects(svalue(filterVB)), 
  #                         items=c("MSnExp", "All classes")))
  #   
  #   layoutVB[2, 1] <- glabel("Filter objects by class: ", container=layoutVB)
  #   layoutVB[2, 2] <- (filter2VB <- gcombobox(container=layoutVB, 
  #                         handler=function(h, ...) 
  #                           tableVB[] <- getObjects(svalue(filter2VB)), 
  #                         items=c("MSnExp", "All classes")))
  
  
  tableVB <- gtable(container=panelVBmiddle, items=getObjects(svalue(filterVB)), expand=TRUE)
  
  openObjectHandler <- function(h, ...) {
    if(length(svalue(tableVB))==0) 
      gmessage(title=" ", message="Please select object!")
    else if(svalue(tableVB)=="No objects found") return(FALSE)
    else openObject(svalue(tableVB))
  }
  
  addHandlerDoubleclick(tableVB, handler=openObjectHandler)
  
  buttonRefresh <- gbutton(text="Refresh", container=panelVBbottom, 
                           handler=function(h, ...) tableVB[] <- getObjects(svalue(filterVB)))
  addSpring(panelVBbottom)
  buttonOpen <- gbutton(text="Open", container=panelVBbottom, 
                        handler=openObjectHandler)
  visible(windowVB) <- TRUE
}