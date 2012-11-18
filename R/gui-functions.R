wrapper <- function(filename=NULL, device="png") {  
  
  env <- environment()
  
  initialiseEnvironment(env)
  zoomWindowClosed <- TRUE
  settings <- list()
  settings$allowMultipleWindows <- FALSE
  
  # Later allow user to save settings via `fix()`
  
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
  
  drawMain(env)
  initialiseGUI()
  
  if(!is.null(filename)) {
    makeMzRrampAccessors(filename, env)
    updateExperiment(env)
  }   
}

updateExperimentInfo <- function() {
  rtRange <- round(expRtRange(), digits=digits)
  precMzRange <- round(expPrecMzRange(), digits=digits)
  svalue(expInfo$rtfrom) <- rtRange[1]
  svalue(expInfo$rtto) <- rtRange[2]
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
  counterReset(env)    
  filterReset(env)    
  filterSwitch(TRUE)
  updateSpectrum()
  
}

updateSpectrum <- function(h=list(action=0), ...) {
  # If called by buttons Previous or Next, then h$action will have value -1 or 1. 
  env$counter <- counter + h$action
  env$index <- currSequence[counter]
  svalue(specInfo$rt) <- spRtime(index)
  svalue(specInfo$ind) <- paste(spIndex(index), "of", length(currSequence))
  svalue(specPrecInfo$pmz) <- round(spPrecMz(index), digits=digits)
  svalue(specPrecInfo$int) <- round(spPrecInt(index), digits=digits)
  svalue(specPrecInfo$charge) <- round(spPrecCharge(index), digits=digits)  
  
  # kill precursor info for MS1
  dispPrec <- spMsLevel(index)==2  
  lapply(specPrecInfo, function(i) enabled(i) <- dispPrec)
  lapply(regularPrec, function(i) enabled(i) <- dispPrec)
  
  # Turn off buttons, click handlers and filters while drawing
  buttonSwitch(FALSE)
#   filterSwitch(FALSE)
  clickSwitch(FALSE)  
  
  plotSpectrum()
  plotXIC() 
  
  clickSwitch(TRUE)
  buttonSwitch(TRUE)
#   filterSwitch(TRUE)
}

clickSwitch <- function(on) {
  if(on) {
    mapply(unblockHandler, list(plotBottom, plotTop), clickHandlerIDs) 
    cat("\nClicking enabled")
  } else {
    mapply(blockHandler, list(plotBottom, plotTop), clickHandlerIDs) 
    cat("\nClicking disabled")   
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
                    filter = file.types, handler=function(h, ...) NULL)
  if(!is.na(filename)) {
    envir <- parent.env(environment())
    makeMzRrampAccessors(filename, envir)
    updateExperiment(envir)      
  }
}

drawMain <- function(env) {
  
  if(device=="tkrplot") {
    options(guiToolkit = "tcltk")
    require(tkrplot)
  }
  
  # Window and structure
  env$msGUIWindow <- gwindow("msGUI", visible=FALSE)
  env$groupMain <- ggroup(container=env$msGUIWindow, horizontal=FALSE, expand=TRUE)
  env$groupUpper <- ggroup(container=env$groupMain)
  env$groupMiddle <- ggroup(container=env$groupMain, expand=TRUE)
  env$groupMiddleLeft <- ggroup(container=env$groupMiddle, horizontal=FALSE, expand=TRUE)
  env$groupMiddleRight <- gframe(container=env$groupMiddle, horizontal=FALSE, expand=TRUE)
  env$groupPlots <- ggroup(container=env$groupMiddleRight, horizontal=FALSE, expand=TRUE)
  
  ## Top group
  env$buttonOpenFile <- gbutton("Open file", container=env$groupUpper, 
                                handler=openFileHandler)
  env$buttonOpenObject <- gbutton("Open R object", container=env$groupUpper)
  addSpring(env$groupUpper)
  env$buttonSettings <- gbutton("Settings", container=env$groupUpper)
  env$buttonHelp <- gbutton("Help", container=env$groupUpper)
  
  # Lists for formatting and retrieval
  env$headings <- list()
  env$regular <- list()
  env$regularPrec <- list() # to distinguish precursor-related fields
  env$deemp <- list()
  env$separator <- list()
  env$expInfo <- list()
  env$specInfo <- list()
  env$specPrecInfo <- list() # to distinguish precursor-related fields
  env$filterInfo <- list()
  env$filterInfoMS <- list()
  
  # Left pane content
  
  env$le <- glayout(container=env$groupMiddleLeft, spacing=2)
  
  # Experiment info
  
  i <- 1
  
  env$le[i, 1:3] <- (env$separator$t1 <- glabel("", container=env$le))
  
  env$le[i + 1, 1:3] <- (env$headings$t1 <- glabel("Experiment information", container=env$le))
  env$le[i + 2, 2, anchor=c(-1,-1)] <- (env$deemp$t1 <- glabel("from", container=env$le))
  env$le[i + 2, 3, anchor=c(-1,-1)] <- (env$deemp$t2 <- glabel("to", container=env$le))
  env$le[i + 3, 1] <- (env$regular$t1 <- glabel("Retention time", container=env$le))
  env$le[i + 3, 2] <- (env$expInfo$rtfrom <- glabel(container=env$le))
  env$le[i + 3, 3] <- (env$expInfo$rtto <- glabel(container=env$le))
  env$le[i + 4, 1] <- (env$regular$t2 <- glabel("Precursor M/Z", container=env$le))
  env$le[i + 4, 2] <- (env$expInfo$pmzfrom <- glabel(container=env$le))
  env$le[i + 4, 3] <- (env$expInfo$pmzto <- glabel(container=env$le))
  env$le[i + 5, 1] <- (env$regular$t3 <- glabel("Spectra count", container=env$le))
  env$le[i + 5, 2] <- (env$regular$t4 <- glabel("MS1", container=env$le))
  env$le[i + 5, 3] <- (env$expInfo$ms1 <- glabel(container=env$le))
  env$le[i + 6, 2] <- (env$regular$t15 <- glabel("MS2", container=env$le))
  env$le[i + 6, 3] <- (env$expInfo$ms2 <- glabel(container=env$le))
  
  env$le[i + 7, 1:3] <- (env$separator$t2 <- glabel("", container=env$le))
  
  # Spectrum Info
  
  i <- 9
  
  env$le[i, 1:3] <- (env$headings$t2 <- glabel("Spectrum information", container=env$le))
  env$le[i + 1, 1] <- (env$regular$t5 <- glabel("Retention time", container=env$le))
  env$le[i + 1, 2:3] <- (env$specInfo$rt <- glabel("", container=env$le))
  env$le[i + 2, 1] <- (env$regular$t6 <- glabel("Index", container=env$le))
  env$le[i + 2, 2:3] <- (env$specInfo$ind <- glabel("", container=env$le, editable=TRUE))
  env$le[i + 3, 1] <- (env$regularPrec$t1 <- glabel("Precursor M/Z", container=env$le))
  env$le[i + 3, 2:3] <- (env$specPrecInfo$pmz <- glabel("", container=env$le))
  env$le[i + 4, 1] <- (env$regularPrec$t2 <- glabel("Precursor intensity", container=env$le))
  env$le[i + 4, 2:3] <- (env$specPrecInfo$int <- glabel("", container=env$le))
  env$le[i + 5, 1] <- (env$regularPrec$t3 <- glabel("Precursor charge", container=env$le))
  env$le[i + 5, 2:3] <- (env$specPrecInfo$charge <- glabel("", container=env$le))
  
  env$le[i + 6, 1:5] <- (env$separator$t3 <- glabel("", container=env$le))
  
  # Filters
  
  i <- 16
  
  env$le[i, 1] <- (env$headings$t3 <- glabel("Filter", container=env$le))
  env$le[i, 2, anchor=c(-1,-1)] <- (env$deemp$t5 <- glabel("from", container=env$le))
  env$le[i, 3, anchor=c(-1,-1)] <- (env$deemp$t6 <- glabel("to", container=env$le))
  
  env$le[i + 1, 1] <- (env$filterInfo$rt$active <- gcheckbox("Retention time", container=env$le))
  env$le[i + 1, 2] <- (env$filterInfo$rt$from <- gedit("", container=env$le, width=5))
  env$le[i + 1, 3] <- (env$filterInfo$rt$to <- gedit("", container=env$le, width=5))
  env$le[i + 2, 1] <- (env$filterInfo$index$active <- gcheckbox("Index", container=env$le))
  env$le[i + 2, 2] <- (env$filterInfo$index$from <- gedit("", container=env$le, width=5))
  env$le[i + 2, 3] <- (env$filterInfo$index$to <- gedit("", container=env$le, width=5))
  env$le[i + 3, 1] <- (env$filterInfo$pmz$active <- gcheckbox("Precursor M/Z", container=env$le))
  env$le[i + 3, 2] <- (env$filterInfo$pmz$from <- gedit("", container=env$le, width=5))
  env$le[i + 3, 3] <- (env$filterInfo$pmz$to <- gedit("", container=env$le, width=5))
  env$le[i + 4, 1] <- (env$filterInfo$spi$active <- gcheckbox("Precursor intensity", container=env$le))
  env$le[i + 4, 2] <- (env$filterInfo$spi$from <- gedit("", container=env$le, width=5))
  env$le[i + 4, 3] <- (env$filterInfo$spi$to <- gedit("", container=env$le, width=5))
  env$le[i + 5, 1] <- (env$filterInfo$pc$active <- gcheckbox("Precursor charge", container=env$le))
  env$le[i + 5, 2] <- (env$filterInfo$pc$from <- gedit("", container=env$le, width=5))
  env$le[i + 5, 3] <- (env$filterInfo$pc$to <- gedit("", container=env$le, width=5))
  env$le[i + 6, 1] <- (env$filterInfo$mass$active <- gcheckbox("Precursor mass", container=env$le))
  env$le[i + 6, 2] <- (env$filterInfo$mass$from <- gedit("", container=env$le, width=5))
  env$le[i + 6, 3] <- (env$filterInfo$mass$to <- gedit("", container=env$le, width=5))
    
  env$le[i + 7, 1:5] <- (env$separator$t3 <- glabel("", container=env$le))
  
  env$le[i + 8, 1, anchor=c(-1,0)] <- (env$regular$t12 <- glabel("Display MS levels", container=env$le))
  env$le[i + 8, 2] <- (env$filterInfoMS$ms1 <- gcheckbox("MS1", checked=TRUE, container=env$le))
  env$le[i + 8, 3] <- (env$filterInfoMS$ms2 <- gcheckbox("MS2", checked=TRUE, container=env$le))
  
  # Buttons
  
  env$groupLeftButtons <- ggroup(container=env$groupMiddleLeft, expand=FALSE, anchor=c(-1, 1))
  env$buttonLeft <- gbutton(text=gettext("Previous"), handler=env$updateSpectrum, 
                            action=-1, cont=env$groupLeftButtons)
  
  env$buttonRight <- gbutton(text=gettext("Next"), handler=env$updateSpectrum, 
                             action=1, cont=env$groupLeftButtons)
  
  if(any(c("cairo", "png")==device)) {
    
    env$plotTop <- ggraphics(container=env$groupPlots, width=500, height=250, ps=12, dpi=75)
    env$plotBottom <- ggraphics(container=env$groupPlots, width=500, height=250, ps=12, dpi=75)
    
  } else if(device=="gimage"){
  
    env$plotTop <- gimage(container=env$groupPlots)
    env$plotBottom <- gimage(container=env$groupPlots)
    
  } else if(device=="tkrplot") {
    
    plotTopDrawn <<- plotBottomDrawn <<- FALSE
    
  } else stop("Unimplemented device!")

  
  lapply(env$headings, function(x) font(x) <- env$textHead)
  lapply(env$regular, function(x) font(x) <- env$textReg)
  lapply(env$deemp, function(x) font(x) <- env$textDeemp)
  lapply(env$expInfo, function(x) font(x) <- env$textReg)
  lapply(env$specInfo, function(x) font(x) <- env$textReg)
#   lapply(env$filterInfo, function(x) font(x) <- env$textReg)
  lapply(env$separator, function(x) font(x) <- list(size=4))
    
  drawZoom <- function(env) {
    env$zoomWindow <- gwindow("Spectrum fragment", visible=TRUE, 
                              handler=function(h, ...) {
                                env$zoomWindowClosed <- TRUE
                                plotSpectrum()
                                })
    env$zoomWindowClosed <- FALSE
    env$groupZoomMain <- ggroup(container=env$zoomWindow, horizontal=FALSE, expand=TRUE)
    env$plotZoom <- ggraphics(width=600, height=600, container=env$groupZoomMain, dpi=75, ps=12)
  }  
  
  fixX <- function(x, lower, upper) {
    if(device=="png") {
      x <- ((x + .04) / 1.08 - 40/500) * 500 / (500 - 40 - 20) * (upper - lower) + lower
      sapply(sapply(x, max, lower), min, upper) # so that coordinates don't exceed limits
    } else x
  }
  
  fixY <- function(x, lower, upper) {
    if(device=="png") {
      x <- ((x + .04) / 1.08 - 40/250) * 250 / (250 - 40 - 22) * (upper - lower) + lower
      sapply(sapply(x, max, lower), min, upper)
    } else x
  }
  
  if(device!="tkrplot") {    
    handlerClickGeneric = function(h,...) {
      cat("\nclicked on:", c(h$x, h$y))
    } 
    
    handlerClickSpectrum = function(h,...) {
      spmin <- min(spLowMZ()[spMsLevel()==spMsLevel(index)])
      spmax <- max(spHighMZ()[spMsLevel()==spMsLevel(index)])
      cat("\nx:", h$x, "y: ", h$y, 
          "  fixed: x:", fixX(h$x, spmin, spmax), 
          "y:", fixY(h$y, 0, 1))
      
      plotSpectrum(zoom=list(x=fixX(h$x, spmin, spmax), 
                             y=fixY(h$y, 0, 1)))
      if(zoomWindowClosed | settings$allowMultipleWindows) drawZoom(env)
      visible(env$plotZoom) <- TRUE 
      plotSpectrumZoom(list(x=fixX(h$x, spmin, spmax), 
                            y=fixY(h$y, 0, 1)))
    }   
    
    handlerClickXIC = function(h,...) {
      xicrange <- range(xic(n=c(1, 2)[c(svalue(filterInfoMS$ms1), svalue(filterInfoMS$ms2))], 
                            FALSE)[, 1])
      cat("\nXIC clicked on:", c(h$x, h$y), "  fixed x:", fixX(h$x, 
                                                             xicrange[1], 
                                                             xicrange[2])) 
      
      prevCounter <- counter
      
      env$counter <- which.min(abs(spRtime(currSequence) - fixX(h$x, 
                                                             xicrange[1], 
                                                             xicrange[2])))
      
      # update graphs if index has changed
      if(prevCounter!=counter) updateSpectrum()      
      
    }
    env$clickHandlerIDs <- list()
    env$clickHandlerIDs[[1]] <- addHandlerClicked(env$plotBottom, handler=handlerClickXIC)
    env$clickHandlerIDs[[2]] <- addHandlerChanged(env$plotTop, handler=handlerClickSpectrum)
  }

  env$filterSpectraHandlerIDs <- lapply(env$filterInfo, 
                                        function(i) lapply(i, addHandlerChanged, 
                                                           handler=filterSpectra))
  env$filterSpectraMSHandlerIDs <- lapply(env$filterInfoMS, addHandlerChanged, 
                                                           handler=filterSpectra)
}
