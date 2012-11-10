wrapper <- function(filename=NULL, device="png") {  
  
  env <- environment()
  
  initialiseEnvironment(env)
  
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
  counter <<- counter + h$action
  index <<- currSequence[counter]
  svalue(specInfo$rt) <- spRtime(index)
  svalue(specInfo$ind) <- paste(spIndex(index), "of", length(currSequence))
  svalue(specInfo$pmz) <- round(spPrecMz(index), digits=digits)
  svalue(specInfo$int) <- round(spPrecInt(index), digits=digits)
  
  # Turn off buttons, click handlers and filters while drawing
  buttonSwitch(FALSE)
  filterSwitch(FALSE)
  clickSwitch(FALSE)  
  
  plotSpectrum()
  plotXIC() 
  
  clickSwitch(TRUE)
  buttonSwitch(TRUE)
  filterSwitch(TRUE)
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
  env$groupMiddleLeft <- ggroup(container=env$groupMiddle, horizontal=FALSE, spacing=0, expand=TRUE)
  env$groupMiddleRight <- gframe(container=env$groupMiddle, horizontal=FALSE, expand=TRUE)
  env$groupPlots <- gpanedgroup(container=env$groupMiddleRight, horizontal=FALSE, expand=TRUE)
  
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
  env$deemp <- list()
  env$separator <- list()
  env$expInfo <- list()
  env$specInfo <- list()
  env$filterInfo <- list()
  
  # Left pane content
  
  env$le <- glayout(container=env$groupMiddleLeft, spacing=2)
  
  # Experiment info
  
  env$le[1, 1:3] <- (env$separator$t1 <- glabel("", container=env$le))
  
  env$le[2, 1:3] <- (env$headings$t1 <- glabel("Experiment information", container=env$le))
  env$le[3, 2, anchor=c(-1,-1)] <- (env$deemp$t1 <- glabel("from", container=env$le))
  env$le[3, 3, anchor=c(-1,-1)] <- (env$deemp$t2 <- glabel("to", container=env$le))
  env$le[4, 1] <- (env$regular$t1 <- glabel("Retention time", container=env$le))
  env$le[4, 2] <- (env$expInfo$rtfrom <- glabel(container=env$le))
  env$le[4, 3] <- (env$expInfo$rtto <- glabel(container=env$le))
  env$le[5, 1] <- (env$regular$t2 <- glabel("Precursor M/Z", container=env$le))
  env$le[5, 2] <- (env$expInfo$pmzfrom <- glabel(container=env$le))
  env$le[5, 3] <- (env$expInfo$pmzto <- glabel(container=env$le))
  env$le[6, 1] <- (env$regular$t3 <- glabel("Spectra count", container=env$le))
  env$le[6, 2] <- (env$regular$t4 <- glabel("MS1", container=env$le))
  env$le[6, 3] <- (env$expInfo$ms1 <- glabel(container=env$le))
  env$le[7, 2] <- (env$regular$t15 <- glabel("MS2", container=env$le))
  env$le[7, 3] <- (env$expInfo$ms2 <- glabel(container=env$le))
  
  env$le[8, 1:3] <- (env$separator$t2 <- glabel("", container=env$le))
  
  # Spectrum Info
  
  env$le[9, 1:3] <- (env$headings$t2 <- glabel("Spectrum information", container=env$le))
  env$le[10, 1] <- (env$regular$t5 <- glabel("Retention time", container=env$le))
  env$le[10, 2:3] <- (env$specInfo$rt <- glabel("", container=env$le))
  env$le[11, 1] <- (env$regular$t6 <- glabel("Index", container=env$le))
  env$le[11, 2:3] <- (env$specInfo$ind <- glabel("", container=env$le, editable=TRUE))
  env$le[12, 1] <- (env$regular$t7 <- glabel("Precursor M/Z", container=env$le))
  env$le[12, 2:3] <- (env$specInfo$pmz <- glabel("", container=env$le))
  env$le[13, 1] <- (env$regular$t8 <- glabel("Precursor intensity", container=env$le))
  env$le[13, 2:3] <- (env$specInfo$int <- glabel("", container=env$le))
  
  env$le[14, 1:5] <- (env$separator$t3 <- glabel("", container=env$le))
  
  # Filters
  
  env$le[15, 1] <- (env$headings$t3 <- glabel("Filter", container=env$le))
  env$le[15, 2, anchor=c(-1,-1)] <- (env$deemp$t5 <- glabel("from", container=env$le))
  env$le[15, 3, anchor=c(-1,-1)] <- (env$deemp$t6 <- glabel("to", container=env$le))
  
  env$le[16, 1] <- (env$filterInfo$rtactive <- gcheckbox("Retention time", container=env$le))
  env$le[16, 2] <- (env$filterInfo$rtfrom <- gedit("", container=env$le, width=5))
  env$le[16, 3] <- (env$filterInfo$rtto <- gedit("", container=env$le, width=5))
  env$le[17, 1] <- (env$filterInfo$indexactive <- gcheckbox("Index", container=env$le))
  env$le[17, 2] <- (env$filterInfo$indexfrom <- gedit("", container=env$le, width=5))
  env$le[17, 3] <- (env$filterInfo$indexto <- gedit("", container=env$le, width=5))
  env$le[18, 1] <- (env$filterInfo$pmzactive <- gcheckbox("Precursor M/Z", container=env$le))
  env$le[18, 2] <- (env$filterInfo$pmzfrom <- gedit("", container=env$le, width=5))
  env$le[18, 3] <- (env$filterInfo$pmzto <- gedit("", container=env$le, width=5))
  env$le[19, 1] <- (env$filterInfo$spiactive <- gcheckbox("Precursor intensity", container=env$le))
  env$le[19, 2] <- (env$filterInfo$spifrom <- gedit("", container=env$le, width=5))
  env$le[19, 3] <- (env$filterInfo$spito <- gedit("", container=env$le, width=5))
  env$le[20, 1] <- (env$filterInfo$pcactive <- gcheckbox("Precursor charge", container=env$le))
  env$le[20, 2] <- (env$filterInfo$pcfrom <- gedit("", container=env$le, width=5))
  env$le[20, 3] <- (env$filterInfo$pcto <- gedit("", container=env$le, width=5))
  env$le[21, 1] <- (env$filterInfo$massactive <- gcheckbox("Precursor mass", container=env$le))
  env$le[21, 2] <- (env$filterInfo$massfrom <- gedit("", container=env$le, width=5))
  env$le[21, 3] <- (env$filterInfo$massto <- gedit("", container=env$le, width=5))
  
  tooltip(env$filterInfo$massactive) <- "sample tooltip"
  
  env$le[21, 1:5] <- (env$separator$t3 <- glabel("", container=env$le))
  
  env$le[22, 1, anchor=c(-1,0)] <- (env$regular$t12 <- glabel("Display MS levels", container=env$le))
  env$le[22, 2] <- (env$filterInfo$ms1 <- gcheckbox("MS1", checked=TRUE, container=env$le))
  env$le[22, 3] <- (env$filterInfo$ms2 <- gcheckbox("MS2", checked=TRUE, container=env$le))
  
  # Buttons
  
  env$groupLeftButtons <- ggroup(container=env$groupMiddleLeft, expand=TRUE)
  env$buttonLeft <- gbutton(text=gettext("Previous"), handler=env$updateSpectrum, 
                            action=-1, cont=env$groupLeftButtons)
  
  env$buttonRight <- gbutton(text=gettext("Next"), handler=env$updateSpectrum, 
                             action=1, cont=env$groupLeftButtons)
  
  if(any(c("cairo", "png")==device)) {
    
    env$plotTop <- ggraphics(container=env$groupPlots, width=500, height=250, ps=12, dpi=72)
    env$plotBottom <- ggraphics(container=env$groupPlots, width=500, height=250, ps=12, dpi=72)
    
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
  lapply(env$filterInfo, function(x) font(x) <- env$textReg)
  lapply(env$separator, function(x) font(x) <- list(size=4))
    
  if(device!="tkrplot") {    
    handlerClickGeneric = function(h,...) {
      cat("\nclicked on:", c(h$x, h$y))
    }
    handlerClickXIC = function(h,...) {
      cat("\nclicked on:", c(h$x, h$y))     
      
      prevCounter <- counter
      counter <<- which.min(abs(spRtime(currSequence) - h$x))
      
      # update graphs if index has changed
      if(prevCounter!=counter) updateSpectrum()      
      
    }
    env$clickHandlerIDs <- list()
    env$clickHandlerIDs[[1]] <- addHandlerClicked(env$plotBottom, handler=handlerClickXIC)
    env$clickHandlerIDs[[2]] <- addHandlerClicked(env$plotTop, handler=handlerClickGeneric)
  }

  env$filterSpectraHandlerIDs <- lapply(env$filterInfo, addHandlerChanged, handler=filterSpectra)
  
}
