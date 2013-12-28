wrapper <- function(filename=NULL, object=NULL, device="png", verbose=FALSE) {
  env <- environment()
  env$settings <- defaultSettings()
  options(guiToolkit = "RGtk2")
  env$zoomWindowClosed <- TRUE
  env$XICWindowClosed <- TRUE
  env$optionsWindowClosed <- TRUE
  env$spectrumZoom <- NULL
  env$spectrumInt <- NULL
  env$XICZoom <- NULL
  env$XICInt <- NULL
  env$anyMS1spectra <- TRUE
  env$XICvalues <- TRUE
  env$experimentLoaded <- FALSE
  env$closingMsGUI <- FALSE
  env$forceRedraw <- FALSE
  env$clickMode <- TRUE
  drawMain(env)
  visible(env$msGUIWindow) <- TRUE  
  visible(env$plotTop) <- TRUE
  plotMsg()
  visible(env$plotBottom) <- TRUE
  plotMsg()
  enabled(env$buttonLeft) <- FALSE
  enabled(env$buttonRight) <- FALSE
  filterSwitch(FALSE, env)
  clickSwitch(FALSE, env)
  
  if(!is.null(filename)) {
    if (verbose) cat("Loading file...   ")
    makeMzRrampAccessors(filename, env)
    if (verbose) cat("done\n")
    updateExperiment(env)
    env$experimentLoaded <- TRUE
  }   
  
  if(!is.null(object)) {
    if (verbose) cat("Loading object...   ")
    makeMSnExpAccessors(object, env)
    if (verbose) cat("done\n")
    updateExperiment(env)
    env$experimentLoaded <- TRUE
  }
}

formatRt2 <- function(x, format, digits) {
  switch(format, 
         "minutes:seconds" = formatRt(x), 
         "minutes" = round(x/60, digits=digits), 
         "seconds" = round(x, digits=digits))
}

deformatRt <- function(x, format) {
  switch(format, 
         "minutes:seconds" = ifelse(grepl("^[0-9]*(:[0-9]{0,2})?$", x), {
            x <- unlist(strsplit(x, ":"))
            as.numeric(x[1])*60 + ifelse(length(x) > 1, as.numeric(x[2]), 0)}, NA), 
         "minutes" = ifelse(grepl("^[0-9]+([.][0-9]*)?$", x), as.numeric(x)*60, NA), 
         "seconds" = ifelse(grepl("^[0-9]+([.][0-9]*)?$", x), as.numeric(x), NA))
}

ident <- function(x, ...) x

updateExperimentInfo <- function(env) {
  precMzRange <- round(env$expPrecMzRange(), digits=env$settings$digits)
  svalue(env$expInfo$rtfrom) <- formatRt2(env$expRtRange()[1], env$settings$RtFormat, env$settings$digits)
  svalue(env$expInfo$rtto) <- formatRt2(env$expRtRange()[2], env$settings$RtFormat, env$settings$digits)
  svalue(env$expInfo$pmzfrom) <- precMzRange[1]
  svalue(env$expInfo$pmzto) <- precMzRange[2]
  svalue(env$expInfo$ms1) <- env$nMS1()
  svalue(env$expInfo$ms2) <- env$nMS2()
}

counterReset <- function(env) {    
  env$currSequence <- env$spIndex()
  env$counter <- 1
}

updateExperiment <- function(env) {  
  updateExperimentInfo(env)
  
  env$filterNames <- c("Retention time", "Index", "Prec MZ", "Prec intensity", 
                       "Prec charge", "Prec mass")
  
  # data cache. 
  env$filterData <- list(
    spRtime=env$spRtime(), 
    spIndex=env$spIndex(), 
    spPrecMz=env$spPrecMz(), 
    spPrecInt=env$spPrecInt(), 
    spPrecCharge=env$spPrecCharge(),
    spPrecMzTimesSpPrecCharge=env$spPrecMz()*env$spPrecCharge())
  
  env$filterTransform <- list(formatRt2, ident, ident,
                              ident, ident, ident)
  env$filterDetransform <- list(deformatRt, as.numeric, as.numeric,
                                as.numeric, as.numeric, as.numeric)
  
  env$nSpectra <- length(env$filterData[[1]])
  env$xLimits <- sapply(
    1:2, 
    function(i) 
      if(any(env$spMsLevel()==i)) 
        c(min(env$spLowMZ()[env$spMsLevel()==i]), 
          max(env$spHighMZ()[env$spMsLevel()==i])) 
      else rep(NA, 2)
    )
  env$anyMS1spectra <- any(env$spMsLevel()==1)
  if(env$anyMS1spectra) {
    dt <- cbind(env$spMsLevel(), env$spPrecMz())
    where <- which(dt[, 1]==1)
    frame <- cbind(where[-length(where)] + 1, where[-1] - 1)
    env$MS2indices <- frame[frame[, 1] < frame[, 2], ]
    
    env$xLimitsXIC <- range(env$xic()[, 1])
  } else {
    visible(env$plotBottom) <- TRUE
    plotMsg("Experiment contains no MS1 spectra.")
  }
  
  resetCache(env)
  counterReset(env)    
  filterReset(env)    
  filterSwitch(TRUE, env)
  clickSwitch(TRUE, env)
  updateSpectrum(h=list(action=list(0, env)))  
}

updateSpectrumInfo <- function(blank=FALSE, env) {
  if(blank) {
    lapply(list(env$specInfo$rt, env$specInfo$ind, env$specInfo$acno, env$specInfo$mslvl,
                env$specPrecInfo$pmz, env$specPrecInfo$int, env$specPrecInfo$charge), 
           function(i) svalue(i) <- "")
  } else {
    svalue(env$specInfo$rt) <- formatRt2(env$spRtime(env$index), env$settings$RtFormat, env$settings$digits)
    svalue(env$specInfo$ind) <- paste(env$counter, " of ", length(env$currSequence))
    svalue(env$specInfo$acno) <- env$spIndex(env$index)
    svalue(env$specInfo$mslvl) <- env$spMsLevel(env$index)  
    svalue(env$specPrecInfo$pmz) <- round(env$spPrecMz(env$index), digits=env$settings$digits)
    svalue(env$specPrecInfo$int) <- round(env$spPrecInt(env$index), digits=env$settings$digits)
    svalue(env$specPrecInfo$charge) <- round(env$spPrecCharge(env$index), digits=env$settings$digits) 
  }
}

updateSpectrum <- function(h, ...) {
  env <- h$action[[2]]
  # If called by buttons Previous or Next, then h$action will have value -1 or 1. 
  env$counter <- env$counter + h$action[[1]]
  env$index <- env$currSequence[env$counter]
  
  updateSpectrumInfo(FALSE, env)
  
  # kill precursor info for MS1
  dispPrec <- env$spMsLevel(env$index)==2  
  lapply(env$specPrecInfo, function(i) enabled(i) <- dispPrec) 
  lapply(env$regularPrec, function(i) enabled(i) <- dispPrec)
  
  # Turn off buttons, click handlers and filters while drawing
  buttonSwitch(FALSE, env)
  clickSwitch(FALSE, env)  
  
  plotSpectrum(zoom=env$spectrumZoom, int=env$spectrumInt, env=env)
  
  if(!env$zoomWindowClosed) {
    visible(env$plotZoom) <- TRUE      
    plotSpectrumZoom(env$spectrumZoom, env$spectrumInt, env)
  }
  
  if(env$anyMS1spectra) {
    plotXIC(env$XICZoom, env$XICInt, env=env) 
    if(!env$XICWindowClosed) {
      visible(env$plotXICw) <- TRUE      
      plotChromaZoom(env)
    }
  }
  
  clickSwitch(TRUE, env)
  if (env$forceRedraw) clickSwitch(TRUE, env)
  buttonSwitch(TRUE, env)
}

clickSwitch <- function(on=TRUE, env) {
  fun <- if (on) unblockHandler else blockHandler
  mapply(fun, list(env$plotTop, env$plotBottom), env$clickHandlerIDs)
  if (!env$zoomWindowClosed) fun(env$plotZoom, env$clickHandlerZoomIDs)
  if (!env$XICWindowClosed) fun(env$plotXICw, env$clickHandlerZoomXICIDs)
  if (!env$anyMS1spectra) blockHandler(env$plotBottom, env$clickHandlerIDs[[2]])
}

buttonSwitch <- function(on=TRUE, env) {
  # (Selectively) turn on Next/Previous buttons
  if(on) {    
    # If current spectra selection has a single spectrum, both buttons disabled
    if(length(env$currSequence)==1) {
      enabled(env$buttonLeft) <- FALSE
      enabled(env$buttonRight) <- FALSE      
    } else {
      # If we have multiple spectra and the first is current
      if(env$counter==1) {
        enabled(env$buttonLeft) <- FALSE
        enabled(env$buttonRight) <- TRUE
      } else {
        # ... the last one is current...
        if(env$counter==length(env$currSequence)) {
          enabled(env$buttonLeft) <- TRUE
          enabled(env$buttonRight) <- FALSE
          # ...or we're in the middle
        } else {
          enabled(env$buttonLeft) <- TRUE
          enabled(env$buttonRight) <- TRUE
        }
      }
    }
  } else {
    enabled(env$buttonLeft) <- FALSE
    enabled(env$buttonRight) <- FALSE
  }
}

openFileHandler <- function(h, ...) {
  env <- h$action
  filename <- gfile(text = "Open file", type = "open", 
                    filter = env$settings$fileTypes, 
                    handler=function(h, ...) NULL)
  if(!is.na(filename)) {
    if(env$experimentLoaded) {      
      visible(env$plotTop) <- TRUE
      plotMsg()
      visible(env$plotBottom) <- TRUE
      plotMsg()
    }
    if (env$verbose) cat("Loading file...   ")
    makeMzRrampAccessors(filename, env)
    if (env$verbose) cat("done\n")
    updateExperiment(env)    
    env$experimentLoaded <- TRUE
  }
}

drawMain <- function(env) {
  # Window and structure
  env$msGUIWindow <- gwindow("msGUI", visible=FALSE, action=env, handler=function(h, ...) {
    env <- h$action
    env$closingMsGUI <- TRUE
    if(!env$zoomWindowClosed) dispose(env$zoomWindow)
    if(!env$XICWindowClosed) dispose(env$XICWindow)
    # Clean-up cache
    if(class(env$cache)!="function") {
      files <- unlist(c(env$cache$spectra[!sapply(env$cache$spectra, is.null)], 
                        env$cache$xic[!sapply(env$cache$xic, is.null)]))
      if(length(files)>0) file.remove(files)      
    }
    # Clean-up graphics devices    
    sapply(dev.list()[grepl("^png", names(dev.list()))], dev.off)
  })
  env$groupMain <- ggroup(container=env$msGUIWindow, horizontal=FALSE) 
  env$groupUpper <- ggroup(container=env$groupMain)
  env$groupMiddle <- ggroup(container=env$groupMain) 
  env$groupMiddleLeft <- ggroup(container=env$groupMiddle, horizontal=FALSE)
  env$groupMiddleRight <- ggroup(container=env$groupMiddle, horizontal=FALSE)
  env$groupPlots <- ggroup(container=env$groupMiddleRight, horizontal=FALSE)
  
  ## Top group
  env$buttonOpenFile <- gbutton("Open file", container=env$groupUpper, 
                                handler=openFileHandler, action=env)
#   env$buttonOpenObject <- gbutton("Open R object", container=groupUpper, 
#                                   handler=drawVarBrowser, action=env)
  addSpring(env$groupUpper)
  env$buttonSettings <- gbutton("Settings", container=env$groupUpper, 
                                handler=drawOptions, action=env)
  env$buttonHelp <- gbutton("Help", container=env$groupUpper)
  
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
  env$tools <- list()
  
  # Left pane content  
  env$le <- glayout(container=env$groupMiddleLeft, spacing=env$settings$uiGridSpacing)
  
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
  env$le[i + 2, 2:3] <- (env$specInfo$ind <- glabel("", container=env$le))
  env$le[i + 3, 1] <- (env$regular$t16 <- glabel("Acquisition number", container=env$le))
  env$le[i + 3, 2:3] <- (env$specInfo$acno <- glabel("", container=env$le))
  env$le[i + 4, 1] <- (env$regular$t17 <- glabel("MS level", container=env$le))
  env$le[i + 4, 2:3] <- (env$specInfo$mslvl <- glabel("", container=env$le))
  env$le[i + 5, 1] <- (env$regularPrec$t1 <- glabel("Precursor M/Z", container=env$le))
  env$le[i + 5, 2:3] <- (env$specPrecInfo$pmz <- glabel("", container=env$le))
  env$le[i + 6, 1] <- (env$regularPrec$t2 <- glabel("Precursor intensity", container=env$le))
  env$le[i + 6, 2:3] <- (env$specPrecInfo$int <- glabel("", container=env$le))
  env$le[i + 7, 1] <- (env$regularPrec$t3 <- glabel("Precursor charge", container=env$le))
  env$le[i + 7, 2:3] <- (env$specPrecInfo$charge <- glabel("", container=env$le))
  
  env$le[i + 8, 1:5] <- (env$separator$t3 <- glabel("", container=env$le))
  
  # Filters  
  i <- 18
  
  env$le[i, 1] <- (env$headings$t3 <- glabel("Filter", container=env$le))
  env$le[i, 2, anchor=c(-1,-1)] <- (env$deemp$t5 <- glabel("from", container=env$le))
  env$le[i, 3, anchor=c(-1,-1)] <- (env$deemp$t6 <- glabel("to", container=env$le))
  
  env$le[i + 1, 1] <- (env$filterInfo$rt$active <- gcheckbox("Retention time", container=env$le))
  env$le[i + 1, 2] <- (env$filterInfo$rt$from <- gedit("", container=env$le, width=5))
  env$le[i + 1, 3] <- (env$filterInfo$rt$to <- gedit("", container=env$le, width=5))
  env$le[i + 2, 1] <- (env$filterInfo$index$active <- gcheckbox("Acquisition number", container=env$le))
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
  
  env$le[i + 8, 1, anchor=c(-1,0)] <- (env$regular$t12 <- glabel("Display MS levels", 
                                                                 container=env$le))
  env$le[i + 8, 2] <- (env$filterInfoMS$ms1 <- gcheckbox("MS1", checked=TRUE, 
                                                         container=env$le, 
                                                         use.togglebutton=TRUE))
  env$le[i + 8, 3] <- (env$filterInfoMS$ms2 <- gcheckbox("MS2", checked=TRUE, 
                                                         container=env$le, 
                                                         use.togglebutton=TRUE))
  
  env$le[i + 9, 1:5] <- (env$separator$t9 <- glabel("", container=env$le))
  
  env$le[i + 10, 1] <- (env$filterInfoXIC$XIC$active <- gcheckbox("Prec M/Z for XIC", container=env$le))
  env$le[i + 10, 2] <- (env$filterInfoXIC$XIC$x <- gedit("", container=env$le, width=5, coerce.with=as.numeric))
  env$le[i + 10, 3, anchor=c(-1, 0)] <- (env$filterInfoXIC$XIC$text <- glabel("+/-.5 Da", container=env$le))
  
  env$le[i + 11, 1:5] <- (env$separator$t19 <- glabel("", container=env$le))
  
  # Tools  
  env$le[i + 12, 1] <- (env$headings$t37 <- glabel("Tools", container=env$le))

  env$groupClickMode <- ggroup(container=env$groupMiddleLeft)
  env$tools$zoom <- gcheckbox("Zoom", checked=env$clickMode, 
                              container=env$groupClickMode, 
                              use.togglebutton=TRUE)
  env$tools$integrate <- gcheckbox("Integrate", checked=!env$clickMode, 
                                   container=env$groupClickMode, 
                                   use.togglebutton=TRUE)
  
  # Buttons  
  addSpring(env$groupMiddleLeft)
  
  env$groupLeftButtons <- ggroup(container=env$groupMiddleLeft)
  env$buttonLeft <- gbutton(text=gettext("Previous"), handler=updateSpectrum, 
                            action=list(-1, env), container=env$groupLeftButtons)
  
  env$buttonRight <- gbutton(text=gettext("Next"), handler=updateSpectrum, 
                             action=list(1, env), container=env$groupLeftButtons)
  
  env$plotTop <- ggraphics(container=env$groupPlots, width=env$settings$width,
                           height=env$settings$spectrumHeight, ps=12, dpi=75)
  env$plotBottom <- ggraphics(container=env$groupPlots, width=env$settings$width, 
                              height=env$settings$chromaHeight, ps=12, dpi=75)
  
  # Styling  
  lapply(env$headings, setFont, env$settings$fontHead)
  lapply(env$regular, setFont, env$settings$fontReg)
  lapply(env$regularPrec, setFont, env$settings$fontReg)
  lapply(env$deemp, setFont, env$settings$fontDeemp)
  lapply(env$expInfo, setFont, env$settings$fontReg)
  lapply(env$specInfo, setFont, env$settings$fontReg)
  lapply(env$separator, setFont, list(size=4))
  lapply(env$specPrecInfo, setFont, env$settings$fontReg)
  lapply(env$filterInfo, function(x) lapply(x, setFont, env$settings$fontReg)) 
  lapply(env$filterInfoXIC, function(x) lapply(x, setFont, env$settings$fontReg)) 
  lapply(env$filterInfo, function(x) setFontGtk(x$active, env$settings$fontReg)) 
  lapply(env$filterInfoXIC, function(x) setFontGtk(x$active, env$settings$fontReg))
  
  # Zoom handlers and GUI functions  
  drawZoom <- function(env) {    
    env$zoomWindowClosed <- FALSE
    env$zoomWindow <- gwindow("Spectrum fragment", visible=FALSE, parent=env$msGUIWindow, 
                              handler=function(h, ...) {
                                env <- h$action
                                env$zoomWindowClosed <- TRUE
                                env$spectrumZoom <- NULL
                                if(!env$closingMsGUI) plotSpectrum(env=env)
                              }, action=env)
    env$groupZoomMain <- ggroup(container=env$zoomWindow, horizontal=FALSE, expand=TRUE)
    env$plotZoom <- ggraphics(width=250, height=250, container=env$groupZoomMain, dpi=96, ps=12)
    visible(env$zoomWindow) <- TRUE
    env$clickHandlerZoomIDs <- addHandlerChanged(env$plotZoom, handler=handlerClickZoom, action=env)
  }   
  
  drawZoomXIC <- function(env) { 
    env$XICWindowClosed <- FALSE
    env$XICWindow <- gwindow("XIC", visible=FALSE, parent=env$msGUIWindow, 
                             handler=function(h, ...) { 
                               env <- h$action
                               env$XICWindowClosed <- TRUE
                               env$XICZoom <- NULL
                               if(!env$closingMsGUI) plotXIC(env=env)
                             }, action=env)
    env$groupXICMain <- ggroup(container=env$XICWindow, horizontal=FALSE, expand=TRUE)
    env$plotXICw <- ggraphics(width=250, height=250, container=env$groupXICMain, dpi=96, ps=12)
    visible(env$XICWindow) <- TRUE
    env$clickHandlerZoomXICIDs <- addHandlerChanged(env$plotXICw, handler=handlerClickZoomXIC, action=env)
  } 
  
  fixX <- function(x, lower, upper, width, device) {
    if (device == "png") {
      x <- ((x + .04) / 1.08 - 50/width) * width / (width - 50 - 25) * (upper - lower) + lower
      sapply(sapply(x, max, lower), min, upper) # so that coordinates don't exceed limits
    } else x
  }
  
  fixY <- function(x, lower, upper, height, device) {
    if (device == "png") {
      x <- ((x + .02) / 1.08 - 40/height) * height / (height - 40 - 22) * (upper - lower) + lower
      sapply(sapply(x, max, lower), min, upper)
    } else x
  }  
  
  handlerClickSpectrum = function(h,...) {
    env <- h$action
    xlim <- env$xLimits[, env$spMsLevel(env$index)];
    coords <- list(x=fixX(h$x, xlim[1], xlim[2], env$settings$width, env$device),
                   y=fixY(h$y, 0, 1.05, env$settings$spectrumHeight, env$device))
    if(env$clickMode) {
      env$spectrumZoom <- coords
    } else {
      env$spectrumInt <- coords
    }
    plotSpectrum(zoom=env$spectrumZoom, int=env$spectrumInt, env=env)
    if(!is.null(env$spectrumZoom)) {
      if(env$zoomWindowClosed) drawZoom(env)
      visible(env$plotZoom) <- TRUE
      plotSpectrumZoom(env$spectrumZoom, env$spectrumInt, env)
    }
  }  
  
  handlerClickXIC <- function(h,...) {
    env <- h$action
    clickSwitch(FALSE, env)
    xicRangeX <- range(env$xic(n=1, FALSE)[, 1])
    xCoord <- fixX(h$x, xicRangeX[1], xicRangeX[2], env$settings$width, env$device)[2]
    
    if (env$verbose) cat("coords: x", h$x, "y", h$y, "recalculated coords", xCoord, "\n")
    
    if (same(h[c("x", "y")])) {
      prevCounter <- env$counter    
      env$counter <- which.min(abs(env$spRtime(env$currSequence)-xCoord))      
      # update graphs if index has changed
      if(prevCounter!=env$counter) updateSpectrum(h=list(action=list(0, env)))       
    } else {
      coords <- list(x=fixX(h$x, xicRangeX[1], xicRangeX[2], env$settings$width, env$device), 
                     y=fixY(h$y, 0, 1.05, env$settings$chromaHeight, env$device))
      if (env$clickMode) {
        env$XICZoom <- coords
      } else {
        env$XICInt <- coords
      }
      if (env$XICWindowClosed & env$clickMode) {
        drawZoomXIC(env) 
        # & env$clickMode prevents the zoom window opening in integration mode. 
        visible(env$plotXICw) <- TRUE
        plotChromaZoom(env)
      } else if (!env$XICWindowClosed) {
        visible(env$plotXICw) <- TRUE
        plotChromaZoom(env)
      }
      plotXIC(env$XICZoom, env$XICInt, env=env)
    }        
    clickSwitch(TRUE, env)   
  } 
  
  handlerClickZoom <- function(h,...) {
    env <- h$action
    coords <- h[c("x", "y")]
    if(env$clickMode) {
      env$spectrumZoom <- coords
    } else {
      env$spectrumInt <- coords
    }
    plotSpectrum(zoom=env$spectrumZoom, int=env$spectrumInt, env=env)
    if(env$zoomWindowClosed) drawZoom(env)
    visible(env$plotZoom) <- TRUE          
    plotSpectrumZoom(env$spectrumZoom, env$spectrumInt, env)
  } 
  
  handlerClickZoomXIC <- function(h, ...) {
    env <- h$action
    coords <- h[c("x", "y")]
    if(env$verbose) cat("coords: x", h$x, "y", h$y, "\n")
    if(env$clickMode) {
      env$XICZoom <- coords
    } else {
      env$XICInt <- coords
    }
#     env$XICZoom <- h    
    plotXIC(zoom=env$XICZoom, int=env$XICInt, env=env)    
    if(env$XICWindowClosed) drawZoomXIC(env)
    visible(env$plotXICw) <- TRUE          
    plotChromaZoom(env)
  }
  
  handlerClickMode <- function(h, ...) {
    mapply(blockHandler, env$tools, env$toolsHandlerIDs)
    env <- h$action
    env$clickMode <- !env$clickMode
    svalue(env$tools$zoom) <- env$clickMode
    svalue(env$tools$integrate) <- !env$clickMode
    mapply(unblockHandler, env$tools, env$toolsHandlerIDs)
  }
  
  env$clickHandlerIDs <- list()
  env$clickHandlerIDs[[1]] <- addHandlerChanged(env$plotTop, handler=handlerClickSpectrum, action=env)
  env$clickHandlerIDs[[2]] <- addHandlerChanged(env$plotBottom, handler=handlerClickXIC, action=env)
  
  env$filterSpectraHandlerIDs <- lapply(env$filterInfo, 
                                        function(i) lapply(i, addHandlerChanged, 
                                                           handler=filterSpectra, action=env))
  
  env$filterSpectraMSHandlerIDs <- lapply(env$filterInfoMS, addHandlerChanged, 
                                          handler=filterSpectra, action=env)
  env$filterXICHandlerIDs <- lapply(env$filterInfoXIC, 
                                    function(i) lapply(i, addHandlerChanged, 
                                                       handler=filterSpectra, action=env))
  env$toolsHandlerIDs <- lapply(env$tools, addHandlerChanged, 
                                handler=handlerClickMode, action=env)
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

openObject <- function(object, env) {
  if(env$verbose) cat("Opening object...    ", object)
  if(env$experimentLoaded) {      
    visible(env$plotTop) <- TRUE
    plotMsg()
    visible(env$plotBottom) <- TRUE
    plotMsg()
  }
  makeMSnExpAccessors(get(object), env)
  if(env$verbose) cat("done!\n")
  updateExperiment(env)
}

drawVarBrowser <- function(h, ...) {
  env <- h$action
  windowVB <- gwindow(title="Browse R objects", visible=FALSE, 
                      width=400, height=300, parent=env$msGUIWindow)
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
  addSpring(env$panelVBbottom)
  buttonOpen <- gbutton(text="Open", container=panelVBbottom, 
                        handler=openObjectHandler)
  visible(windowVB) <- TRUE
}

optsHandlerDefaults <- function(h, ...) {
  env <- h$action
  defaults <- defaultSettings()
  svalue(env$opts$spectrumHeight) <- defaults$spectrumHeight
  svalue(env$opts$chromaHeight) <- defaults$chromaHeight
  svalue(env$opts$width) <- defaults$width
  svalue(env$opts$labelNumber) <- defaults$labelNumber
  svalue(env$opts$MS1PlotType) <- ifelse(defaults$MS1PlotType=="h", "", " ")
  svalue(env$opts$MS2PlotType) <- ifelse(defaults$MS2PlotType=="h", "", " ")
  svalue(env$opts$chromaMode) <- ifelse(defaults$chromaMode, "Base peak intensity", 
                                    "Total ion count")
  svalue(env$opts$Da) <- defaults$Da
  svalue(env$opts$labelDist) <- defaults$labelDist
  svalue(env$opts$RtFormat) <- defaults$RtFormat
}

drawOptions <- function (h, ...) {
  env <- h$action
  
  if(!env$optionsWindowClosed) return(NULL)
  
  env$optionsWindowClosed <- FALSE
  
  env$optsWindow <- gwindow("Settings", visible=FALSE, height=50, width=50, 
                            parent=env$msGUIWindow, handler=function(h, ...) {
                              h$action$optionsWindowClosed <- TRUE
                            }, action=env)
  env$optsGroup <- ggroup(container=env$optsWindow, horizontal=FALSE)
  env$l <- glayout(container=env$optsGroup, homogeneous=TRUE, spacing=30)
  env$l[1, 1] <- (env$l1 <- glayout(container=env$l, spacing=env$settings$uiGridSpacing))
  env$l[1, 2] <- (env$l2 <- glayout(container=env$l, spacing=env$settings$uiGridSpacing))
  
  env$l1[1, 1:3] <- (env$opts$headings$t1 <- glabel("Graph settings", container=env$l1))
  env$l1[2, 1, anchor=c(-1, 0)] <- (env$opts$text$t19 <- glabel("Size", container=env$l1)) 
  env$l1[2, 2, anchor=c(0, 0)] <- (env$opts$text$t1 <- glabel("height", container=env$l1))  
  env$l1[2, 3, anchor=c(0, 0)] <- (env$opts$text$t2 <- glabel("width", container=env$l1))
  env$l1[3, 2] <- (env$opts$spectrumHeight <- gspinbutton(from=250, to=500, by=10, 
                                                          value=env$settings$spectrumHeight, 
                                                          digits=0, container=env$l1))
  env$l1[3, 3] <- (env$opts$width <- gspinbutton(from=500, to=800, by=10, 
                                                 value=env$settings$width, digits=0, 
                                                 container=env$l1))
  env$l1[3, 1, anchor=c(-1, 0)] <- (env$opts$text$t1 <- glabel("Spectrum", 
                                                               container=env$l1))
  env$l1[4, 1, anchor=c(-1, 0)] <- (env$opts$text$t2 <- glabel("Chromatogram", 
                                                               container=env$l1))
  env$l1[4, 2] <- (env$opts$chromaHeight <- gspinbutton(from=250, to=500, by=10, 
                                                        value=env$settings$chromaHeight, 
                                                        digits=0, container=env$l1))  
  i <- 5
  
  env$l1[i, 1:3] <- (env$opts$separator$t1 <- glabel("", container=env$l1))
  env$l1[i + 1, 1:2, anchor=c(-1, 0)] <- (env$opts$text$t3 <- glabel("# of peak labels", 
                                                                     container=env$l1))
  env$l1[i + 1, 3] <- (env$opts$labelNumber <- gedit(env$settings$labelNumber, 
                                                     coerce.with=as.numeric, 
                                                     width=2, container=env$l1))
  env$l1[i + 2, 1:2, anchor=c(-1, 0)] <- (env$opts$text$t16 <- glabel("Distance labels by", 
                                                                      container=env$l1))
  env$l1[i + 2, 3] <- (env$opts$labelDist <- gedit(env$settings$labelDist, 
                                            coerce.with=as.numeric, 
                                            width=2, container=env$l1))
  env$l1[i + 3, 1:3, anchor=c(-1, 0)] <- (env$opts$text$t17 <- glabel("Retention time in", 
                                                                    container=env$l1))
  env$l1[i + 4, 2:3] <- (env$opts$RtFormat <- gdroplist(env$settings$RtFormats, 
                                                        match(env$settings$RtFormat, env$settings$RtFormats),  
                                                        container=env$l1))  
  env$l1[i + 5, 1:3] <- (env$opts$separator$t2 <- glabel("", container=env$l1))
  
  env$l2[1, 1:3] <- (env$opts$headings$t2 <- glabel("MS mode", container=env$l2))
  env$l2[2, 2, anchor=c(1, 0)] <- (env$opts$text$t4 <- glabel("MS1     ", container=env$l2))
  env$l2[2, 3, anchor=c(1, 0)] <- (env$opts$text$t5 <- glabel("MS2     ", container=env$l2))
  env$l2[3:4, 2] <- (env$opts$MS1PlotType <- gradio(c("", " "), 
                                                    ifelse(env$settings$MS1PlotType=="h", 1, 2), 
                                                    container=env$l2))
  env$l2[3:4, 3] <- (env$opts$MS2PlotType <- gradio(c("", " "), 
                                                    ifelse(env$settings$MS2PlotType=="h", 1, 2), 
                                                    container=env$l2))
  env$l2[3, 1, anchor=c(-1, 0)] <- (env$opts$text$t6 <- glabel("Centroided             ", container=env$l2))
  env$l2[4, 1, anchor=c(-1, 0)] <- (env$opts$text$t7 <- glabel("Profile", container=env$l2))
  
  i <- 5
  
  env$l2[i, 1:3] <- (env$opts$separator$t3 <- glabel("", container=env$l2))
  env$l2[i + 1, 1:3] <- (env$opts$headings$t4 <- glabel("Chromatogram mode", 
                                                        container=env$l2))
  env$l2[i + 2, 1:3] <- (env$opts$chromaMode <- gradio(c("Total ion count", 
                                                         "Base peak intensity"), 
                                                       ifelse(env$settings$chromaMode, 2, 1), 
                                                       container=env$l2))
  env$l2[i + 3, 1:3] <- (env$opts$separator$t5 <- glabel("", container=env$l2))
  env$l2[i + 4, 1:3] <- (env$opts$headings$t7 <- glabel("Filtering", container=env$l2))
  env$l2[i + 5, 1:2, anchor=c(-1, 0)] <- (env$opts$text$t16 <- glabel("Da size", 
                                                                      container=env$l1))
  env$l2[i + 5, 3] <- (env$opts$Da <- gedit(env$settings$Da, 
                                            coerce.with=as.numeric, 
                                            width=2, container=env$l1))
  env$l2[i + 6, 1:3] <- (env$opts$separator$t4 <- glabel("", container=env$l2))
  
  env$optsButtons <- ggroup(container=env$optsGroup, horizontal=TRUE)
  env$opts$defaults <- gbutton("Restore defaults", container=env$optsButtons, 
                               width=50, handler=optsHandlerDefaults, action=env)
  
  addSpring(env$optsButtons)
  
  env$opts$ok <- gbutton("OK", container=env$optsButtons, 
                         width=env$settings$minButtonWidth, 
                         handler=function(h, ...) {
    if(env$experimentLoaded & 
         any(c(env$settings$spectrumHeight!=svalue(env$opts$spectrumHeight), 
             env$settings$chromaHeight!=svalue(env$opts$chromaHeight), 
             env$settings$width!=svalue(env$opts$width), 
             env$settings$labelNumber!=svalue(env$opts$labelNumber), 
             env$settings$MS1PlotType!=ifelse(svalue(env$opts$MS1PlotType)=="", "h", "l"), 
             env$settings$MS2PlotType!=ifelse(svalue(env$opts$MS2PlotType)=="", "h", "l"), 
             env$settings$chromaMode!=(svalue(env$opts$chromaMode)=="Base peak intensity"), 
             env$settings$Da!=svalue(env$opts$Da), 
             env$settings$RtFormat!=svalue(env$opts$RtFormat), 
             env$settings$labelDist!=svalue(env$opts$labelDist)
    ))) {
      if(env$verbose) cat("Applying changes...\n")
      env$settings$spectrumHeight <- svalue(env$opts$spectrumHeight)
      env$settings$chromaHeight <- svalue(env$opts$chromaHeight)
      env$settings$width <- svalue(env$opts$width)
      env$settings$labelNumber <- svalue(env$opts$labelNumber)
      env$settings$MS1PlotType <- ifelse(svalue(env$opts$MS1PlotType)=="", "h", "l")
      env$settings$MS2PlotType <- ifelse(svalue(env$opts$MS2PlotType)=="", "h", "l")
      env$settings$chromaMode <- (svalue(env$opts$chromaMode)=="Base peak intensity")
      env$settings$Da <- svalue(env$opts$Da)
      prevRtFormat <- env$settings$RtFormat
      env$settings$RtFormat <- svalue(env$opts$RtFormat)
      env$settings$labelDist <- svalue(env$opts$labelDist)
      
      updateExperimentInfo(env)
      blockFilters(TRUE, env)
      svalue(env$filterInfo$rt$from) <-formatRt2(deformatRt(svalue(env$filterInfo$rt$from), 
                                                          prevRtFormat), 
                                                 env$settings$RtFormat, env$settings$digits)
      svalue(env$filterInfo$rt$to) <-formatRt2(deformatRt(svalue(env$filterInfo$rt$to), 
                                                          prevRtFormat), 
                                               env$settings$RtFormat, env$settings$digits)
      filterSpectra(h=list(action=env))
      blockFilters(FALSE, env)
      
      size(env$plotTop) <- c(env$settings$width, env$settings$spectrumHeight)
      size(env$plotBottom) <- c(env$settings$width, env$settings$chromaHeight)
      if(!env$zoomWindowClosed) {
        visible(env$plotZoom) <- TRUE      
        plotSpectrumZoom(env$spectrumZoom, env$spectrumInt, env)       
      }
      if(!env$XICWindowClosed) {
        visible(env$plotXICw) <- TRUE      
        plotChromaZoom(env)        
      }        
      resetCache(env)
      updateSpectrum(h=list(action=list(0, env)))
    }
    dispose(env$optsWindow)
  })
  env$opts$cancel <- gbutton("Cancel", container=env$optsButtons, width=50, 
                             handler=function(h, ...) dispose(h$action$optsWindow), action=env)
  addSpace(env$optsButtons, 18)
  
  lapply(env$opts$headings, setFont, env$settings$fontHead2)
  
  visible(env$optsWindow) <- TRUE
}

setFont <- function(x, style) font(x) <- style

same <- function(h) abs(h$x[1]-h$x[2]) < 1/200 & abs(h$y[1]-h$y[2]) < 1/200

# Workaround for widget styling. Hat tip to J Verzani: http://stackoverflow.com/questions/14940349/how-to-set-font-for-a-gcheckbox-object/
setFontGtk <- function(object, spec) {
  require(RGtk2)
  widget <- getToolkitWidget(object)$getChildren()[[1]]
  font_descr <- pangoFontDescriptionNew()
  if(!is.null(spec$weight))
    font_descr$setWeight(PangoWeight[spec$weight])
  if(!is.null(spec$style))
    font_descr$setStyle(PangoStyle[spec$style])
  if(!is.null(spec$size))
    font_descr$setSize(spec$size * PANGO_SCALE)
  if(!is.null(spec$family))
    font_descr$setFamily(spec$family)
  widget$modifyFont(font_descr)
  
  if(!is.null(spec$color))
    widget$modifyFg(GtkStateType[1], spec$color)
}
