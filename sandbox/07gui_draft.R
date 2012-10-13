library(gWidgets)
library(mzR)
source("../R/middle-layer.R")

# filename <- "d:/Dropbox/Documents/cambridge/r project/Data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"

options(guiToolkit = "RGtk2")
# options(guiToolkit = "tcltk")
 
.msGUIenv <- new.env(parent = emptyenv(), hash = TRUE)
assign("file.types", 
       list("Mass spectrometry data files" = list(patterns = c("*.netCDF", "*.mzXML", "*.mzData", "*.mzML")), 
            "R data files" = list(patterns = c("*.RData","*.rda")), 
            "All files" = list(patterns = c("*"))), 
       pos=.msGUIenv)
assign("digits", 2, pos=.msGUIenv)
assign("textHead", list(weight="bold", family="sans", 
                                size=12, color="grey04"), pos=.msGUIenv)
assign("textReg", list(weight="normal", family="sans", 
                                size=8, color="grey10"), pos=.msGUIenv)
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
#   )

drawMain <- function() {
  
  mget <- function(what) get(what, pos=.msGUIenv)
  mset <- function(what, value) assign(what, value, .msGUIenv)
    
  
  updateSpectrum <- function(h=list(action=0), ...) {
    counter <- mget("counter") + h$action
    mset("counter", counter)
    index <- mget("currSequence")[counter]
    svalue(specInfo$rt) <- spRtime()[index]
    svalue(specInfo$ind) <- spIndex()[index]
    svalue(specInfo$pmz) <- spPrecMz()[index]
    svalue(specInfo$int) <- mget("currSequence")[index]
    
    t <- proc.time()
    obj <- peaks(index)
    cat("\ndata: ", (proc.time()-t)[3])
    t <- proc.time()
    visible(plotTop) <- TRUE
#     par(mar=rep(0, 4))
    par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
    plot(obj, type = "l", frame.plot=FALSE, axes=TRUE)#, xlab="", ylab="")
    cat(" plot: ", (proc.time()-t)[3])
  }
  
  updateExperiment <- function() {
    rtRange <- round(expRtRange(), digits=mget("digits"))
    precMzRange <- round(expPrecMzRange(), digits=mget("digits"))
    svalue(expInfo$rtfrom) <- rtRange[1]
    svalue(expInfo$rtto) <- rtRange[2]
    svalue(expInfo$pmzfrom) <- precMzRange[1]
    svalue(expInfo$pmzto) <- precMzRange[2]
    svalue(expInfo$ms1) <- nMS1()
    svalue(expInfo$ms2) <- nMS2()
    
    index <- spIndex()
    mset("currSequence", index)
    mset("counter", index[1])
    enabled(buttonLeft) <- FALSE
    
    visible(plotBottom) <- TRUE
    par(mar=rep(0, 4))
    plot(xic(NULL), type = "l", frame.plot=FALSE, axes=TRUE) #, xlab="", ylab="")
  }
  
  # Window and structure
  msGUIWindow <- gwindow("msGUI", visible=FALSE)
  groupMain <- ggroup(container=msGUIWindow, horizontal=FALSE, expand=TRUE)
  groupUpper <- ggroup(container=groupMain, horizontal=TRUE)
  groupMiddle <- ggroup(container=groupMain, horizontal=TRUE)
  groupMiddleLeft <- ggroup(container=groupMiddle, horizontal=FALSE, spacing=5)
  groupMiddleRight <- gframe(container=groupMiddle, horizontal=FALSE, expand=TRUE)
  groupPlots <- gpanedgroup(container=groupMiddleRight, horizontal=FALSE, expand=TRUE)
  
  ## Top group
  buttonOpenFile <- gbutton("Open file", container=groupUpper, 
                            handler=function(h, ...) {
                              gfile(text = "Open file", type = "open", 
                                    initialfilename = "", 
                                    filter = get("file.types", .msGUIenv), 
                                    multi=FALSE, handler=function(h, ...) {
                                      makeMzRrampAccessors(h$file)
                                      updateExperiment() 
                                      updateSpectrum()  
                                    })
                            })
  
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
  
  le <- glayout(container=groupMiddleLeft, spacing=4)
  
  #Experiment info
  
  le[1, 1:5] <- (separator$t1 <- glabel("", container=le))
  
  le[2, 1:5] <- (headings$t1 <- glabel("Experiment information", container=le))
  le[3, 1] <- (regular$t1 <- glabel("Retention time:", container=le))
  le[3, 2] <- (deemp$t1 <- glabel("from", container=le))
  le[3, 3] <- (expInfo$rtfrom <- glabel(container=le))
  le[3, 4] <- (deemp$t2 <- glabel("to", container=le))
  le[3, 5] <- (expInfo$rtto <- glabel(container=le))
  le[4, 1] <- (regular$t2 <- glabel("Precursor M/Z:", container=le))
  le[4, 2] <- (deemp$t3 <- glabel("from", container=le))
  le[4, 3] <- (expInfo$pmzfrom <- glabel(container=le))
  le[4, 4] <- (deemp$t4 <- glabel("to", container=le))
  le[4, 5] <- (expInfo$pmzto <- glabel(container=le))
  le[5, 1] <- (regular$t3 <- glabel("Spectra count", container=le))
  le[5, 2] <- (regular$t4 <- glabel("MS1", container=le))
  le[5, 3] <- (expInfo$ms1 <- glabel(container=le))
  le[6, 2] <- (regular$t15 <- glabel("MS2", container=le))
  le[6, 3] <- (expInfo$ms2 <- glabel(container=le))
  
  le[7, 1:5] <- (separator$t2 <- glabel("", container=le))
  
  # Spectrum Info
  
  le[8, 1:5] <- (headings$t2 <- glabel("Spectrum information", container=le))
  le[9, 1] <- (regular$t5 <- glabel("Retention time:", container=le))
  le[9, 2] <- (specInfo$rt <- glabel("", container=le))
  le[10, 1] <- (regular$t6 <- glabel("Index", container=le))
  le[10, 2] <- (specInfo$ind <- glabel("", container=le))
  le[11, 1] <- (regular$t7 <- glabel("Precursor M/Z", container=le))
  le[11, 2] <- (specInfo$pmz <- glabel("", container=le))
  le[12, 1] <- (regular$t8 <- glabel("Precursor intensity", container=le))
  le[12, 2] <- (specInfo$int <- glabel("", container=le))
  
  le[13, 1:5] <- (separator$t3 <- glabel("", container=le))
  
#   labelFilter <- (headings$Filter <- glabel("Filter", container=groupMiddleLeft))
#   layoutFilter <- glayout(container=groupMiddleLeft, spacing=4)
  le[14, 1] <- (headings$t3 <- glabel("Filter", container=le))
  le[15, 2] <- (deemp$t5 <- glabel("from", container=le))
  le[15, 3] <- (deemp$t6 <- glabel("to", container=le))
  le[15, 5] <- (deemp$t7 <- glabel("active", container=le))
  
  le[16, 1] <- (regular$t9 <- glabel("Retention time:", container=le))
  le[16, 2] <- (filterInfo$rtfrom <- gedit(3, container=le, width=5))
  le[16, 3] <- (filterInfo$rtto <- gedit(5, container=le, width=5))
  le[16, 5] <- (filterInfo$rtactive <- gcheckbox(text="", checked=FALSE, 
                                                          container=le))
  le[17, 1] <- (regular$t10 <- glabel("Precursor M/Z", container=le))
  le[17, 2] <- (filterInfo$pmzfrom <- gedit(3, container=le, width=5))
  le[17, 3] <- (filterInfo$pmzto <- gedit(5, container=le, width=5))
  le[17, 5] <- (filterInfo$pmzactive <- gcheckbox(text="", checked=FALSE, 
                                                           container=le))
  le[18, 1] <- (regular$t11 <- glabel("Precursor intensity", container=le))
  le[18, 2] <- (filterInfo$spifrom <- gedit(3, container=le, width=5))
  le[18, 3] <- (filterInfo$spito <- gedit(5, container=le, width=5))
  le[18, 5] <- (filterInfo$spiactive <- gcheckbox(text="", checked=FALSE, 
                                                           container=le))
  le[19, 1] <- (regular$t12 <- glabel("MS levels", container=le))
  le[19, 2:4] <- (filterInfo$ms1 <- gcheckboxgroup(c("MS1", "MS2"), 
                                                            checked=TRUE, 
                                                            horizontal=TRUE, 
                                                            container=le))
  
  groupLeftButtons <- ggroup(container=groupMiddleLeft, horizontal=TRUE, expand=TRUE)
  buttonLeft <- gbutton(text=gettext("Previous"), 
                         handler=updateSpectrum, 
                         action=-1,
                         cont=groupLeftButtons)
  
  buttonRight <- gbutton(text=gettext("Next"), 
                             handler=updateSpectrum, 
                             action=1,
                             cont=groupLeftButtons)
  
  plotTop <- ggraphics(container=groupPlots, width=400, height=250, dpi=72)
  plotBottom <- ggraphics(container=groupPlots, width=400, height=250, dpi=72)

  lapply(headings, function(x) font(x) <- mget("textHead"))
  lapply(regular, function(x) font(x) <- mget("textReg"))
  lapply(deemp, function(x) font(x) <- mget("textDeemp"))
  lapply(expInfo, function(x) font(x) <- mget("textReg"))
  lapply(specInfo, function(x) font(x) <- mget("textReg"))
  lapply(filterInfo, function(x) font(x) <- mget("textReg"))
  lapply(separator, function(x) font(x) <- list(size=4))
  
  visible(msGUIWindow) <- TRUE
    
}

drawMain()