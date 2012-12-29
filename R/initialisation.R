defaultSettings <- function() {
  
  fileTypes <- list("Mass spectrometry data files" = list(patterns = c("*.netCDF", "*.mzXML", "*.mzData", "*.mzML")), 
#                     "R data files" = list(patterns = c("*.RData","*.rda")), 
                    "All files" = list(patterns = c("*")))
  digits <- 2
  
  labelNumber <- 5
  labelThreshold <- 2
  
  MS1PlotType <- "l"
  MS2PlotType <- "h"
  
  chromaHeight <- 250
  spectrumHeight <- 250
  width <- 500
  
  chromaMode <- FALSE # i.e. total ion count is the default
  
  # Fonts
  fontHead <- list(weight="bold", family="sans", size=12, color="grey25")
  fontReg <- list(weight="light", family="sans", size=8, color="grey05")
  fontRegGrey <- list(weight="light", family="sans", size=8, color="grey50")
  fontDeemp <- list(weight="light", family="sans", size=8, color="grey50")
  
  return(environment())
}

initialiseGUI <- function() {  
  
  visible(msGUIWindow) <- TRUE  
  visible(plotTop) <- TRUE
  plotMsg()
  visible(plotBottom) <- TRUE
  plotMsg()
  enabled(buttonLeft) <- FALSE
  enabled(buttonRight) <- FALSE
  filterSwitch(FALSE)
  clickSwitch(FALSE)
}