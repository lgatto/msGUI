initialiseEnvironment <- function(env) {
  assign("file.types", 
         list("Mass spectrometry data files" = list(patterns = c("*.netCDF", "*.mzXML", "*.mzData", "*.mzML")), 
              "R data files" = list(patterns = c("*.RData","*.rda")), 
              "All files" = list(patterns = c("*"))), 
         pos=env)
  assign("digits", 2, pos=env)
  assign("textHead", list(weight="bold", family="sans", 
                          size=12, color="grey40"), pos=env)
  assign("textReg", list(weight="light", family="sans", 
                         size=8, color="grey10"), pos=env)
  assign("textRegGrey", list(weight="light", family="sans", 
                             size=8, color="grey50"), pos=env)
  assign("textDeemp", list(weight="light", family="sans", 
                           size=8, color="grey50"), pos=env)
}

initialiseGUI <- function() {  
  visible(msGUIWindow) <- TRUE  
  enabled(buttonLeft) <- FALSE
  enabled(buttonRight) <- FALSE
  filterSwitch(FALSE)
}