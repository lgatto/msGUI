library(gWidgets)
library(mzR)
library(png)
source("../R/middle-layer.R")
source("../R/initialisation.R")
source("../R/gui-functions.R")
source("../R/filtering-functions.R")
source("../R/plot-functions.R")

options(guiToolkit = "RGtk2")
# options(guiToolkit = "tcltk")


# Simply
wrapper()

# # or
# filename <- "d:/Dropbox/Documents/cambridge/r project/Data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"
# wrapper(filename)