library(gWidgets)
library(MSnbase)
library(mzR)
library(png)
library(gWidgetsRGtk2)
library(cairoDevice)

source("R/middle-layer.R")
source("R/initialisation.R")
source("R/gui-functions.R")
source("R/filtering-functions.R")
source("R/plot-functions.R")
source("R/cache-functions.R")
source("sandbox/middle-layer-msnExp.R")

options(guiToolkit = "RGtk2")

filename <- "e:/GitHub/msGUI/sandbox/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"
wrapper(filename, device="png")
# wrapper(filename, device="cairo")
