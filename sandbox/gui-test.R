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
source("R/msGUI.R")

filename <- "/media/Data/GitHub/msGUI/sandbox/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"

# generic
msGUI()

# MSnExp object 
msGUI(MSnbase::itraqdata)

# file
msGUI(filename)

# using options
msGUI(MSnbase::itraqdata, device="cairo", verbose=TRUE)