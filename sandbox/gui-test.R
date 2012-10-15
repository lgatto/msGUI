filename <- "/home/lgatto/Data/Velos/2010_07_TMTspikes/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"
# filename <- "d:/Dropbox/Documents/cambridge/r project/Data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"

source("../R/middle-layer-simple.R")
source("07gui_draft.R")
makeMzRrampAccessors(filename)
drawMain(TRUE)