source("../R/middle-layer.R")

filename <- "/home/lgatto/Data/Velos/2010_07_TMTspikes/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"
populateMsGuiEnv(filename)

ls(all=TRUE)
ls(.msGUIenv)
names(mzRrampAccessors)


par(mfrow = c(3,1))
plot(mzRrampAccessors$xic(1), type = "l")
plot(mzRrampAccessors$xic(2), type = "l")
plot(mzRrampAccessors$xic(NULL), type = "l") 

plot(mzRrampAccessors$peaks(100), type = "l")
