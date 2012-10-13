require("mzR")  
source("../R/middle-layer.R")

filename <- "/home/lgatto/Data/Velos/2010_07_TMTspikes/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"

makeMzRrampAccessors(filename)

ls()

expRtRange()
nMS1()
nMS2()
expPrecMzRange()
spRtime(7297)
spRtime(7304)
spIndex(7297)
spIndex(7304)
spPrecMz(7304)
spPrecInt(7304)
spMsLevel(7297)
spMsLevel(7304)
spPrecCharge(7304)
spPrecScanNum(7304)
spBasePeakMz(7297)
spBasePeakMz(7304)
spBasePeakInt(7297)
spBasePeakInt(7304)
spPeaksCount(7297)
spPeaksCount(7304)

par(mfrow = c(4,1))
plot(xic(1), type = "l")
plot(xic(2), type = "l")
plot(xic(NULL), type = "l") 
plot(peaks(100), type = "l")


##
makeMSnExpAccessors()
try(expRtRange())
try(peaks(100))
