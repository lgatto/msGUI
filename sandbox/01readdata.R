file <- "d:/Dropbox/Documents/cambridge/r project/Data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"

require(MSnbase)

m <- openMSfile(file)
object.size(m)
runInfo(m)

exp_info <- runInfo(m)

# m1 <- readMSData(file, msLevel=1)

# object.size(m1)
# length(m1)
# m1h <- header(m1)

m2 <- readMSData(file, msLevel=2)

object.size(m2)
length(m2)

# m2

m2h <- header(m2)
# View(m2h)


#### Later. 
# import_file <- function(file) {
#   extension <- rev(unlist(strsplit(file, "\\.")))[[1]]
#   extension <- tolower(extension)
#   # If RData, need to save current state of objects
#   if(extension%in%c("rda", "rdata")) {
#     
#   }
# }