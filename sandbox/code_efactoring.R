source.dir <- "R"
output.dir <- "new"
make.backup <- TRUE

changes <- list(c("drawMain", "getWindowMain"), 
                c("validityCheck", "validate"))

###############################################################################
               
file.names <- list.files(source.dir)

file.names <- file.names[!grepl("\\.bak$", file.names)]
if (make.backup) lapply(file.names, 
                        function(file) file.copy(file.path(source.dir, file), 
                                                 paste0(file.path(source.dir, file), ".bak")))

file.contents <- sapply(file.names, function(x) readLines(file.path(source.dir, x), warn=FALSE))

for (word in changes) file.contents <- lapply(file.contents, 
                                              function(file) {
                                                file <- gsub(paste0("(\\W)", word[1], "(\\W)"), 
                                                     paste0("\\1", word[2], "\\2"), 
                                                     file)
                                                gsub(paste0("^", word[1], "(\\W)"), 
                                                     paste0(word[2], "\\1"), 
                                                     file)
                                                })

if(!file.exists(output.dir)) dir.create(output.dir)
mapply(writeLines, file.contents, file.path(output.dir, file.names))

library(formatR)

tidy.dir(output.dir)

# Test

library(gWidgets)
library(MSnbase)
library(mzR)
library(png)
library(gWidgetsRGtk2)
library(cairoDevice)

sapply(file.path(output.dir, file.names), source)

options(guiToolkit = "RGtk2")

filename <- "e:/GitHub/msGUI/sandbox/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"
wrapper(filename, device="png")

