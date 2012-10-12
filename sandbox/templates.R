## middle layer definition

templates <- list("mzRramp" = list(),
                  "MSnExp" = list()
                  ## possibly more in the future
                  )

## a list of accessorts we want to define
## this would basically implement what you have compiled in the spread sheet

## Below:
## fh is the raw file handle,
## info, hd, ... would be stored/cached in an environment to avoid fetching them repeatedly
## These need to global variables, accessible in a enviroment created when the GUI is opened (say .msGUIenv)

templates[["mzRramp"]] <-
  list("rtime" = function() .msGUIenv$hd$retentionTime,
       "rtRange" = function(info) c(.msGUIenv$info$dStartTime, .msGUIenv$info$dEndTime)
       "peak" = function(fh, i) .msGUIenv$peaks(fh, i)
       ## ...
       )

templates[["MSnbase"]] <-
  list("rtime" = function() stop("Not yet implemented"),
       "rtRange" = function() stop("Not yet implemented")
       "peak" = function() stop("Not yet implemented")
       ## ...
       )

## the specific template would be set when the data is opened.
accessors <- templates[["mzRramp"]] ## for example if we work with a file

##them, to get the retention time, you use
accessors$rtime()


## The implementation above is a bit rough, but should work.
