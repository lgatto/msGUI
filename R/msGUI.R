setGeneric("msGUI", function(object,...) standardGeneric("msGUI"))

setMethod("msGUI", "character",
          function(object, ...) {
            ## use mzR::openMSfile
            TRUE
          })

setMethod("msGUI", "MSnExp",
          function(object, ...) {
            ## use MSnbase infrastructure
            TRUE
          })
