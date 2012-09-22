setGeneric("msGUI", function(object,...) standardGeneric("msGUI"))

setMethod("msGUI", "missing",
          function(object, ...) {
            ## opens GUI
            TRUE
          })


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
