setGeneric("msGUI", function(object,...) standardGeneric("msGUI"))

setMethod("msGUI", "missing",
          function(object, ...) {
            wrapper(...)
          })

setMethod("msGUI", "character",
          function(object, ...) {
            wrapper(filename=object, ...)
          })

## setMethod("msGUI", "MSnExp",
##           function(object, ...) {
##             wrapper(object=object, ...)
##           })
