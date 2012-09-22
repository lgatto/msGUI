\name{msGUI-methods}
\docType{methods}
\alias{msGUI-methods}
\alias{msGUI,character-method}
\alias{msGUI,MSnExp-method}
\alias{msGUI,missing-method}

\title{ Starts the msGUI software }

\description{
  The \code{msGUI} interface allows one to browser MS spectra and
  display extracted ion chromatograms from raw data files or
  \code{"\linkS4class{MSnExp}"} instances. 
}


\section{Methods}{
  \describe{
    \item{\code{signature(object = "missing")}}{
      Opens the graphical interface without displaying any data.
    }
    \item{\code{signature(object = "character")}}{
      \code{object} is the name of a \code{mzData}, \code{mzXML} or
      \code{mzML} file to be open and dispayed.
    }
    \item{\code{signature(object = "MSnExp")}}{
      \code{object} is the \code{"\linkS4class{MSnExp}"} instance to be
      displayed. 
    }
  }
}


\seealso{
  \code{mzR} package for access to raw data and the
  \code{"\linkS4class{MSnExp}"} class from the \code{MSnbase} package
  for a convenient container for such data.
}

\references{
  Cite \code{mzData}, \code{mzXML} and \code{mzML} refs. 
}

\keyword{methods}

