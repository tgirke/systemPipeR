\name{tryCMD}
\alias{tryCMD}
\alias{tryCL}
\title{
Collect information about the third-party software
}
\description{
Function to check if third-party software or utility is installed and set in the PATH.
}
\usage{
tryCMD(command)
}
\arguments{
  \item{command}{
a character vector containing the command line name to be tested.
}
}
\value{
It will return a positive message if the software is set on the PATH or an error message if the software is not set it. 
}
\note{
Please note that not necessary the software is not installed if the message indicates an error,  but it has not been exported on the current PATH.
}

\author{
Danela Cassol
}
\examples{
tryCMD(command="R") 
\dontrun{
tryCMD(command="blastp") 
tryCMD(command="hisat2") 
  }
}
\keyword{ utilities }