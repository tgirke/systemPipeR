\name{targets.as.df}
\alias{targets.as.df}
\alias{targets.as.list}
\title{
Convert targets to list or data.frame
}
\description{
Convert targets files to \code{list} or \code{data.frame} object. 
}
\usage{
targets.as.df(x)
targets.as.list(x)
}
\arguments{
  \item{x}{
  An object of the class \code{"list"} or \code{"data.frame"} that stores data 
  from each target file, as \code{targets(WF)}.
  }
}
\value{
\code{data.frame} or \code{list} containing all the \code{targets} file information. 
}
\seealso{
\code{showClass("SYSargs2")}
}
\author{
Daniela Cassol
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targetslist <- targets.as.list(x=targets)
targets.as.df(x=targetslist)
}
\keyword{ utilities }
