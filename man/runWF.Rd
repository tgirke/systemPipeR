\name{runWF}
\alias{runWF}
\title{
Execute SYSargsList
}
\description{
Function to execute all the code list specified in SYSargsList object.
}
\usage{
runWF(args, warning.stop=FALSE, error.stop=TRUE, silent=FALSE)
}
\arguments{
  \item{args}{
  object of class \code{SYSargsList}.
}
  \item{warning.stop}{
If set to \code{TRUE}, the process will be interrupted when a warning is detected.
}
  \item{error.stop}{
If set to \code{TRUE}, the process will be interrupted when a error is detected.
}
  \item{silent}{
If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
}
\value{
It will return an \code{SYSargsList} updated.
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
See also as \code{SYSargsList-class}.
}
\examples{
\dontrun{
## Construct SYSargsList object from Rmd file
library(systemPipeRdata)
script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- tempdir()
SYSconfig <- initProject(projPath=dir_path, targets=targets, script=script, overwrite = TRUE)
sysargslist <- initWF(sprconfig ="SPRconfig.yml")
sysargslist <- configWF(x=sysargslist, input_steps = "1:3")
sal <- runWF(args = sysargslist)
}
}
\keyword{ utilities }
