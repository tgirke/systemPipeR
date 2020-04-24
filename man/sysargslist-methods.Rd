\name{sysargslist}
\docType{methods}
\alias{sysargslist-method}
\alias{sysargslist}
\alias{sysconfig}
\alias{sysconfig-methods}
\alias{codeSteps}
\alias{codeSteps-methods}
\alias{stepsWF}
\alias{stepsWF-methods}
\alias{SYSargs2_steps}
\alias{SYSargs2_steps-methods}
\alias{statusWF}
\alias{statusWF-methods}
\alias{projectWF}
\alias{projectWF-methods}
\title{ SYSargsList accessor methods }
\description{
	Methods to access information from \code{SYSargsList} object.
}
\usage{
	sysargslist(x)
}
\arguments{
  \item{x}{
	object of class \code{SYSargsList}.
}
}
\value{
	various outputs
}
\author{
Daniela Cassol and Thomas Girke
}
\examples{
\dontrun{
## Construct SYSargsList object from Rmd file
library(systemPipeRdata)
script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- tempdir()
SYSconfig <- initProject(projPath=dir_path, targets=targets, script=script, overwrite = TRUE)
sysargslist <- initWF(sysconfig ="SYSconfig.yml")
  }
}
\keyword{ utilities }