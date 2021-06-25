\name{runWF}
\alias{runWF}
\title{
Execute SYSargsList
}
\description{
Function to execute all the code list specified in SYSargsList object.
}
\usage{
runWF(args, dir = TRUE, dir.name = NULL, force = FALSE,
      warning.stop=FALSE, error.stop=TRUE, silent=FALSE, ...)
}
\arguments{
  \item{args}{
  object of class \code{SYSargsList}.
}
  \item{dir}{
  This option allows creating an exclusive results folder for each step in the workflow and a sub-folder for each sample defined in the \code{targets} file. All the outputs and log files for the particular step will be created in the respective folders. Default is \code{dir=TRUE}. Option available only for an object of class \code{SYSargs2}.
}
  \item{dir.name}{
   Name of the workflow directory. Default is \code{dir.name=FALSE}. Note: This argument is required when the \code{dir=TRUE}.
}
  \item{force}{
   Internally, the function checks if the expected \code{output} files exist, and it skips the command lines when the respective files exist. If the argument \code{force} is set to \code{TRUE}, the command line will be executed and the files overwrite. Default is \code{force=FALSE}. 
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
  \item{\dots}{
	Additional arguments to pass on from \code{runCommandline()}.
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
sysargslist <- SPRproject()
sal <- runWF(sal)
}
}
\keyword{ utilities }
