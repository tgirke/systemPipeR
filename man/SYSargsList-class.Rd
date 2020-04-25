\name{SYSargsList-class}
\Rdversion{1.1}
\docType{class}
\alias{SYSargsList-class}
\alias{[,SYSargsList,ANY,ANY,ANY-method}
\alias{[[,SYSargsList,ANY,ANY-method}
\alias{[[<-,SYSargsList,ANY,ANY,ANY-method}
\alias{$,SYSargsList-method}
\alias{codeSteps,SYSargsList-method}
\alias{coerce,list,SYSargsList-method}
\alias{coerce,SYSargsList,list-method}
\alias{length,SYSargsList-method}
\alias{names,SYSargsList-method}
\alias{projectWF,SYSargsList-method}
\alias{show,SYSargsList-method}
\alias{statusWF,SYSargsList-method}
\alias{dataWF,SYSargsList-method}
\alias{stepsWF,SYSargsList-method}
\alias{summaryWF,SYSargsList-method}
\alias{SYSargs2_steps,SYSargsList-method}
\alias{sysargslist,SYSargsList-method}
\alias{sysconfig,SYSargsList-method}

\title{Class \code{"SYSargsList"}}
\description{
S4 class container for storing parameters and infrmation of command-line- or R-based workflow.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SYSargsList", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{sysconfig}:}{Object of class \code{"list"} storing all the config information for the project }
    \item{\code{codeSteps}:}{Object of class \code{"list"} storing all the \code{code steps} of the workflow }
    \item{\code{stepsWF}:}{Object of class \code{"numeric"} storing all the index steps levels in the workflow }
    \item{\code{dataWF}:}{Object of class \code{"data.frame"} storing all the data steps in the workflow }
    \item{\code{SYSargs2_steps}:}{Object of class \code{"list"} storing all the \code{SYSargs2} objects of the workflow }
    \item{\code{statusWF}:}{Object of class \code{"list"} storing all the success and failure of each step in the workflow}
    \item{\code{projectWF}:}{Object of class \code{"list"} storing all the directories from the project }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY", drop = "ANY")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY")}: subsetting of class with bracket operator }
    \item{[[<-}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY", value = "ANY")}: replacement method for SYSargs2 class }
    \item{$}{\code{signature(x = "SYSargsList")}: extracting slots elements by name }
    \item{codeSteps}{\code{signature(x = "SYSargsList")}: ... }
    \item{coerce}{\code{signature(from = "list", to = "SYSargsList")}: \code{as(list, "SYSargsList")} }
    \item{coerce}{\code{signature(from = "SYSargsList", to = "list")}: \code{as(SYSargsList, "list")} }
    \item{length}{\code{signature(x = "SYSargsList")}: extracts number of \code{SYSargsList} objects }
    \item{names}{\code{signature(x = "SYSargsList")}: extracts slot names }
    \item{projectWF}{\code{signature(x = "SYSargsList")}: extract data from \code{projectWF} slot }
    \item{show}{\code{signature(object = "SYSargsList")}: summary view of \code{SYSargsList} steps }
    \item{statusWF}{\code{signature(x = "SYSargsList")}: extract data from \code{statusWF} slot }
    \item{dataWF}{\code{signature(x = "SYSargsList")}: extract data from \code{dataWF} slot }
    \item{stepsWF}{\code{signature(x = "SYSargsList")}: extract data from \code{stepsWF} slot }
    \item{SYSargs2_steps}{\code{signature(x = "SYSargsList")}: extract data from \code{SYSargs2_steps} slot }
    \item{sysargslist}{\code{signature(x = "SYSargsList")}: Coerce back to list \code{as(SYSargsList, "list")} }
    \item{sysconfig}{\code{signature(x = "SYSargsList")}: extract data from \code{sysconfig} slot }
  }
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{SYSargs2} and \code{initWF}
}
\examples{
showClass("SYSargsList")
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
\keyword{ classes }

