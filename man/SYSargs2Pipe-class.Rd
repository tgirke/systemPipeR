\name{SYSargs2Pipe-class}
\Rdversion{1.1}
\docType{class}
\alias{SYSargs2Pipe-class}
\alias{[,SYSargs2Pipe,ANY,ANY,ANY-method}
\alias{[[,SYSargs2Pipe,ANY,ANY-method}
\alias{[[<-,SYSargs2Pipe,ANY,ANY,ANY-method}
\alias{$,SYSargs2Pipe-method}
\alias{coerce,list,SYSargs2Pipe-method}
\alias{coerce,SYSargs2Pipe,list-method}
\alias{length,SYSargs2Pipe-method}
\alias{names,SYSargs2Pipe-method}
\alias{show,SYSargs2Pipe-method}
\alias{summaryWF,SYSargs2Pipe-method}
\alias{SYSargs2Pipe_ls,SYSargs2Pipe-method}
\alias{track,SYSargs2Pipe-method}
\alias{WF_steps,SYSargs2Pipe-method}

\title{Class \code{"SYSargs2Pipe"}}
\description{
\code{SYSargs2Pipe class} stores a list \code{SYSargs2} objects. Each \code{SYSargs2} objects stores all the information and instructions needed for processing a set of input files with a specific command-line or a series of command-line within a workflow. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SYSargs2Pipe", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{WF_steps}:}{Object of class \code{"list"} storing all the \code{SYSargs2} objects }
    \item{\code{track}:}{Object of class \code{"list"} storing all the output files from each \code{SYSargs2} objects }
    \item{\code{summaryWF}:}{Object of class \code{"list"} storing the summary of all the expected files exists and how many were missing for each \code{SYSargs2} objects }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SYSargs2Pipe", i = "ANY", j = "ANY", drop = "ANY")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SYSargs2Pipe", i = "ANY", j = "ANY")}: subsetting of class with bracket operator }
    \item{[[<-}{\code{signature(x = "SYSargs2Pipe", i = "ANY", j = "ANY", value = "ANY")}: replacement method for SYSargs2 class }
    \item{$}{\code{signature(x = "SYSargs2Pipe")}: extracting slots elements by name }
    \item{coerce}{\code{signature(from = "list", to = "SYSargs2Pipe")}: \code{as(list, "SYSargs2Pipe")} }
    \item{coerce}{\code{signature(from = "SYSargs2Pipe", to = "list")}: \code{as(SYSargs2Pipe, "list")} }
    \item{length}{\code{signature(x = "SYSargs2Pipe")}: extracts number of \code{SYSargs2} objects }
    \item{names}{\code{signature(x = "SYSargs2Pipe")}: extracts slot names }
    \item{show}{\code{signature(object = "SYSargs2Pipe")}: summary view of \code{SYSargs2} objects }
    \item{summaryWF}{\code{signature(x = "SYSargs2Pipe")}: extract data from \code{targets} slot }
    \item{SYSargs2Pipe_ls}{\code{signature(x = "SYSargs2Pipe")}: Coerce back to list \code{as(SYSargs2Pipe, "list")} }
    \item{track}{\code{signature(x = "SYSargs2Pipe")}: extract data from \code{track} slot }
    \item{WF_steps}{\code{signature(x = "SYSargs2Pipe")}: extract data from \code{WF_steps} slot }
	 }
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{loadWorkflow} and \code{renderWF} and \code{runCommandline} and \code{clusterRun}
}
\examples{
showClass("SYSargs2Pipe")

## Construct SYSargs2 object from CWl param, CWL input, and targets files 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", 
                  input_file="hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF

## Keep track 
WF_set <- run_track(WF_ls = c(WF))
WF_steps(WF_set)
track(WF_set)
summaryWF(WF_set)[1]
}
\keyword{ classes }
