\name{SYSargs2Pipe_ls}
\alias{SYSargs2Pipe_ls}
\alias{SYSargs2Pipe-method}
\alias{WF_steps}
\alias{WF_steps-methods}
\alias{track}
\alias{track-methods}
\alias{summaryWF}
\alias{summaryWF-methods}
\docType{data}
\title{Deprecated \code{"SYSargs2Pipe"} accessor methods}
\description{
These methods are deprecated and should be replaced by \code{SYSargsList-methods}.
Methods to access information from \code{SYSargs2Pipe} object.
}
\usage{
	SYSargs2Pipe_ls(x)
}
\arguments{
  \item{x}{
	object of class \code{SYSargs2Pipe}.
}
}
\value{
	various outputs
}
\details{
The following methods are deprecated and will be made defunct; use the replacement indicated below:
  \itemize{
    \item{\code{SYSargs2Pipe}: \code{\link{SYSargs2-class}}}
  }
}
\author{
Daniela Cassol and Thomas Girke
}
\examples{
\dontrun{
## Construct SYSargs2 object number 1
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", 
                  input_file="hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF

## Construct SYSargs2 object number 2
targetsPE <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-pe", package="systemPipeR")
WF1 <- loadWorkflow(targets=targetsPE, wf_file="hisat2-mapping-pe.cwl", 
                  input_file="hisat2-mapping-pe.yml", dir_path=dir_path)
WF1 <- renderWF(WF1, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
WF1

## Keep track 
WF_set <- run_track(WF_ls = c(WF1, WF))
WF_steps(WF_set)
track(WF_set)
summaryWF(WF_set)[1]
}
}
\keyword{ utilities }
