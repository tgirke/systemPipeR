\name{SYSargs2list}
\alias{SYSargs2list}
\alias{SYSargs2list-method}
\alias{targets}
\alias{targets-methods}
\alias{wf}
\alias{wf-methods}
\alias{clt}
\alias{clt-methods}
\alias{yamlinput}
\alias{yamlinput-methods}
\alias{cmdlist}
\alias{cmdlist-methods}
\alias{input}
\alias{input-methods}
\alias{output}
\alias{output-methods}
\alias{cwlfiles}
\alias{cwlfiles-methods}
\alias{inputvars}
\alias{inputvars-methods}
\docType{data}
\title{ SYSargs2 accessor methods }                                  
\description{
	Methods to access information from \code{SYSargs2} object.
}
\usage{
	SYSargs2list(x)
}
\arguments{
  \item{x}{
	object of class \code{SYSargs2}
}
}
\value{
	various outputs
}
\author{
Daniela Cassol and Thomas Girke
}
\examples{
## Construct SYSargs2 object from CWl param, CWL input, and targets files 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", 
                  input_file="hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
names(WF); modules(WF); targets(WF)[1]; cmdlist(WF)[1:2]; output(WF)
}
\keyword{ utilities }
