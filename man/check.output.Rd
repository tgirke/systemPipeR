\name{check.output}
\alias{check.output}
\title{
Checking if the output files exist
}
\description{
This function returns a \code{data.frame} indicating the number of existing files
and how many files are missing.
}
\usage{
check.output(args)
}

\arguments{
  \item{args}{
	object of class \code{SYSargs2} or \code{SYSargsList}.
}
}
\author{
Daniela Cassol and Thomas Girke
}
\examples{
## Construct SYSargs2 object 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", 
                  input_file="hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF

check.output(WF)
}
\keyword{ utilities }
