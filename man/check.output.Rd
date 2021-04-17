\name{check.output}
\alias{check.output}
\title{
Checking if the output files exist
}
\description{
This function returns a logical vector indicating whether the expected output files 
for each workflow step exist. 
}
\usage{
check.output(args, subset = 1, index = 1)
}

\arguments{
  \item{args}{
	object of class \code{SYSargs2}.
}
  \item{subset}{
  name or numeric position of the values to be subsetting in the \code{slot}. The \code{subset} is the path for the expected output files for all the steps in the workflow. Default is \code{subset=1}
}
  \item{index}{
  A numeric index positions of the file in \code{SYSargs2} object, slot \code{output}. It requires a \code{subset} to be defined. Default is \code{index=1}.
}
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{subsetWF}
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
