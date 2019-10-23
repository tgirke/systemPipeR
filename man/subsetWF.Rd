\name{subsetWF}
\alias{subsetWF}

\title{
Subsetting SYSargs2 class slots  
}
\description{
Return subsets of \code{character} for the \code{input}, \code{output} or the list of command-line for each workflow \code{step}.
}
\usage{
subsetWF(args, slot, subset=NULL, index=NULL, delete=FALSE)
}
\arguments{
  \item{args}{
	object of class \code{SYSargs2}.
}
  \item{slot}{
  three options available: \code{type="input"} returns input slot from \code{SYSargs2} object; \code{type="output"} returns output slot from \code{SYSargs2} object; and \code{type="step"} returns all the command-line for each workflow step from \code{SYSargs2} object.
}
  \item{subset}{
  name or numeric position of the values to be subsetting in the \code{slot}. If \code{slot="input"}, the \code{subset} are the variables defined in the \code{param.yml} file. If \code{slot="step"}, the \code{subset} is the command line defined on the \code{SYSargs2} object for all the steps of the workflow. If \code{slot="output"}, the \code{subset} is the path for the expected output files for all the steps in the workflow. Default is \code{subset=NULL}
}
  \item{index}{
  A numeric index positions of the file in \code{SYSargs2} object, slot \code{output}. It requires a \code{subset} to be defined. Default is \code{index=NULL}.
}
  \item{delete}{
  allows to delete a subset of files in the case of \code{slot="output"}. Default is \code{delete=NULL}.
}
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{loadWorkflow}
\code{renderWF}
}
\examples{
## Construct SYSargs2 object 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", 
                  input_file="hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF

## Testing subset_wf function
input <- subsetWF(WF, slot="input", subset='FileName')
output <- subsetWF(WF, slot="output", subset=1, index=1)
step.cmd <- subsetWF(WF, slot="step", subset=1) ## subset all the HISAT2 commandline 
# subsetWF(WF, slot="output", subset=1, index=1, delete=TRUE) ## in order to delete the subset files list
}
\keyword{ utilities }