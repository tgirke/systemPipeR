\name{check.output}
\alias{check.output}
\alias{check.outfiles}
\title{
Checking if the outfiles files exist
}
\description{
This function returns a \code{data.frame} indicating the number of existing files
and how many files are missing.
}
\usage{
check.output(sysargs, type="data.frame")
check.outfiles(sysargs, type="data.frame")
}

\arguments{
  \item{sysargs}{
    object of class \code{SYSargs2} or \code{SYSargsList}.
}
  \item{type}{
    return object option. It can be \code{data.frame} or \code{list}.
}
}
\value{
\code{data.frame} or \code{list} containing all the \code{outfiles} file information. 
}
\seealso{
- \code{\link{SYSargs2-class}}

- \code{\link{SYSargsList-class}}
}
\author{
Daniela Cassol and Thomas Girke
}
\examples{
## Construct SYSargs2 object 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2/hisat2-mapping-se.cwl", 
                  input_file="hisat2/hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
## Check output
check.output(WF)
check.output(WF, "list")

## Construct SYSargsList object 
sal <- SPRproject(overwrite=TRUE) 
targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")
appendStep(sal) <- SYSargsList(step_name = "echo", 
                      targets=targetspath, dir=TRUE,
                      wf_file="example/workflow_example.cwl", input_file="example/example.yml", 
                      dir_path = system.file("extdata/cwl", package="systemPipeR"),
                      inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
## Check outfiles
check.outfiles(sal)
}
\keyword{ utilities }
