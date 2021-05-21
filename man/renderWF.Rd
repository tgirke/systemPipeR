\name{renderWF}
\alias{renderWF}
\title{
Populate all the command-line in an \code{SYSargs2} object
}
\description{
The \code{SYSargs2} S4 class object is constructed from the \code{loadWorkflow}, which stores all the information and instructions needed for processing a set of input files with a specific command-line or a series of command-line within a workflow. The \code{renderWF} function populates all the command-line for each sample in each step of the particular workflow. Each sample level input/outfile operation uses its own \code{SYSargs2} instance. The output of \code{SYSargs2} define all the expected output files for each step in the workflow, which usually it is the sample input for the next step in an \code{SYSargs2} instance. Between different instances, this connectivity is established by writing the subsetting output with the \code{writeTargetsout} function to a new targets file that serves as input to the next \code{loadWorkflow} and \code{renderWF} call. By chaining several \code{SYSargs2} steps together one can construct complex workflows involving many sample-level input/output file operations with any combination of command-line or R-based software. 
}
\usage{
renderWF(WF, inputvars = NULL)
}
\arguments{
  \item{WF}{
  Object of class \code{SYSargs2}.
}
  \item{inputvars}{
  variables list defined in the \code{input} file that matches the column names defined in the \code{targets} file.
}
}
\value{
\code{SYSargs2} object
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{showClass("SYSargs2")}
\code{loadWorkflow}
\code{writeTargetsout}
}
\examples{
## Construct SYSargs2 object from CWl param, CWL input, and targets files 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", 
                  input_file="hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
}
\keyword{ utilities }