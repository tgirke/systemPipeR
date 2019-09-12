\name{loadWorkflow}
\alias{loadWorkflow}
\title{
Constructs SYSargs2 object from CWL param and targets files
}
\description{
The constructor functions create an \code{SYSargs2} S4 class object from three input files: a \code{CWL param} and \code{input} files, and one simple tabular files, a \code{targets} file. The latter is optional for workflow steps lacking input files. The\code{CWL param} provides all the parameters required for running command-line software, following the standard and specification defined on \href{https://www.commonwl.org/}{Common Workflow Language (CWL)}. The \code{input} file provides additional information for the command-line, allowing each sample level input/outfile operation uses its own \code{SYSargs2} instance. In the \code{targets} file users could provide the paths to the initial sample input files (e.g. FASTQ) along with sample labels, and if appropriate biological replicate and contrast information for controlling differential abundance analyses.
}
\usage{
loadWorkflow(targets = NULL, wf_file, input_file, dir_path = ".")
}
\arguments{
  \item{targets}{
  path to \code{targets} file. Assign \code{NULL} to run the pipeline without 'targets' file. This can be useful for running specific workflows which do not require input files.
}
  \item{wf_file}{
  path to \code{CWL param} file.
}
  \item{input_file}{
  path to \code{input} file.
}
  \item{dir_path}{
  full path to the directory with the \code{CWL param} and \code{input} files.
}
}
\details{
Please note that the \code{CWL param} and \code{input} files need to be in the same directory.
}
\value{
\code{SYSargs2} object
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{renderWF}
\code{showClass("SYSargs2")}
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
