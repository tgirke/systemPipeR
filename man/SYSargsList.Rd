\name{SYSargsList}
\alias{SYSargsList}
\title{
SYSargsList
}
\description{
S4 class container for storing parameters of command-line- or R-based software. 
\code{SYSargsList} instances are constructed by the \code{SYSargsList} function.
}
\usage{
SYSargsList(args = NULL, targets = NULL, wf_file = NULL, input_file = NULL, 
dir_path = ".", inputvars = NULL, step_name = "default", silent = FALSE,  
rm_targets_col = NULL, dependency = NULL)
}
\arguments{
  \item{args}{
  \code{SYSargs2} object
}
  \item{targets}{
  the path to \code{targets} file. The \code{targets} file can be either a simple tabular or yml file. Also, it is possible to assign \code{NULL} to run the pipeline without the 'targets' file. This can be useful for running specific workflows that do not require input files.
}
  \item{wf_file}{
  name and path to \code{CWL param} file.
}
  \item{input_file}{
  name and path to \code{input} file.
}
  \item{dir_path}{
  full path to the directory with the \code{CWL param} and \code{input} files.
}
  \item{inputvars}{
  variables list defined in the \code{input} file that matches the column names defined in the \code{targets} file.
}
  \item{step_name}{
  character with the step index name.
}
  \item{silent}{
  If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{rm_targets_col}{
  targets file colunms to be removed.
}
  \item{dependency}{
  step previous dependency.
}
}
\value{
\code{SYSargsList} object
}
\author{
Daniela Cassol
}
\seealso{
\code{showClass("SYSargsList")}
}
\examples{
dir_path <- system.file("extdata/cwl/example/", package="systemPipeR")
sal <- SYSargsList(wf_file="example.cwl", input_file="example_single.yml", dir_path = dir_path)
sal
cmdlist(sal)
}
\keyword{ utilities }


