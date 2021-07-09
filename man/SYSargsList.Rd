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
SYSargsList(args = NULL, step_name = "default", 
            targets = NULL, wf_file = NULL, input_file = NULL, dir_path = ".", 
            inputvars = NULL, rm_targets_col = NULL, 
            dir = TRUE, dependency = NA, silent = FALSE)
}
\arguments{
  \item{args}{
  \code{SYSargs2} object
}
  \item{step_name}{
  character with the step index name.
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
  \item{rm_targets_col}{
  targets file colunms to be removed.
}
  \item{dir}{
  This option allows creating an exclusive results folder for each step in the 
  workflow and a sub-folder for each sample defined in the targets file. All the 
  outputs and log files for the particular step will be created in the respective
  folders. Default is \code{dir=TRUE}.
}
  \item{dependency}{
  step previous dependency.
}
  \item{silent}{
  If set to \code{TRUE}, all messages returned by the function will be suppressed.
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


