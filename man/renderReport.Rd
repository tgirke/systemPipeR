\name{renderLogs}
\alias{renderLogs}
\title{
Render RMarkdown Logs Report
}
\description{
Render the logs report file to the specified output format using pandoc.
}
\usage{
renderLogs(sysargs, type = "html_document", fileName = "default", silent = FALSE)
}
\arguments{
  \item{sysargs}{
  object of class \code{SYSargsList}.
}
  \item{type}{
  	The R Markdown output format to convert to. Default is \code{"ALL"}.The option 
  	can be the name of a format (e.g. \code{"pdf_document"} or \code{"html_document"}).
}
  \item{fileName}{
  haracter string naming a file output. Default is \code{"logs_<date>.Rmd"}.
}
  \item{silent}{
If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
}
\value{
It will return an \code{SYSargsList} updated.
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
See also as \code{SYSargsList-class}.
}
\examples{
## Construct SYSargsList object from Rmd file
sal <- SPRproject(overwrite=TRUE)
targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")

## Constructor and `appendStep<-`
appendStep(sal) <- SYSargsList(step_name = "echo", 
                               targets=targetspath, dir=TRUE,
                               wf_file="example/workflow_example.cwl", input_file="example/example.yml", 
                               dir_path = system.file("extdata/cwl", package="systemPipeR"),
                               inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
appendStep(sal) <- LineWise(code = {
                            hello <- lapply(getColumn(sal, step=1, 'outfiles'), function(x) yaml::read_yaml(x))
                            }, 
                            step_name = "R_read", 
                            dependency = "echo")
\dontrun{
sal <- runWF(sal)
sal <- renderLogs(sal)
}
}
\keyword{ utilities }