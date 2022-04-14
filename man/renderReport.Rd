\name{renderReport}
\alias{renderReport}
\title{
Render RMarkdown Report
}
\description{
Render the technical report file to the specified output format using pandoc.
}
\usage{
renderReport(sysargs, fileName ="SPR_Report",
              rmd_title = "SPR workflow Template - Report",
              rmd_author = "Author",
              rmd_date= "Last update: `r format(Sys.time(), '\%d \%B, \%Y')`",
              type = c("html_document"),
              desc = "This is a workflow template.",
              quiet = FALSE, open_file = TRUE)
}
\arguments{
  \item{sysargs}{
  object of class \code{SYSargsList}.
}
  \item{fileName}{
  character string naming a file output. Default is \code{"spr_report.Rmd"}.
}
  \item{rmd_title}{
  string, title of the Rmd.
}
  \item{rmd_author}{
  string, author(s) of the Rmd, put all authors in a single character string.
}
  \item{rmd_date}{
  string, date header of Rmd.
}
  \item{type}{
  The R Markdown output format to convert to. The option 
  can be the name of a format (e.g. \code{"pdf_document"} or \code{"html_document"}).
}
  \item{desc}{
  string, or character vector of strings, some description text in format 
  Rmarkdown that will be added to the document before the workflow steps start. 
  It can be a single line or multiple lines by providing a character vector,
  each item is one line.
}
  \item{quiet}{
  If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{open_file}{
  Default is \code{TRUE}. 
}
}
\value{
It will return an \code{SYSargsList} updated, with the file path location.
}
\author{
Daniela Cassol 
}
\seealso{
See also as \code{SYSargsList-class}.
}
\examples{
sal <- SPRproject(overwrite = TRUE)
file_path <- system.file("extdata", "spr_simple_wf.Rmd", package = "systemPipeR")
sal <- importWF(sal, file_path = file_path, verbose = FALSE)
targetspath <- system.file("extdata/cwl/example/targets_example.txt", package = "systemPipeR")
appendStep(sal) <- SYSargsList(step_name = "echo", 
                              targets = targetspath, dir = TRUE,
                              wf_file = "example/workflow_example.cwl", input_file = "example/example.yml", 
                              dir_path = system.file("extdata/cwl", package = "systemPipeR"),
                              inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
sal <- renderReport(sal, open_file = FALSE)
}
\keyword{ utilities }