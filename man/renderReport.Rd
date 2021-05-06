\name{renderReport}
\alias{renderReport}
\title{
Render RMarkdown Report
}
\description{
Render the input file to the specified output format using pandoc.
}
\usage{
renderReport(sysargslist, type = c("html_document"),  silent = FALSE)
}
\arguments{
  \item{sysargslist}{
  object of class SYSargsList.
}
  \item{type}{
  	The R Markdown output format to convert to. Default is \code{"ALL"}.The option 
  	can be the name of a format (e.g. \code{"pdf_document"} or \code{"html_document"}).
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
\dontrun{
## Construct SYSargsList object from Rmd file
library(systemPipeRdata)
script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- tempdir()
SYSconfig <- initProject(projPath=dir_path, targets=targets, script=script, overwrite = TRUE)
sysargslist <- initWF(sysconfig ="SYSconfig.yml")
sysargslist <- configWF(x=sysargslist, input_steps = "1:3")
sysargslist <- runWF(sysargslist = sysargslist, steps = "ALL")
renderReport(x=sysargslist)
}
}
\keyword{ utilities }