\name{sal2rmd}
\alias{sal2rmd}
\title{
Translate SYSargsList back to a workflow template Rmarkdown file
}
\description{
This function takes a SYSargsList object and translate it to SPR workflow 
template Rmarkdown format.
}
\usage{
sal2rmd(sal, out_path = "spr_template.Rmd", rmd_title = "SPR workflow template", 
        rmd_author = "my name", 
        rmd_date = "Last update: `r format(Sys.time(), '\%d \%B, \%Y')`", 
        rmd_output = "html_document", 
        desc = "This is a workflow template.", verbose = TRUE)
}
\arguments{
  \item{sal}{
  \code{SYSargsList} object.
}
  \item{out_path}{
   string, output file name.
}
  \item{rmd_title}{
  string, title of the Rmd
}
  \item{rmd_author}{
  string, author(s) of the Rmd, put all authors in a single character string
}
  \item{rmd_date}{
  string, date header of Rmd
}
  \item{rmd_output}{
  string, output format of Rmd, used in header
}
  \item{desc}{
  string, or character vector of strings, some description text in format 
  Rmarkdown that will be added to the document before the workflow steps start. 
  It can be a single line or multiple lines by providing a character vector,
  each item is one line.
}
  \item{verbose}{
  bool, open verbose mode? This will show you more information as 
  the function runs 
}
}
\value{
    no return
}

\author{
    Le Zhang and Daniela Cassol
}
\examples{
file_path <- system.file("extdata/spr_simple_wf.Rmd", package="systemPipeR")
sal <- SPRproject(overwrite = TRUE)
sal <- importWF(sal, file_path)
sal2rmd(sal)
}
\keyword{ utilities }
