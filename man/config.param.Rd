\name{config.param}
\alias{config.param}
\title{
Adding param file
}
\description{
Replace or adding a input configuration setting at "YML" param file
}
\usage{
config.param(input_file = NULL, param, file = "default", silent = FALSE)
}
\arguments{
  \item{input_file}{
  a \code{list} of parameters, \code{param file} path, or \code{SYSargs2} object.
}
  \item{param}{
  object of class \code{list}, expressing the values and names to be added or 
  replace at the command-line definition or the configuration files.
}
  \item{file}{
   name and path of the new file. If set to \code{default}, the name of the 
   write file will have the pattern: \code{<Date-Time>_<OriginalName>.yml}.
   If set to \code{append}, the \code{param} information it will be append on the 
   same file.
}
  \item{silent}{
  if set to TRUE, all messages returned by the function will be suppressed.
}
}
\author{
Daniela Cassol
}
\examples{
\dontrun{
input_file <- system.file("extdata", "cwl/hisat2/hisat2-se/hisat2-mapping-se.yml", package="systemPipeR")
param <- list(thread=10, fq=list(class="File", path="./results2"))
input <- config.param(input_file=input_file, param, file="default")
 }
}
\keyword{ utilities }
