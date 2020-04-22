\name{tryPath}
\alias{tryPath}
\title{
Validation of the files or directories
}
\description{
Function to check if the full path (file or directory) exists. 
}
\usage{
tryPath(path)
}
\arguments{
  \item{path}{
a character vector of full path name.
}
}
\details{
This function produces a character vector of the file or directory name defined on the \code{path} argument.
}
\value{
A character vector containing the name of the file or directory. If the path does not exist, it will return an error message. 
}
\author{
Daniela Cassol
}
\examples{
file <- system.file("extdata/workflows/rnaseq/", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
tryPath(path=file)
}
\keyword{ utilities }