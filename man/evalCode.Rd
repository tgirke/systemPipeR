\name{evalCode}
\alias{evalCode}
\title{
Toggles option \code{eval} on the RMarkdown files 
}
\description{
Function to evaluate (\code{eval=TRUE}) or not evaluate (\code{eval=FALSE}) R chunk codes in the Rmarkdown file. This function does not run the code, just toggles between \code{TRUE} or \code{FALSE} the option \code{eval} and writes out a new file with the chosen option.
}
\usage{
evalCode(infile, eval = TRUE, output)
}
\arguments{
  \item{infile}{
name and path of the infile file, format \code{Rmd}.
}
  \item{eval}{
whether to evaluate the code and include its results. The default is \code{TRUE}. 
}
  \item{output}{
name and path of the output file. File format \code{Rmd}.
}
}
\value{
Writes Rmarkdown file containing the exact copy of the \code{infile} file with the option choose on the \code{eval} argument. 
It will be easy to toggle between run all the R chunk codes or not. 
}
\author{
Daniela Cassol
}
\examples{
\dontrun{
library(systemPipeRdata)
file <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
evalCode(infile=file, eval=FALSE, output="test.Rmd")
 }
}
\keyword{ utilities }
