\name{importRmd}
\alias{importRmd}
\title{
Import RMarkdown file as workflow
}
\description{
Import RMarkdown file as workflow
}
\usage{
importRmd(file_path, ignore_eval = TRUE, verbose = TRUE)
}

\arguments{
  \item{file_path}{
  string, file path of the workflow file
}
  \item{ignore_eval}{
  logical, treat all R chunks' `eval` option as `TRUE` in workflow Rmd file even if some chunks have `eval=FALSE`
}
  \item{verbose}{
  logical, print out verbose message while function running
}
}
\value{
data.frame
}
\author{
Le Zhang and Daniela Cassol
}
\seealso{
\code{\link{SYSargsList}}
}
\examples{
library(systemPipeRdata)
script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
#importRmd(script)
}
\keyword{ utilities }
