\name{importWF}
\alias{importWF}
\title{
Import RMarkdown file as workflow
}
\description{
Import RMarkdown file as workflow
}
\usage{
importWF(sal, file_path, ignore_eval = TRUE, verbose = TRUE, ...)
}

\arguments{
  \item{sal}{
  SYSargsList empty object
}
  \item{file_path}{
  string, file path of the workflow file
}
  \item{ignore_eval}{
  logical, treat all R chunks' `eval` option as `TRUE` in workflow Rmd file even if some chunks have `eval=FALSE`
}
  \item{verbose}{
  logical, print out verbose message while function running
}
  \item{...}{
  Additional arguments to pass on to \code{SPRproject()}.
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
rmdpath <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
#importWF(rmdpath)
}
\keyword{ utilities }
