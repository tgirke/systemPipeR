\name{LineWise}
\alias{LineWise}
\title{
LineWise
}
\description{
S4 class container for storing parameters of R-based software.
}
\usage{
LineWise(code, stepName="default", codeChunkStart=integer(), rmdPath=character(), dependency=character())
}
\arguments{
  \item{code}{
R code
}
  \item{codeChunkStart}{
R Markdown line start
}
  \item{rmdPath}{
path of R Markdown file
}
  \item{stepName}{
Step Name
}
  \item{dependency}{
dependency
}
}
\value{
\code{LineWise} object
}
\author{
Daniela Cassol
}
\seealso{
\code{showClass("LineWise")}
}
\examples{
LineWise("1+1", stepName="NewStep")
}
\keyword{ utilities }