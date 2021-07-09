\name{LineWise}
\alias{LineWise}
\title{
LineWise
}
\description{
S4 class container for storing parameters of R-based software.
}
\usage{
LineWise(code, step_name="default", codeChunkStart=integer(), rmdPath=character(), 
        dependency=NA, status = list())
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
  \item{step_name}{
Step Name
}
  \item{dependency}{
dependency
}
  \item{status}{
status
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
LineWise({
          1+1
          },
          step_name="NewStep")
}
\keyword{ utilities }