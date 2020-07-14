\name{moduleload}
\alias{moduleload}
\alias{modulelist}
\alias{moduleAvail}
\alias{moduleClear}
\alias{moduleInit}
\alias{moduleUnload}
\title{
Interface to module system
}
\description{
Functions to list and load software from a module system in R. The functions are
the equivalent of \code{module avail} and \code{module load} on the Linux command-line, respectively.
}
\usage{
moduleload(module_name)
moduleUnload(module_name)
modulelist()
moduleAvail()
moduleClear()
moduleInit()

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{module_name}{
	Name of software to load character vector. 
}
}
\author{
Tyler Backman and Thomas Girke
}
\examples{
\dontrun{
## List all software from module system
modulelist()

## Examples for loading software from module system
moduleload(module="bowtie2/2.0.6", envir="PATH")
moduleload(module="python", envir=c("PATH", "LD_LIBRARY_PATH", "PYTHONPATH"))
}
}
\keyword{ utilities }
