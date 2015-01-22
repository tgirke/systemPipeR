\name{sysargs}
\alias{targetsin}
\alias{targetsin-methods}
\alias{targetsout}
\alias{targetsout-methods}
\alias{targetsheader}
\alias{targetsheader-methods}
\alias{modules}
\alias{modules-methods}
\alias{software}
\alias{software-methods}
\alias{cores}
\alias{cores-methods}
\alias{other}
\alias{other-methods}
\alias{reference}
\alias{reference-methods}
\alias{results}
\alias{results-methods}
\alias{infile1}
\alias{infile1-methods}
\alias{infile2}
\alias{infile2-methods}
\alias{outfile1}
\alias{outfile1-methods}
\alias{SampleName}
\alias{SampleName-methods}
\alias{sysargs}
\alias{sysargs-methods}
\alias{outpaths}
\alias{outpaths-methods}
\docType{data}
\title{ SYSargs accessor methods}

\description{
	Methods to access information from \code{SYSargs} object.
}
\usage{
	sysargs(x)
}
\arguments{
  \item{x}{
	object of class \code{SYSargs}
}
}
\value{
	various outputs
}
\author{
Thomas Girke
}
\examples{
## Construct SYSargs object from param and targets files 
param <- system.file("extdata", "tophat.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args <- systemArgs(sysma=param, mytargets=targets)
args
names(args); modules(args); cores(args); outpaths(args); sysargs(args)
}
\keyword{ utilities }

