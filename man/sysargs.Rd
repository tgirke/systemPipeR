\name{sysargs}
\alias{targetsin}
\alias{targetsin-methods}
\alias{targetsout}
\alias{targetsout-methods}
%\alias{targetsheader}
\alias{targetsheader-methods}
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
%\alias{SampleName}
\alias{SampleName-methods}
\alias{sysargs}
\alias{sysargs-methods}
\alias{outpaths}
\alias{outpaths-methods}
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
	The content of the slot of the same name in \code{x}:
	\code{targetsin} and \code{targetsout} return a \code{data.frame},
	\code{cores} a \code{numeric}, and \code{modules}, \code{software},
	\code{other}, \code{reference}, \code{results}, \code{infile1},
	\code{infile2}, \code{outfile1}, \code{sysargs} and \code{outpaths} a
	\code{character} vector. \code{SampleName} returns the sample names, taken
	from the names of the \code{sysargs} slot.
}
\author{
Thomas Girke
}
\examples{
## Construct SYSargs object from param and targets files 
param <- system.file("extdata", "hisat2.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args <- systemArgs(sysma=param, mytargets=targets)
args
names(args); modules(args); cores(args); outpaths(args); sysargs(args)
}
\keyword{ utilities }

