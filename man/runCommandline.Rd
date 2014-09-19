\name{runCommandline}
\alias{runCommandline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Execute SYSargs
}
\description{
Function to execute system parameters specified in \code{SYSargs} object
}
\usage{
runCommandline(args, runid = "01", usemodule=TRUE)
}
\arguments{
  \item{args}{
	object of class \code{SYSargs}
}
  \item{runid}{
	Run identifier used for log file to track system call commands. Default is \code{"01"}.
}
  \item{usemodule}{
	Determines whether the external executable software should be loaded from a module system. If set to \code{TRUE} then the software version specified in the \code{*.param} file will be used. The function \code{modules} returns this parameter from a \code{SYSargs} instance. On systems where a module system is not available, \code{usemodule} should be set to \code{FALSE}.
}
}
\value{
	Output files, their paths can be obtained with \code{outpaths()} from \code{SYSargs} container. 
	In addition, a \code{character vector} is returned containing the same paths.
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

\dontrun{
## Execute SYSargs on single machine
runCommandline(args=args)

## Execute SYSargs on multiple machines  
qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=cores(tophat), memory="mem=10gb", time="walltime=20:00:00")
qsubRun(args=args, qsubargs=qsubargs, Nqsubs=1, package="systemPipeR")
## Alignment stats
read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }

