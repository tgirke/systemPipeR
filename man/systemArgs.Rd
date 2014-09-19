\name{systemArgs}
\alias{systemArgs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructs SYSargs object from param and targets files
}
\description{
Constructs \code{SYSargs} object from 'param' and 'targets' files.
The 'param' file sets the parameters of the command-line software
and the 'targets' file defines the samples (e.g. FASTQ files).  
}
\usage{
systemArgs(sysma, mytargets, type = "SYSargs")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sysma}{ path to 'param' file; file structure follows a simple name/value syntax that
	converted into JSON format; for details about the file structure see sample files provided
	by package.
}
  \item{mytargets}{ path to targets file
}
  \item{type}{ \code{type="SYSargs"} returns \code{SYSargs}, \code{type="json"} returns param file 
		content in JSON format (requires \code{rjson} library)
}
}
\value{
\code{SYSargs} object or \code{character} string in JSON format
}
\author{
Thomas Girke
}
\seealso{
\code{showClass("SYSargs")}
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
qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=1", cores=cores(args), memory="mem=10gb", time="walltime=20:00:00")
qsubRun(appfct="runCommandline(args=args)", appargs=args, qsubargs=qsubargs, Nqsubs=1, submitdir="results", package="systemPipeR")
## Alignment stats
read_statsDF <- alignStats(fqpaths=infile1(args), bampaths=outpaths(args), fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{utilities}
