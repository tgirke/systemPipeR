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
	by package. To run the pipeline without 'param' file, \code{NULL} can be assigned. This 
        can be userful to run partial workflow, e.g. with pregenerated bam files.
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

## Execute SYSargs on multiple machines of a compute cluster 
resources <- list(walltime="00:25:00", nodes=paste0("1:ppn=", cores(args)), memory="2gb")
reg <- clusterRun(args, conffile=".BatchJobs.R", template="torque.tmpl", Njobs=18, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
showStatus(reg)
file.exists(outpaths(args))
sapply(1:length(args), function(x) loadResult(reg, x)) # Works once all jobs have completed successfully.

## Alignment stats
read_statsDF <- alignStats(fqpaths=infile1(args), bampaths=outpaths(args), fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{utilities}
