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
runCommandline(args, runid = "01", make_bam=TRUE, ...)
}
\arguments{
  \item{args}{
	object of class \code{SYSargs}
}
  \item{runid}{
	Run identifier used for log file to track system call commands. Default is \code{"01"}.
}
  \item{make_bam}{
	Autodetects SAM file outputs and converts them to sorted and indexed BAM files. Default is \code{make_bam=TRUE}. 
}
  \item{\dots}{
	Additional arguments to pass on to \code{runCommandline()}.
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
param <- system.file("extdata", "hisat2.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args <- systemArgs(sysma=param, mytargets=targets)
args
names(args); modules(args); cores(args); outpaths(args); sysargs(args)

\dontrun{
## Execute SYSargs on single machine
runCommandline(args=args)

## Execute SYSargs on multiple machines of a compute cluster 
resources <- list(walltime=120, ntasks=1, ncpus=cores(args), memory=1024) 
reg <- clusterRun(args, conffile=".batchtools.conf.R", template="batchtools.slurm.tmpl", Njobs=18, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
getStatus(reg=reg)
file.exists(outpaths(args))
sapply(1:length(args), function(x) loadResult(reg, x)) # Works once all jobs have completed successfully.

## Alignment stats
read_statsDF <- alignStats(args) 
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }

