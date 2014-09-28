\name{getQsubargs}
\alias{getQsubargs}
\title{
Arguments for qsub 
}
\description{
Note: This function as been deprecated. Please use \code{clusterRun} instead. \code{getQsubargs} defines arguments to submit runX job(s) to queuing system (e.g. Torque) via \code{qsub}.
}
\usage{
getQsubargs(software = "qsub", queue = "batch", Nnodes = "nodes=1", cores = as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory = "mem=10gb", time = "walltime=20:00:00")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{software}{
	Software to use for submission to queuing system. Default is \code{qsub}.
}
  \item{queue}{
	Name of queue to to use. Default is \code{batch}.
}
  \item{Nnodes}{
	Number of compute nodes to use for processing. Default is \code{nodes=1}.
}
  \item{cores}{
	Number of CPU cores to use per compute node. Default will use what is provided by under \code{-p} in \code{myargs} of \code{systemArgs()} output.
}
  \item{memory}{
	Amount of RAM to reserve per node.
}
  \item{time}{
	Walltime limit each job is allowed to run per node.
}
}
\value{
	\code{list}
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
