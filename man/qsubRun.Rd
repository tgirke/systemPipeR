\name{qsubRun}
\alias{qsubRun}
\title{
Submit command-line tools to cluster
}
\description{
Submits command-line tools to queue (e.g. Torque) of compute cluster using run specifications defined by \code{runX} and \code{getQsubargs} functions. 
}
\usage{
qsubRun(appfct="runCommandline(args=args, runid='01')", args, qsubargs, Nqsubs = 1, package = "systemPipeR", shebang="#!/bin/bash")
}
\arguments{
  \item{appfct}{
	Accpets runX functions, such as \code{appfct="runCommandline(args, runid)"}
}
  \item{args}{
	Argument list returned by \code{systemArgs()}.
}
  \item{qsubargs}{
	Argument list returned by \code{getQsubargs()}.
}
  \item{Nqsubs}{
	Interger defining the number of qsub processes. Note: the function will not assign more qsub processes than there are FASTQ files. E.g. if there are 10 FASTQ files and \code{Nqsubs=20} then the function will generate only 10 qsub processes. To increase the number of CPU cores used by each process, one can increase the \code{p} value under \code{systemArgs()}.
}
  \item{package}{
	Package to load. Name provided as character vector of length one. Default is \code{sytemPipeR}.
}
  \item{shebang}{
	defines \code{shebang} (fist line) used in submission shell script; default is set to \code{#!/bin/bash}.
}
}
\value{
Returns list where list components contain FASTQ file names and their names are the qsub process IDs assiged by the queuing system. In addition, three files will be generated for each qsub submission process: \code{submitargs0X} (R object containing appargs), \code{submitargs0X.R} (R script using appargs) and \code{submitargs0X.sh} (shell submission script). In addition, the chosen \code{runX} function will output a \code{submitargs0X_log} file for each qsub process containing the executable commands processed by each qsub instance.
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
