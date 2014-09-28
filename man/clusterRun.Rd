\name{clusterRun}
\alias{clusterRun}
\title{
Submit command-line tools to cluster
}
\description{
Submits non-R command-line software to queueing/scheduling systems of compute clusters using run specifications defined by \code{runCommandline}. \code{runCluster} can be used with most queueing systems since it is based on utilities from the \code{BatchJobs} package which supports the use of template files (\code{*.tmpl}) for defining the run parameters of the different schedulers. The path to the \code{*.tmpl} file needs to be specified in a conf file provided under the \code{conffile} argument.
}
\usage{
clusterRun(args, conffile = ".BatchJobs.R", template = "torque.tmpl", Njobs, runid = "01", resourceList)
}
\arguments{
  \item{args}{
	Object of class \code{SYSargs}. 
}
  \item{conffile}{
	Path to conf file (default location \code{./.BatchJobs.R}). This file contains in its simplest form just one command, such as this line for the Torque scheduler: 
        \code{cluster.functions <- makeClusterFunctionsTorque("torque.tmpl")}.
For more detailed information visit this page: https://code.google.com/p/batchjobs/wiki/DortmundUsage
}
  \item{template}{
	The template files for a specific queueing/scheduling systems can be downloaded from here: 
https://github.com/tudo-r/BatchJobs/blob/master/examples/cfTorque/simple.tmpl
}
  \item{Njobs}{
	Interger defining the number of cluster jobs. For instance, if \code{args} contains 18 command-line jobs and \code{Njobs=9}, then the function will distribute them accross 9 cluster jobs each running 2 command-line jobs. To increase the number of CPU cores used by each process, one can do this under the corresonding argument of the command-line tool, e.g. \code{-p} argument for Tophat.
}
  \item{runid}{
	Run identifier used for log file to track system call commands. Default is \code{"01"}.
}
  \item{resourceList}{
	\code{List} for reserving for each cluster job sufficient computing resources including memory, number of nodes, CPU cores, walltime, etc. For more details, one can consult the template file for each queueing/scheduling system. 
}
}
\value{Object of class \code{Registry}, as well as files and directories created by the executed command-line tools.
}
\references{
For more details on \code{BatchJobs}, please consult the following pages:
http://sfb876.tu-dortmund.de/PublicPublicationFiles/bischl_etal_2012a.pdf
https://github.com/tudo-r/BatchJobs
http://goo.gl/k3Tu5Y
}
\author{
Thomas Girke
}
\seealso{
	\code{clusterRun} replaces the older functions \code{getQsubargs} and \code{qsubRun}.
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

## Execute SYSargs on multiple machines of a compute cluster. The following
## example uses the conf and template files for the Torque scheduler. Please
## read instructions how to obtain the corresponding files for other schedulers. 
file.copy(system.file("extdata", ".BatchJobs.R", package="systemPipeR"), ".")
file.copy(system.file("extdata", "torque.tmpl", package="systemPipeR"), ".")
resources <- list(walltime="00:25:00", nodes=paste0("1:ppn=", cores(args)), memory="2gb")
reg <- clusterRun(args, conffile=".BatchJobs", template="torque.tmpl", Njobs=18, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
showStatus(reg)
file.exists(outpaths(args))
sapply(1:length(args), function(x) loadResult(reg, x)) # Works once all jobs have completed successfully.

## Alignment stats
read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }
