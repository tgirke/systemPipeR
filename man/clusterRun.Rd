\name{clusterRun}
\alias{clusterRun}
\title{
Submit command-line tools to cluster
}
\description{
Submits non-R command-line software to queueing/scheduling systems of compute clusters using run specifications defined by functions similar to \code{runCommandline}. \code{runCluster} can be used with most queueing systems since it is based on utilities from the \code{batchtools} package which supports the use of template files (\code{*.tmpl}) for defining the run parameters of the different schedulers. The path to the \code{*.tmpl} file needs to be specified in a conf file provided under the \code{conf.file} argument.
}
\usage{
clusterRun(args, FUN=runCommandline, conf.file=".batchtools.conf.R", template="batchtools.slurm.tmpl", runid="01",resourceList)   
}
\arguments{
  \item{args}{
	Object of class \code{SYSargs}. 
}
  \item{FUN}{
	Accpets functions such as \code{runCommandline(args, ...)} where the \code{args} argument is mandatory and needs to be of class \code{SYSargs}.
}
  \item{conf.file}{
	Path to conf file (default location \code{./.batchtools.conf.R}). This file contains in its simplest form just one command, such as this line for the Slurm scheduler: 
        \code{cluster.functions <- makeClusterFunctionsSlurm(template="batchtools.slurm.tmpl")}.
For more detailed information visit this page: https://mllg.github.io/batchtools/index.html
}
  \item{template}{
	The template files for a specific queueing/scheduling systems can be downloaded from here: 
https://github.com/mllg/batchtools/tree/master/inst/templates
}
  \item{runid}{
	Run identifier used for log file to track system call commands. Default is \code{"01"}.
}
  \item{resourceList}{
	\code{List} for reserving for each cluster job sufficient computing resources including memory (Megabyte), number of nodes, CPU cores, walltime (minutes), e-mail, etc. For more details, one can consult the template file for each queueing/scheduling system. 
}
}
\value{Object of class \code{Registry}, as well as files and directories created by the executed command-line tools.
}
\references{
For more details on \code{batchtools}, please consult the following page:
https://github.com/mllg/batchtools/
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
## read the instructions above how to obtain the corresponding files for other schedulers. 
file.copy(system.file("extdata", ".batchtools.conf.R", package="systemPipeR"), ".")
file.copy(system.file("extdata", "batchtools.slurm.tmpl", package="systemPipeR"), ".")
resources <- list(walltime=1200, ntasks=1, ncpus=cores(args), memory=10240)  
reg <- clusterRun(args, conf.file=".batchtools.conf.R", template="batchtools.slurm.tmpl", runid="01",resourceList=resources)
waitForJobs(reg=reg) 

## Monitor progress of submitted jobs
getStatus(reg=reg)
file.exists(outpaths(args))
sapply(1:length(args), function(x) loadResult(reg, x)) # Works once all jobs have completed successfully.

## Alignment stats
read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=TRUE) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }
