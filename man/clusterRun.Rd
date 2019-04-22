\name{clusterRun}
\alias{clusterRun}
\title{
Submit command-line tools to cluster
}
\description{
Submits non-R command-line software to queueing/scheduling systems of compute clusters using run specifications defined by functions similar to \code{runCommandline}. \code{clusterRun} can be used with most queueing systems since it is based on utilities from the \code{batchtools} package which supports the use of template files (\code{*.tmpl}) for defining the run parameters of the different schedulers. The path to the \code{*.tmpl} file needs to be specified in a conf file provided under the \code{conffile} argument.
}
\usage{
clusterRun(args, FUN = runCommandline, more.args = list(args = args, make_bam = TRUE), conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", Njobs, runid = "01", resourceList)
}
\arguments{
  \item{args}{
	Object of class \code{SYSargs} or \code{SYSargs2}.
}
  \item{FUN}{
	Accepts functions such as \code{runCommandline(args, ...)} where the \code{args} argument is mandatory and needs to be of class \code{SYSargs} or \code{SYSargs2}.
}
  \item{more.args}{
  Object of class \code{list}, which provides the arguments that control the \code{FUN} function. 
}
  \item{conffile}{
	Path to conf file (default location \code{./.batchtools.conf.R}). This file contains in its simplest form just one command, such as this line for the Slurm scheduler: \code{cluster.functions <- makeClusterFunctionsSlurm(template="batchtools.slurm.tmpl")}. For more detailed information visit this page: https://mllg.github.io/batchtools/index.html
}
  \item{template}{
	The template files for a specific queueing/scheduling systems can be downloaded from here: 
https://github.com/mllg/batchtools/tree/master/inst/templates. Slurm, PBS/Torque, and Sun Grid Engine (SGE) templates are provided.
}
  \item{Njobs}{
	Interger defining the number of cluster jobs. For instance, if \code{args} contains 18 command-line jobs and \code{Njobs=9}, then the function will distribute them accross 9 cluster jobs each running 2 command-line jobs. To increase the number of CPU cores used by each process, one can do this under the corresonding argument of the command-line tool, e.g. \code{-p} argument for Tophat.
}
  \item{runid}{
	Run identifier used for log file to track system call commands. Default is \code{"01"}.
}
  \item{resourceList}{
	\code{List} for reserving for each cluster job sufficient computing resources including memory (Megabyte), number of nodes, CPU cores, walltime (minutes), etc. For more details, one can consult the template file for each queueing/scheduling system.
}
}
\value{Object of class \code{Registry}, as well as files and directories created by the executed command-line tools.
}
\references{
For more details on \code{batchtools}, please consult the following page: https://github.com/mllg/batchtools/
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
	\code{clusterRun} replaces the older functions \code{getQsubargs} and \code{qsubRun}.
}
\examples{
#########################################
## Examples with \code{SYSargs} object ##
#########################################
## Construct SYSargs object from param and targets files 
param <- system.file("extdata", "hisat2.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args <- systemArgs(sysma=param, mytargets=targets)
args
names(args); modules(args); cores(args); outpaths(args); sysargs(args)

\dontrun{
## Execute SYSargs on multiple machines of a compute cluster. The following
## example uses the conf and template files for the Slurm scheduler. Please 
## read the instructions on how to obtain the corresponding files for other schedulers. 
file.copy(system.file("extdata", ".batchtools.conf.R", package="systemPipeR"), ".")
file.copy(system.file("extdata", "batchtools.slurm.tmpl", package="systemPipeR"), ".")
resources <- list(walltime=120, ntasks=1, ncpus=cores(args), memory=1024) 
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, make_bam = TRUE), conffile=".batchtools.conf.R", template="batchtools.slurm.tmpl", Njobs=18, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
getStatus(reg=reg)
file.exists(outpaths(args))
}

##########################################
## Examples with \code{SYSargs2} object ##
##########################################
## Construct SYSargs2 object from CWl param, CWL input, and targets files 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-se/hisat2-mapping-se.cwl",
                   input_file="hisat2-se/hisat2-mapping-se.yml",
                   dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH_", SampleName="_SampleName_"))
WF
names(WF); modules(WF); targets(WF)[1]; cmdlist(WF)[1:2]; output(WF)

\dontrun{
## Execute SYSargs2 on multiple machines of a compute cluster. The following
## example uses the conf and template files for the Slurm scheduler. Please 
## read the instructions on how to obtain the corresponding files for other schedulers.  
file.copy(system.file("extdata", ".batchtools.conf.R", package="systemPipeR"), ".")
file.copy(system.file("extdata", "batchtools.slurm.tmpl", package="systemPipeR"), ".")
resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024) 
reg <- clusterRun(WF, FUN = runCommandline, more.args = list(args = WF, make_bam = FALSE), conffile=".batchtools.conf.R", template="batchtools.slurm.tmpl", Njobs=18, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
getStatus(reg=reg)

## Updates the path in the object \code{output(WF)}
WF <- output_update(WF, dir=FALSE, replace = ".bam")

## Alignment stats
read_statsDF <- alignStats(WF) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }