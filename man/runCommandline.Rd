\name{runCommandline}
\alias{runCommandline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Execute SYSargs and SYSargs2
}
\description{
Function to execute system parameters specified in \code{SYSargs} and \code{SYSargs2} object.
}
\usage{
runCommandline(args, runid = "01", make_bam = FALSE, del_sam=TRUE, dir = TRUE, dir.name = NULL, force=FALSE, ...)
}
\arguments{
  \item{args}{
	object of class \code{SYSargs} or \code{SYSargs2}.
}
  \item{runid}{
	Run identifier used for log file to track system call commands. Default is \code{"01"}.
}
  \item{make_bam}{
	Auto detects SAM file outputs and converts them to sorted and indexed BAM files. Default is \code{make_bam=FALSE}.
}
  \item{del_sam}{
	This option allows deleting the SAM files created when the \code{make_BAM} converts the SAM files to sorted and indexed BAM files. Default is \code{del_sam=TRUE}.
}
  \item{dir}{
  This option allows creating an exclusive results folder for each step in the workflow and a sub-folder for each sample defined in the \code{targets} file. All the outputs and log files for the particular step will be created in the respective folders. Default is \code{dir=TRUE}. Option available only for an object of class \code{SYSargs2}.
}
  \item{dir.name}{
   Name of the workflow directory. Default is \code{dir.name=FALSE}. Note: This argument is required when the \code{dir=TRUE}.
}
  \item{force}{
   Internally, the function checks if the expected \code{output} files exist, and it skips the command lines when the respective files exist. If the argument \code{force} is set to \code{TRUE}, the command line will be executed and the files overwrite. Default is \code{force=FALSE}. 
}
  \item{\dots}{
	Additional arguments to pass on to \code{runCommandline()}.
}
}
\value{
	Output files, their paths can be obtained with \code{outpaths()} from \code{SYSargs} container or \code{output()} from \code{SYSargs2}. In addition, a \code{character vector} is returned containing the same paths.
}
\author{
Daniela Cassol and Thomas Girke
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
## Execute SYSargs on single machine
runCommandline(args=args)

## Execute SYSargs on multiple machines of a compute cluster. 
file.copy(system.file("extdata", ".batchtools.conf.R", package="systemPipeR"), ".")
file.copy(system.file("extdata", "batchtools.slurm.tmpl", package="systemPipeR"), ".")
resources <- list(walltime=120, ntasks=1, ncpus=cores(args), memory=1024) 
reg <- clusterRun(args, FUN = runCommandline, conffile=".batchtools.conf.R", template="batchtools.slurm.tmpl", Njobs=18, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
getStatus(reg=reg)
file.exists(outpaths(args))

## Alignment stats
read_statsDF <- alignStats(args) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}

##########################################
## Examples with \code{SYSargs2} object ##
##########################################
## Construct SYSargs2 object from CWl param, CWL input, and targets files 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2/hisat2-mapping-se.cwl", 
                  input_file="hisat2/hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
names(WF); modules(WF); targets(WF)[1]; cmdlist(WF)[1:2]; output(WF)

\dontrun{
## Execute SYSargs2 on single machine
WF <- runCommandline(args=WF)

## Execute SYSargs on multiple machines of a compute cluster. 
file.copy(system.file("extdata", ".batchtools.conf.R", package="systemPipeR"), ".")
file.copy(system.file("extdata", "batchtools.slurm.tmpl", package="systemPipeR"), ".")
resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024) 
reg <- clusterRun(WF, FUN = runCommandline, more.args = list(args = WF, make_bam = TRUE), conffile=".batchtools.conf.R", template="batchtools.slurm.tmpl", Njobs=18, runid="01", resourceList=resources)

## Monitor progress of submitted jobs
getStatus(reg=reg)

## Updates the path in the object \code{output(WF)}
WF <- output_update(WF, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam"))

## Alignment stats
read_statsDF <- alignStats(WF) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }