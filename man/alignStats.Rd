\name{alignStats}
\alias{alignStats}
\title{
Alignment statistics 
}
\description{
Generate data frame containing important read alignment statistics such as the
total number of reads in the FASTQ files, the number of total alignments, as
well as the number of primary alignments in the corresponding BAM files. 
}
\usage{
alignStats(args, output_index = 1)
}
\arguments{
  \item{args}{
	Object of class \code{SYSargs} or \code{SYSargs2}. 
}
  \item{index}{
  A numeric index positions of the file in \code{SYSargs2} object, slot \code{output}. Default is \code{output_index=1}.
}
}
\value{
	\code{data.frame} with alignment statistics.
}
\seealso{
\code{clusterRun} and \code{runCommandline} and \code{output_update}
}
\author{
Thomas Girke
}
\examples{
#########################################
## Examples with \code{SYSargs} object ##
#########################################
## Construct SYSargs object from param and targets files 
param <- system.file("extdata", "tophat.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args <- systemArgs(sysma=param, mytargets=targets)
args
names(args); modules(args); cores(args); outpaths(args); sysargs(args)

\dontrun{
## Execute SYSargs on single machine
runCommandline(args=args)

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
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", 
                  input_file="hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
targets <- system.file("extdata", "targets.txt", package="systemPipeR")

names(WF); modules(WF); targets(WF)[1]; cmdlist(WF)[1:2]; output(WF)

\dontrun{
## Execute SYSargs2 on single machine
WF <- runCommandline(args=WF)

## Alignment stats
read_statsDF <- alignStats(WF) 
read_statsDF <- cbind(read_statsDF[targets$FileName,], targets)
write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
}
}
\keyword{ utilities }