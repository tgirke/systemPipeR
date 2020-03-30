\name{output_update}
\alias{output_update}

\title{
Updates the output files paths in the \code{SYSargs2} object
}
\description{
After executing all the command-lines by the \code{runCommadline} function, the output files can be created in specific directories rather then \code{results} in a particular directory. Also, the \code{runCommadline} function allows converting the SAM file outputs to sorted and indexed BAM files. Thus, the \code{output_update} function allows updating the location of these files in the output of the \code{SYSargs2}  object.
}
\usage{
output_update(args, dir = FALSE, dir.name = NULL, replace = FALSE, extension = NULL, make_bam=FALSE)
}
\arguments{
  \item{args}{
	object of class \code{SYSargs2}.
}
  \item{dir}{
  assign \code{TRUE} to update the location of the output files in the \code{args} object accordingly with the workflow name directory. Default is \code{dir=FALSE}.
}
  \item{dir.name}{
  if the results directory name is not specified in the \code{input} file, it is possible to specify here the name. This argument is required if the path name return \code{NULL} from the \code{input} file. Default is \code{dir.name=NULL}.
}
  \item{replace}{
  replace the extension for selected output files in the \code{args} object. Default is \code{replace=FALSE}.
}
  \item{extension}{
  object of class \code{"character"} storing the current extension of the files and the respective replacement. For example, \code{runCommandline} function by default autodetects SAM file outputs in the \code{args} object and create the BAM files. In order to update the output of \code{args} object, the extension argument should be set: extension = c(\code{".sam"}, \code{".bam"}).
}
  \item{make_bam}{
	Auto detects SAM file outputs and update them on the \code{SYSargs2} object for sorted and indexed BAM files. Default is \code{make_bam=FALSE}. This argument should be used in integration with \code{runCommandline} function. 
}
}
\value{
\code{SYSargs2} object with output location files updated.
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
To check directory name in the \code{input} file: \code{yamlinput(WF)$results_path$path}.
}
\examples{ 
## Construct SYSargs2 object from CWl param, CWL input, and targets files 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2-mapping-se.cwl", 
                  input_file="hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
output(WF)

\dontrun{
runCommandline(args=WF, make_bam=TRUE)
## Output paths update
WF <- output_update(WF, dir=FALSE, replace=TRUE, extension=c(".sam", ".bam"))

runCommandline(args=WF, make_bam=TRUE, dir=TRUE)
## Output paths update
WF <- output_update(WF, dir=TRUE, replace=TRUE, extension=c(".sam", ".bam"))
}
}
\keyword{ utilities }
