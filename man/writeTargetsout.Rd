\name{writeTargetsout}
\alias{writeTargetsout}
\title{
Write updated targets out to file
}
\description{
Convenience write function for generating targets files with updated
\code{FileName} columns containing the paths to files generated
by input/output processes. These processes can be commandline- or R-based software.
Typically, the paths to the inputs are stored in the targets infile 
(\code{targetsin(args)} for \code{SYSargs} objects or \code{targets.as.df(targets(WF))} 
for \code{SYSargs2} objects) and the outputs are stored in the targets outfile 
(\code{targetsout(args)} for \code{SYSargs} objects or \code{output(WF)}) for 
\code{SYSargs2} objects. Note: by default the function cannot overwrite 
any existing files. If a file exists then the user has to explicitly 
remove it or set \code{overwrite=TRUE}. 
}
\usage{
writeTargetsout(x, file = "default", silent = FALSE, overwrite = FALSE, step = NULL, new_col=NULL, new_col_output_index=NULL, remove=FALSE, ...)
}
\arguments{
  \item{x}{
Object of class \code{SYSargs} or \code{SYSargs2}.
}
  \item{file}{
Name and path of the output file. If set to "default" then the name of the output file will have the pattern 'targets_<software>.txt', where <software> will be what \code{software(x)} returns, when \code{x} is an object of class \code{SYSargs}. For an object of class \code{SYSargs2}, the output file will have the pattern 'targets_<software>_<step>.txt', where <software> will be the workflow name (\code{files(x)$cwl}) and <step> will be the step chosen in the argument \code{step}.
}
  \item{silent}{
If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{overwrite}{
If set to \code{TRUE}, existing files of same name will be overwritten.
}
  \item{step}{
Name or numeric position of the step in the workflow to write the output files. The names can be check running \code{names(x$clt)}.
}
  \item{new_col}{
A vector of character strings of the new column name to add to target file.
}
  \item{new_col_output_index}{
A vector of numeric index positions of the file in \code{SYSargs2} class \code{output}. It requires \code{new_col} definition.
}
  \item{remove}{
To remove the original columns. Default is \code{FALSE}. 
}
  \item{\dots}{
To pass on additional arguments.
}
}
\value{
Writes tabular targes files containing the header/comment lines from \code{targetsheader(x)}
and the columns from \code{targetsout(x)}. 
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{writeTargetsRef}
}
\examples{
##########################################
## Examples with \code{SYSargs2} object ##
##########################################
## Construct SYSargs2 object 
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl", package="systemPipeR")
WF <- loadWorkflow(targets=targets, wf_file="hisat2/hisat2-mapping-se.cwl", 
                  input_file="hisat2/hisat2-mapping-se.yml", dir_path=dir_path)
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
WF
## Write targets out file
names(WF$clt)
writeTargetsout(x=WF, file="default", step=1, new_col=c("sam_file"), new_col_output_index=c(1))

}
\keyword{ utilities }
