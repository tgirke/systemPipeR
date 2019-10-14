\name{createWF}
\alias{createWF}
\title{
Constructs SYSargs2 object and creates CWL param files
}
\description{
The constructor function creates an \code{SYSargs2} S4 class object from three input files: a \code{CWL param} and \code{input} files, and one simple tabular files, a \code{targets} file. The latter is optional for workflow steps lacking input files. Also, the function creates and saves the\code{CWL param} files. The latest storages all the parameters required for running command-line software, following the standard and specification defined on \href{https://www.commonwl.org/}{Common Workflow Language (CWL)}.
}
\usage{
createWF(targets = NULL, commandLine, results_path = "./results", module_load = "baseCommand", file = "default", overwrite = FALSE, cwlVersion = "v1.0", class = "CommandLineTool")
}
\arguments{
  \item{targets}{
  Path to \code{targets} file. Assign \code{NULL} to run the pipeline without 'targets' file. This can be useful for running specific workflows which do not require input files.
}
  \item{commandLine}{
Object of class \code{list}, which provides all the parameters required for running command-line software.
}
  \item{results_path}{
Path to the \code{results} folder.
}
  \item{module_load}{
Name of software to load as \code{character}. Default is \code{"default"}, which creates a subfolder and two files: \code{*.cwl} and \code{*.yml} at \code{./param/cwl/}.
}
  \item{file}{
Default is \code{"default"}, which creates a subfolder and two files: \code{*.cwl} and \code{*.yml} at \code{./param/cwl/}. If a different location and names are required, the names of the files can be specified as a character vector. 
}
  \item{overwrite}{
If set to \code{TRUE}, existing files of the same name will be overwritten.
}
  \item{cwlVersion}{
version of the Common Workflow Language. More information here: https://www.commonwl.org/.
}
  \item{class}{
Names of Common Workflow Language Specification. The following switches are accepted: \code{CommandLineTool} and \code{Workflow}.
}
}
\value{
\code{SYSargs2} object
}
\references{
For more details on \code{CWL}, please consult the following page: https://www.commonwl.org/
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{loadWorkflow}
\code{renderWF}
\code{showClass("SYSargs2")}
}
\examples{
# "hisat2 -S ./results/_SampleName_.sam  -x ./data/tair10.fasta  -k 1  --min-intronlen 30  --max-intronlen 3000 --threads 4 -U _FASTQ_PATH1_ "
## Provide a list with all the arguments
baseCommand <- list(baseCommand="hisat2")
inputs <- list(
  "S"=list(type="File", preF="-S", yml="./results/_SampleName_.sam"),
  "x"=list(type="File", preF="-x", yml="./data/tair10.fasta"),
  "k"= list(type="int", preF="-k", yml= 1L),
  "threads"= list(type="int", preF="-threads", yml=4L),
  "min-intronlen"= list(type="int", preF="--min-intronlen", yml= 30L),
  "max-intronlen"= list(type="int", preF="--max-intronlen", yml=3000L),
  "U"=list(type="File", preF="-U", yml="./data/_FASTQ_PATH1_") )
outputs <- list("hisat2_sam"=list(type="File", "./results/_SampleName_.sam"))
commandLine <- list(baseCommand=baseCommand, inputs=inputs, outputs=outputs)
\dontrun{
## Create a SYSargs2 object and populate all the command-line
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
WF <- createWF(targets=targets, commandLine, results_path="./results", module_load="baseCommand", file = "default", overwrite = FALSE, cwlVersion = "v1.0", class = "CommandLineTool")
WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
  }
}
\keyword{ utilities }
