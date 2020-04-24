\name{initWF}
\alias{initWF}
\alias{SYSproject}
\alias{initProject}
\title{
Workflow Project Initiation
}
\description{
Collection of functions to construct \code{SYSargsList} workflow control environment 
(S4 object). 
\code{initWF} function constructs \code{SYSargsList} S4 class objects from a script file. 
The latter is an \code{RMarkdown} file containing all the workflow steps. In the future, 
it will be accepting other formats, as \code{Rscript} and \code{bash}. This file 
is defined either by the \code{sysconfig} file or \code{script} argument. 
Internally, \code{initWF} function calls \code{initProject} function, which creates 
a storage folder for all the logs, and the output to stdout or stderr. This function 
also creates the \code{SYSconfig.yml} on the fly to register all the parameters of 
the workflow. This function checks the project structure, which by default is:

\itemize{
  \item{\code{SYSproject/}}
    \itemize{
      \item{\code{data/}}
      \item{\code{param/}}
      \item{\code{results/}}
  }
}

The default structure of the project directories can be modified here. 
}
\usage{
initWF(sysconfig = NULL, subProj = FALSE, dir_str = "level0", dirName = "default", 
script, targets = NULL, silent = FALSE, overwrite = FALSE)

initProject(projPath = "./", data = "data", param = "param", results = "results", 
script = NULL, targets = NULL, subProj = FALSE, filename = "SYSconfig.yml", 
overwrite = FALSE, silent = FALSE)

SYSproject(projPath = "./", overwrite = FALSE, silent = FALSE, tempdir = FALSE)
}
\arguments{
  \item{sysconfig}{
  a character vector of a full \code{sysconfig} file path name. Default is \code{NULL}, 
which the function understand that the \code{Rmd} will be provide.
}
  \item{subProj}{
  is set TRUE, the project will be initialized and a new directory will be generated,
  at the same level as the current project, maintaining the directory structure.
  It will create an empty \code{results} folder, coping the \code{param} from the 
  current project, and creating a symbolic link for the \code{data} folder.
}
  \item{dir_str}{
  if \code{subProj} is set TRUE, it will create a directory on the same level of
  the current location.Default is \code{level0}.
}
  \item{dirName}{
  if \code{subProj} is set TRUE, the directory name can be defined here. 
}
  \item{script}{
  a character vector of a full \code{script} file path name.
This file containing all the information to run all the steps on the workflow analysis.
The current format accepted is \code{Rmd}.
}
  \item{targets}{
  a character vector of a full \code{sample annotation file} path name. Default is \code{NULL}.
The targets file can be either a simple tabular or yml file. Also, it is possible 
to assign NULL to run the pipeline without the 'targets' file. This can be useful
for running specific workflows that do not require input files.
}
  \item{silent}{
if set to TRUE, all messages returned by the function will be suppressed.
}
  \item{overwrite}{
  if set to TRUE, existing files or directories of same name will be overwritten.
}
  \item{projPath}{
  a character vector of a full project path name. Default is the current path. 
}
  \item{data}{
  a character vector of a \code{data} directory name. Default is \code{data}. 
This subdirectory in the project stores all the raw data, reference, and annotation files.
}
  \item{param}{
  a character vector of a \code{param} directory  name. Default is \code{param}.
This subdirectory in the project stores all the parameter and configuration files.
}
  \item{results}{
  a character vector of a \code{results} directory name. Default is \code{results}. 
This subdirectory in the project stores all the analysis results, including but 
not limited: alignment, variant, and peak files (BAM, VCF, BED); tabular result files; 
and image/plot files.
}
  \item{filename}{
  name and path of the output file which will store all the project configuration information. Default is \code{SYSconfig.yml}.
}
  \item{tempdir}{
if set to TRUE, the temporary directory for the workflow management files will be created.
}
}
\details{
\code{SYSproject} function detects and creates the project directory to store all 
the flat files from the workflow management. If \code{overwrite} is TRUE, the previous 
folder will be deleted and replace for a new one. Also, the user can choose to store 
all this information on a temporary folder, setting \code{tempdir} = TRUE. 

\code{initProject} function detects an existing project or creates a project 
structure on the path provide.
}
\value{
\code{SYSproject} will return a character vector containing the path name.
\code{initProject} will return a list containing all the information in the 
  \code{sysconfig} file.
\code{initWF} will return a \code{SYSargsList} object.
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
See also as \code{SYSargsList-class}.
}
\examples{
\dontrun{
## Construct SYSargsList object from Rmd file
library(systemPipeRdata)
script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- tempdir()
SYSconfig <- initProject(projPath=dir_path, targets=targets, script=script, overwrite = TRUE)
sysargslist <- initWF(sysconfig ="SYSconfig.yml")
}
}
\keyword{ utilities }
