\name{SPRproject}
\alias{SPRproject}
\title{
Workflow Project Initiation
}
\description{
Function to construct \code{SYSargsList} workflow control environment (S4 object). 
This function creates and checks the directory structure. If the expected directories
are not available, it is possible to create those. 
The project directory default structure expected is: 

\itemize{
  \item{\code{SPRproject/}}
    \itemize{
      \item{\code{data/}}
      \item{\code{param/}}
      \item{\code{results/}}
  }
}
The project working directory names can be modified, but users need to edit
the code accordingly.   
}
\usage{
SPRproject(projPath = getwd(), data = "data", param = "param", results = "results",
                      logs.dir= ".SPRproject", sys.file="SYSargsList.yml", 
                      envir = new.env(), 
                      restart = FALSE, resume=FALSE, 
                      load.envir = FALSE,
                      overwrite = FALSE, silent = FALSE)
}
\arguments{
  \item{projPath}{
  a character vector of a full project path name. Default is the current path. 
}
  \item{data}{
  a character vector of a \code{data} directory name. Default is \code{data}. 
This subdirectory in the project stores all the raw data, reference, and annotation files.
}
  \item{param}{
  a character vector of a \code{param} directory name. Default is \code{param}.
This subdirectory in the project stores all the parameter and configuration files.
}
  \item{results}{
  a character vector of a \code{results} directory name. Default is \code{results}. 
This subdirectory in the project stores all the analysis results, including but 
not limited: alignment, variant, and peak files (BAM, VCF, BED); tabular result files; 
and image/plot files.
}
  \item{logs.dir}{
  a character vector of a \code{logs} directory name. Default is \code{.SPRproject}.
}
  \item{sys.file}{
  a character vector of \code{SYSargsList} file name which will store all the 
  project configuration information. 
  Default is \code{SYSargsList.yml} and it will be save in the \code{logs.dir}
  folder.
}
  \item{envir}{
  the environment in which workflow will be evaluated. Default will create
  a \code{new.env()}. 
}
  \item{restart}{
  if set to \code{TRUE}, existing \code{SYSargsList} object saved in the 
  \code{logs.dir} directory will be used to restart the workflow. This option
  will delete all the log files. 
}
  \item{resume}{
  if set to \code{TRUE}, existing \code{SYSargsList} object saved in the 
  \code{logs.dir} directory will be used to resume the workflow.
}
  \item{load.envir}{
  This argument allows to load the environment and recover all the objects saved
  during the workflow execution (please check \code{runWF} function for more 
  details). This argument can be set as \code{TRUE} when \code{restart = TRUE} 
  or \code{resume = TRUE} are used.
}
  \item{overwrite}{
  if set to \code{TRUE}, existing \code{logs.dir} directory and all the content,
  as logs files and \code{SYSargsList} file will be removed, and a new and 
  empty \code{SYSargsList} object will be created. 
  This option should be used with caution!
}
  \item{silent}{
  if set to TRUE, all messages returned by the function will be suppressed.
}
}
\details{
  If an \code{SYSargsList} instance was created before or independent of the 
  project initialization, it is possible to append this instance after the 
  project is created. Please see check \code{appendStep<-} function.
}
\value{
\code{SPRproject} will return a \code{SYSargsList} object.
}
\author{
Daniela Cassol 
}
\seealso{
See also as \code{SYSargsList-class}.
}
\examples{
sal <- SPRproject(projPath = tempdir())
sal
}
\keyword{ utilities }
