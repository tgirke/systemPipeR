\name{SPRproject}
\alias{SPRproject}
\title{
Workflow Project Initiation
}
\description{
Function to construct \code{SYSargsList} workflow control environment (S4 object). 
This function creates the directory structure, with it is not available on the 
fly to register all the parameters of the workflow.
This function checks the project structure, which by default is:

\itemize{
  \item{\code{SPRproject/}}
    \itemize{
      \item{\code{data/}}
      \item{\code{param/}}
      \item{\code{results/}}
  }
}

The default structure of the project directories can be modified here. 
}
\usage{
SPRproject(projPath = getwd(), data = "data", param = "param", results = "results",
                      logs.dir= ".SPRproject",
                       restart = FALSE, sys.file=".SPRproject/SYSargsList.yml", silent = FALSE)
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
  a character vector of a \code{param} directory  name. Default is \code{param}.
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
  \item{restart}{
  if set to TRUE, existing SYSargsList save in the directory will be used to restart the 
  workflow.
}
  \item{sys.file}{
  name and path of the SYSargsList file which will store all the project configuration 
  information. Default is \code{.SPRproject/SYSargsList.yml}.
}
  \item{silent}{
if set to TRUE, all messages returned by the function will be suppressed.
}
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
## Construct SYSargsList object from Rmd file
library(systemPipeRdata)
sal <- SPRproject()
sal

}
\keyword{ utilities }
