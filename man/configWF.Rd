\name{configWF}
\alias{configWF}
\title{
Workflow Steps Selection
}
\description{
This function allows the user to control of which step of the workflow will be 
run and the generation of the new RMarkdown.
}
\usage{
configWF(x = sysargslist, input_steps = "ALL", exclude_steps = NULL, silent = FALSE, ...)
}
\arguments{
  \item{x}{
  object of class \code{SYSargsList}.
}
  \item{input_steps}{
  a character vector of all steps desires to preserve on the output file. 
Default is \code{ALL}. It can be used the symbol ":" to select many steps in 
sequence, for example, input_steps=1:5.2, from step 1 to step 5.2.
The symbol "." represents the substeps and symbol "," is used to separate selections.
Jump from a major step to sub-step is supported, but if a major step is selected/excluded, 
all sub-steps of this major step will be selected/excluded. Repeatedly selected 
steps will only result in a unique step.It is recommended to put major steps in 
\code{input_steps}, like \code{1:4, 6:8, 10}; and unwanted sub-steps in 
\code{exclude_step}, like \code{1.1, 3.1.1-3.1.3, 6.5}. Reverse selecting is 
supported e.g. \code{10:1}. 
}
  \item{exclude_steps}{
  a character vector of all steps desires to remove on the output file.
}
  \item{silent}{
  if set to TRUE, all messages returned by the function will be suppressed.
}
  \item{\dots}{
  Additional arguments to pass on to \code{configWF()}.
}
}
\author{
Daniela Cassol
}
\examples{
\dontrun{
library(systemPipeRdata)
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
SYSconfig <- initProject(projPath="./", targets=targets, script=script, overwrite=TRUE, silent=TRUE)

sysargslist <- configWF(x=sysargslist)
 }
}
\keyword{ utilities }
