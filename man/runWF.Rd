\name{runWF}
\alias{runWF}
\title{
Execute SYSargsList
}
\description{
Function to execute all the code list specified in SYSargsList object.
}
\usage{
runWF(sysargslist, steps = "ALL")
}
\arguments{
  \item{sysargslist}{
  object of class SYSargsList.
}
  \item{steps}{
  	a character vector of all steps desires to preserve on the output file. 
  	Default is \code{"ALL"}. It can be used the symbol ":" to select many steps in sequence, 
  	for example, input_steps=1:5.2, from step 1 to step 5.2. The symbol "." represents 
  	the substeps and symbol "," is used to separate selections. Jump from a 
  	major step to sub-step is supported, but if a major step is selected, all sub-steps
  	of this major step will be selected. Repeatedly selected steps will only result 
  	in a unique step. It is recommended to put major steps in input_steps, like 1:4, 6:8, 10; 
  	and unwanted sub-steps in exclude_step, like 1.1, 3.1.1-3.1.3, 6.5. Reverse selecting 
  	is supported e.g. 10:1.
}
}
\value{
It will return an \code{SYSargsList} updated.
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
sysargslist <- configWF(x=sysargslist, input_steps = "1:3")
sysargslist <- runWF(sysargslist = sysargslist, steps = "ALL")
sysargslist <- runWF(sysargslist = sysargslist, steps = "1:2")
}
}
\keyword{ utilities }
