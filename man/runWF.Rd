\name{runWF}
\alias{runWF}
\title{
Execute SYSargsList
}
\description{
Function to execute all the code list specified in SYSargsList object.
}
\usage{
runWF(sysargs, steps = NULL, force = FALSE, saveEnv = TRUE, warning.stop =
                 FALSE, error.stop = TRUE, silent = FALSE, ...)
}
\arguments{
  \item{sysargs}{
  object of class \code{SYSargsList}.
}
  \item{steps}{
  character or numeric. Step name or index. If \code{NULL}, all the step will be 
  executed.
}
  \item{force}{
   Internally, the function checks if the expected \code{output} files exist, and it skips the command lines when the respective files exist. If the argument \code{force} is set to \code{TRUE}, the command line will be executed and the files overwrite. Default is \code{force=FALSE}. 
}
  \item{saveEnv}{
If set to \code{TRUE}, the environment will be saved to an RDS file. 
}
  \item{warning.stop}{
If set to \code{TRUE}, the process will be interrupted when a warning is detected.
}
  \item{error.stop}{
If set to \code{TRUE}, the process will be interrupted when a error is detected.
}
  \item{silent}{
If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{\dots}{
	Additional arguments to pass on from \code{runCommandline()}.
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
## Construct SYSargsList object from Rmd file
sal <- SPRproject(overwrite=TRUE)
targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")

## Constructor and `appendStep<-`
appendStep(sal) <- SYSargsList(step_name = "echo", 
                               targets=targetspath, dir=TRUE,
                               wf_file="example/workflow_example.cwl", input_file="example/example.yml", 
                               dir_path = system.file("extdata/cwl", package="systemPipeR"),
                               inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
appendStep(sal) <- LineWise(code = {
                            hello <- lapply(getColumn(sal, step=1, 'outfiles'), function(x) yaml::read_yaml(x))
                            }, 
                            step_name = "R_read", 
                            dependency = "echo")
sal <- runWF(sal)
}
\keyword{ utilities }
