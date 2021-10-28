\name{runWF}
\alias{runWF}
\title{
Execute SYSargsList
}
\description{
Function to execute all the code list specified in SYSargsList object.
}
\usage{
runWF(sysargs, steps = NULL, targets = NULL,
                  force = FALSE, saveEnv = TRUE,
                  run_step = "ALL", ignore.dep = FALSE,
                  warning.stop = FALSE, error.stop = TRUE, silent = FALSE, ...)
}
\arguments{
  \item{sysargs}{
    object of class \code{SYSargsList}.
}
  \item{steps}{
    character or numeric. Step name or index. If \code{NULL}, all the steps
    will be executed.
}
  \item{targets}{
    This option allows selecting which targets file and, by consequence which 
    command line will be executed for each \code{SYSargs2} class step.
    Default is \code{NULL}, in which all command lines will be executed.
    
}
  \item{force}{
    Internally, the option checks if the expected \code{output} files exist, 
    and it skips the command lines when the respective files exist. 
    If the argument \code{force} is set to \code{TRUE}, the command line will 
    be executed and the files overwrite. Default is \code{force=FALSE}. 
}
  \item{saveEnv}{
    If set to \code{TRUE}, the environment will be saved to an RDS file. 
    To check the RDS file location, please use 
    \code{projectInfo(sysargs)[["envir"]]}.
}
  \item{run_step}{
    character. If the step has "mandatory" or "optional" flag for the execution.
    When \code{ALL}, all the steps will be executed. 
}
  \item{ignore.dep}{
    logical. This option allow to igonore the dependency tree, when \code{TRUE}.
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
\dontrun{
sal <- runWF(sal)
}
}
\keyword{ utilities }
