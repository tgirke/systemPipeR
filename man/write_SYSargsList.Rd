\name{write_SYSargsList}
\alias{write_SYSargsList}
\title{
Writeout SYSargsList object
}
\description{
Function to writeout \code{SYSargsList} workflow control environment (S4 object) 
object.

}
\usage{
write_SYSargsList(sysargs, sys.file=".SPRproject/SYSargsList.yml", silent=TRUE)
}
\arguments{
  \item{sysargs}{
  object of class \code{SYSargsList}.
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
\code{write_SYSargsList} will return a \code{sys.file} path.
}
\author{
Daniela Cassol 
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

sal <- write_SYSargsList(sal)
sal
}
\keyword{ utilities }
