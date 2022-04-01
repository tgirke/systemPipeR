\name{SYSargsList}
% Constructor
\alias{SYSargsList}

\title{SYSargsList Constructor}
\description{
\code{SYSargsList} instances are constructed by the \code{SYSargsList} function.
}
\usage{
SYSargsList(sysargs = NULL, step_name = "default", 
            targets = NULL, wf_file = NULL, input_file = NULL, dir_path = ".", 
            id = "SampleName",
            inputvars = NULL, rm_targets_col = NULL, 
            dir = TRUE, dependency = NA, 
            run_step = "mandatory",
            run_session = "management",
            run_remote_resources = NULL,
            silent = FALSE, 
            projPath = getOption("projPath", getwd()))
}

\arguments{
  \item{sysargs}{
  \code{SYSargs2} object. If the object already exists, it can be used to 
  construct the \code{SYSargsList} object. 
}
  \item{step_name}{
  character with the step index name.
}
  \item{targets}{
  the path to \code{targets} file. The \code{targets} file can be either a simple 
  tabular or yml file. Also, it is possible to assign \code{NULL} to run the 
  pipeline without the 'targets' file. This can be useful for running specific 
  workflows that do not require input files.
}
  \item{wf_file}{
  name and path to \code{CWL parameters} file.
}
  \item{input_file}{
  name and path to \code{input parameters} file.
}
  \item{dir_path}{
  full path to the directory with the \code{CWL parameters} and \code{input} 
  files.
}
  \item{id}{
  A column from \code{targets} file, which will be used as an id for each one
  of the samples. It is required to be unique.
}
  \item{inputvars}{
  Each vector element is required to be defined in the \code{input} file, and the 
  names of the elements are needed to match the column names defined in the 
  \code{targets} file.
}
  \item{rm_targets_col}{
  targets file colunms to be removed.
}
  \item{dir}{
  This option allows creating an exclusive results folder for each step in the 
  workflow. All the outfiles and log files for the particular step will be 
  created in the respective folders. Default is \code{dir=TRUE}.
}
  \item{dependency}{
  character. Dependency tree, required when appending this step to 
  the workflow. Character name of a previous step in the workflow. 
  Default is \code{NA}.
}
  \item{run_step}{
  character. If the step has "mandatory" or "optional" flag for the execution.
}
  \item{run_session}{
  character. If the step has "management" or "compute" flag for the execution.
} 
  \item{run_remote_resources}{
  \code{List} for reserving for each cluster job sufficient computing resources
  including memory (Megabyte), number of nodes (Int), CPU cores, walltime (Minutes),
  etc. It is necessary two additional files: \code{conffile} and \code{template}. 
  \code{conffile} is the path to configuration file (default location 
  \code{./.batchtools.conf.R}). This file contains in its simplest form just one 
  command, such as this line for the scheduler: 
  \code{cluster.functions <- makeClusterFunctionsSlurm(template="batchtools.slurm.tmpl")}. 
  For more detailed information visit this page: https://mllg.github.io/batchtools/index.html
  \code{template} The template files for a specific queueing/scheduling systems 
  can be downloaded from here: 
  https://github.com/mllg/batchtools/tree/master/inst/templates. Slurm, 
  PBS/Torque, and Sun Grid Engine (SGE) templates are provided within the package
  demo data.
}
  \item{silent}{
  If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{projPath}{
  a character vector of a full project path name. Default is the current path. 
}
}
\author{
Daniela Cassol
}
\seealso{
\code{SYSargs2}, \code{LineWise}, and \code{SPRproject}
}
\examples{
sal <- SPRproject(overwrite=TRUE)
targetspath <- system.file("extdata/cwl/example/targets_example.txt", 
                            package="systemPipeR")

## Constructor and `appendStep<-`
appendStep(sal) <- SYSargsList(step_name = "echo", 
                               targets=targetspath, dir=TRUE,
                               wf_file="example/workflow_example.cwl", 
                               input_file="example/example.yml", 
                               dir_path = system.file("extdata/cwl", 
                               package="systemPipeR"),
                               inputvars = c(Message = "_STRING_",
                               SampleName = "_SAMPLE_"))
                               
appendStep(sal) <- LineWise(code = {
                            hello <- lapply(getColumn(sal, step=1, 'outfiles'), function(x) yaml::read_yaml(x))
                            }, 
                            step_name = "R_read", 
                            dependency = "echo")
                            
sal
}
\keyword{ utilities }

