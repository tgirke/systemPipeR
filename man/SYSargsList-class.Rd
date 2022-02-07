\name{SYSargsList-class}
\Rdversion{1.1}
\docType{class}
\alias{SYSargsList-class}
% Accessors
\alias{stepsWF}
\alias{stepsWF,SYSargsList-method}
\alias{statusWF}
\alias{statusWF,SYSargsList-method}
\alias{targetsWF}
\alias{targetsWF,SYSargsList-method}
\alias{outfiles}
\alias{outfiles,SYSargsList-method}
\alias{SE}
\alias{SE,SYSargsList-method}
\alias{dependency}
\alias{dependency,SYSargsList-method}
\alias{targets_connection}
\alias{targets_connection,SYSargsList-method}
\alias{projectInfo}
\alias{projectInfo,SYSargsList-method}
\alias{runInfo}
\alias{runInfo,SYSargsList-method}
% Subsetting and Quick access
\alias{[,SYSargsList,ANY,ANY,ANY-method}
\alias{[[,SYSargsList,ANY,missing-method}
\alias{$,SYSargsList-method}
\alias{names,SYSargsList-method}
\alias{show,SYSargsList-method}
\alias{length,SYSargsList-method}
% Coercion
\alias{coerce,list,SYSargsList-method}
\alias{coerce,SYSargsList,list-method}
\alias{sysargslist}
\alias{sysargslist,SYSargsList-method}
% Methods 
\alias{subset}
\alias{subset,SYSargsList-method}
\alias{getColumn}
\alias{getColumn,SYSargsList-method}
\alias{baseCommand}
\alias{baseCommand,SYSargsList-method}
\alias{SampleName}
\alias{SampleName,SYSargsList-method}
\alias{stepName}
\alias{stepName,SYSargsList-method}
\alias{targetsheader}
\alias{targetsheader,SYSargsList-method}
\alias{yamlinput}
\alias{yamlinput,SYSargsList-method}
\alias{addResources}
\alias{addResources,SYSargsList-method}
% Methods LineWise
\alias{codeLine}
\alias{codeLine,SYSargsList-method}
\alias{copyEnvir}
\alias{copyEnvir,SYSargsList-method}
\alias{viewEnvir}
\alias{viewEnvir,SYSargsList-method}
\alias{cmdlist}
\alias{cmdlist,SYSargsList-method}
% Replacement Methods
\alias{[[<-,SYSargsList,ANY,ANY,ANY-method}
\alias{[[<-,SYSargsList,ANY,ANY-method}
\alias{appendStep<-}
\alias{appendStep<-,SYSargsList-method}
\alias{yamlinput<-}
\alias{yamlinput<-,SYSargsList-method}
\alias{replaceStep<-}
\alias{replaceStep<-,SYSargsList-method}
\alias{renameStep<-}
\alias{renameStep<-,SYSargsList-method}
\alias{dependency<-}
\alias{dependency<-,SYSargsList-method}
\alias{appendCodeLine<-}
\alias{appendCodeLine<-,SYSargsList-method}
\alias{replaceCodeLine<-}
\alias{replaceCodeLine<-,SYSargsList-method}
\alias{updateColumn<-}
\alias{updateColumn<-,SYSargsList-method}

\title{Class \code{"SYSargsList"}}
\description{
\code{SYSargsList} S4 class is a list-like container where each instance stores 
all the input/output paths and parameter components required for a particular 
data analysis step based on command-line- or R-based software.

\code{SYSargsList} instances are constructed by the \code{SYSargsList} function.
}
\usage{
## Accessors

stepsWF(x)
statusWF(x)
targetsWF(x)
outfiles(x)
SE(x)
dependency(x)
projectInfo(x)
runInfo(x)

## Methods

cmdlist(x, ...)
codeLine(x, ...)
SampleName(x, ...)
stepName(x)
baseCommand(x, ...)
targetsheader(x, ...)
yamlinput(x, ...)
viewEnvir(x, silent = FALSE)
copyEnvir(x, list = character(), new.env = globalenv(), silent = FALSE)
addResources(x, step, resources)

## Subset Methods

subset(x, ...)
getColumn(x, step, position = c("outfiles", "targetsWF"), column = 1, 
          names = SampleName(x, step))

## Replacement

appendStep(x, after = length(x), ...) <- value
yamlinput(x, paramName, ...) <- value
replaceStep(x, step, step_name = "default") <- value
renameStep(x, step, ...) <- value
dependency(x, step, ...) <- value
appendCodeLine(x, after = length(x), ...) <- value
replaceCodeLine(x, line, ...) <- value
updateColumn(x, step, position = c("outfiles", "targetsWF")) <- value
}

\arguments{
  \item{x}{
  An instance of class \code{SYSargsList}.
}
  \item{step}{
  character or numeric. Workflow step name or position index.
}
  \item{silent}{
  If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{list}{
  a character vector naming objects to be copyied from the enviroment.
}
  \item{new.env}{
  An environment to copy to. Default is \code{globalenv()}.
}
  \item{resources}{
  \code{List} for reserving for each cluster job sufficient computing resources
  including memory (Megabyte), number of nodes (Int), CPU cores, walltime (Minutes),
  etc. It is necessary two additional files: \code{conffile} and \code{template}. 
  \code{conffile} is the path to configuration file (default location 
  \code{./.batchtools.conf.R}). This file contains in its simplest form just one 
  command, such as this line for the scheduler: 
  \code{cluster.functions <- makeClusterFunctionsSlurm(template="batchtools.slurm.tmpl")}. 
  For more detailed information visit this page: https://mllg.github.io/batchtools/index.html.
  \code{template} The template files for a specific queueing/scheduling systems 
  can be downloaded from here: 
  https://github.com/mllg/batchtools/tree/master/inst/templates. Slurm, 
  PBS/Torque, and Sun Grid Engine (SGE) templates are provided within the package
  demo data.
}
  \item{position}{
  character. Options are \code{"outfiles" or "targetsWF"} slots.
}
  \item{column}{
  character or numeric. Which column will be subset from the \code{position}
  argument.
}
  \item{names}{
  character vector. Names of the workflow step. 
}
  \item{after}{
  A subscript, after which the values are to be appended.
}
  \item{paramName}{
  character. Input name from \code{yamlinput(x)}, which value should be replaced. 
}
  \item{step_name}{
  character with the new step name. Default value will automatically give a name:
  \code{step_<step index position>}.
}
  \item{line}{
  numeric. Index position of the code line to be added or replaced.
}
  \item{value}{
  object containing the values to be replaced to \code{SYSargsList}. Values may
  be of the same class as the orginal values. For \code{updateColumn<-} a 
  \code{DataFrame} must have the same rows as the original \code{DataFrame}.
  However, if there is no column/rows in the original \code{DataFrame}, the new
  \code{DataFrame} will replace the empty one. If there is a non-empty 
  \code{DataFrame}, any existing columns with the same name as this new 
  \code{DataFrame} will be replaced. Any columns that do not exist in the 
  original \code{DataFrame} will be added to the original \code{DataFrame}.
}
  \item{...}{
  Further arguments to be passed to or from other methods.
}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SYSargsList", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{stepsWF}:}{Object of class \code{"list"} storing all the 
    steps objects of the workflow. }
    \item{\code{statusWF}:}{Object of class \code{"list"} storing all the 
    success and failure of each step in the workflow. }
    \item{\code{targetsWF}:}{Object of class \code{"list"} storing all the 
    targets \code{DataFrame} for each step in the workflow. For the 
    \code{LineWise} steps, a \code{DataFrame} with 0 rows and 0 columns will
    be displayed. }
    \item{\code{outfiles}:}{Object of class \code{"list"} storing all the 
    output \code{DataFrame} for each step in the workflow. For the 
    \code{LineWise} steps, a \code{DataFrame} with 0 rows and 0 columns will be
    displayed. }
    \item{\code{SE}:}{Object of class \code{"list"} storing all the 
    \code{SummarizedExperiment} objects in the workflow. }
    \item{\code{dependency}:}{Object of class \code{"list"} storing all the 
    dependency graphs in the workflow. }
    \item{\code{projectInfo}:}{Object of class \code{"list"} storing all the 
    \code{projectInfo} information of the workflow. }
    \item{\code{runInfo}:}{Object of class \code{"list"} storing all the
    \code{runInfo} information of each step in the workflow. }
    \item{\code{targets_connection}:}{Object of class \code{"list"} storing all 
    targets files connection in the workflow. }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY", drop = "ANY")}: 
              subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY")}:
              subsetting of class with bracket operator }
    \item{$}{\code{signature(x = "SYSargsList")}: 
              extracting slots elements by name }
    \item{coerce}{\code{signature(from = "list", to = "SYSargsList")}:
              \code{as(list, "SYSargsList")} }
    \item{coerce}{\code{signature(from = "SYSargsList", to = "list")}: 
              \code{as(SYSargsList, "list")} }
    \item{sysargslist}{\code{signature(x = "SYSargsList")}: Coerce back to 
              list \code{as(SYSargsList, "list")} }
    \item{length}{\code{signature(x = "SYSargsList")}: extracts number of 
              \code{SYSargsList} steps }
    \item{names}{\code{signature(x = "SYSargsList")}: extracts slot names }
    \item{show}{\code{signature(object = "SYSargsList")}: summary view of 
              \code{SYSargsList} steps }
    \item{stepsWF}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{stepsWF} slot }
    \item{statusWF}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{statusWF} slot }
    \item{targetsWF}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{targetsWF} slot }
    \item{outfiles}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{outfiles} slot }
    \item{SE}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{SE} slot }
    \item{dependency}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{dependency} slot }
    \item{projectInfo}{\code{signature(x = "SYSargsList")}: extract data from
              \code{projectInfo} slot }
    \item{runInfo}{\code{signature(x = "SYSargsList")}: extract data from
              \code{runInfo} slot }
    \item{cmdlist}{\code{signature(x = "SYSargsList", ...)}: 
              extracts data from \code{cmdlist} slot for each \code{SYSargs2} 
              step }
    \item{codeLine}{\code{signature(x = "SYSargsList", step)}: 
              extracts data from \code{codeLine} slot for \code{LineWise} 
              step }
    \item{SampleName}{\code{signature(x = "SYSargsList", step)}: 
              extracts Sample ID from \code{SYSargs2} instance step }
    \item{stepName}{\code{signature(x = "SYSargsList")}: 
              extracts steps names from workflow instance }
    \item{baseCommand}{\code{signature(x = "SYSargsList", step)}: 
              extracts \code{baseCommand} from \code{SYSargs2} instance step }
    \item{targetsheader}{\code{signature(x = "SYSargsList", step)}: 
               extracts \code{targetsheader} from \code{SYSargs2} instance step }
    \item{yamlinput}{\code{signature(x = "SYSargsList", step)}: 
              extracts data from \code{yamlinput} slot for each \code{SYSargs2} 
              step }
    \item{viewEnvir}{\code{signature(x = "SYSargsList", silent = FALSE)}: 
              return a vector of character strings giving the names of the 
              objects in the \code{SYSargsList} environment }
    \item{copyEnvir}{\code{signature(x = "SYSargsList", list = character(), 
                            new.env = globalenv(), silent = FALSE)}:
            copy of the contents or select objects from \code{SYSargsList}
            environment and place them into new.env }
    \item{addResources}{\code{signature(x = "SYSargsList", step, resources)}:
            Adds the computing resources for one or multiple steps in the workflow.
            If the particular step(s) is set to be executed "management section,"
            when the resources is added, the step(s) will be executed on the
            "compute section." }
    \item{[[<-}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY", value = "ANY")}: 
              replacement method for \code{SYSargsList} class }
    \item{appendStep<-}{\code{signature(x = "SYSargsList", after = length(x))}:
              insert the \code{SYSargsList} or \code{LineWise} object onto x at 
              the position given by \code{after} }
    \item{yamlinput<-}{\code{signature(x = "SYSargsList", step, paramName )}: 
              replace a value in the \code{yamlinput} slot for a specific step instance }
    \item{replaceStep<-}{\code{signature(x = "SYSargsList", step, step_name = "default")}:
              replace a specific step in the workflow instance }
    \item{renameStep<-}{\code{signature(x = "SYSargsList")}: 
              rename a \code{stepName} in the workflow instance }
    \item{dependency<-}{\code{signature(x = "SYSargsList", step)}: 
              replace dependency graph for a specific step instance }
    \item{appendCodeLine<-}{\code{signature(x = "SYSargsList", step, after = length(x))}:
              insert the R code in a specific step at the position given by \code{after} }
    \item{replaceCodeLine<-}{\code{signature(x = "SYSargsList", step, line)}:
              replace the R code in a specific step at the position given by \code{line} }
    \item{updateColumn<-}{\code{signature(x = "SYSargsList", step, 
              position = c("outfiles", "targetsWF"))}:
              update or add a new column in targetsWF or outfiles slots }
  }
}
\author{
Daniela Cassol and Thomas Girke
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

length(sal)
names(sal)

## Accessors

stepsWF(sal)
statusWF(sal)
targetsWF(sal)
outfiles(sal)
SE(sal)
dependency(sal)
projectInfo(sal)
runInfo(sal)

## Methods

cmdlist(sal, step=1, targets=1:2) ## SYSargs2 step
codeLine(sal, step=2) ## LineWise step
SampleName(sal, step="echo")
stepName(sal)
baseCommand(sal, 1) ## SYSargs2 step
targetsheader(sal, step=1) ## SYSargs2 step
yamlinput(sal, step=1) ## SYSargs2 step
viewEnvir(sal)
copyEnvir(sal, list = character(), new.env = globalenv())

resources <- list(conffile= system.file("extdata/.batchtools.conf.R", 
                            package="systemPipeR"),
                  template= system.file("extdata/batchtools.slurm.tmpl", 
                            package="systemPipeR"),
                  Njobs=3, ## Usually, the samples number
                  walltime=60, ## minutes
                  ntasks=1,
                  ncpus=4, 
                  memory=1024 ## Mb
)
addResources(sal, 1, resources= resources)

## Subset Methods

sal_sub <- subset(sal, subset_steps=1, input_targets=1:2, keep_steps = TRUE)
sal_sub
targetsIn <- getColumn(sal, step=1, position = c("outfiles"))
targetsIn

## Replacement
renameStep(sal, step=1) <- "new_echo"
dependency(sal, step=2) <- "new_echo"
updateColumn(sal, step=2, position = c("targetsWF")) <- data.frame(targetsIn)
targetsWF(sal)

replaceStep(sal, step=2) <- LineWise(code = {
                            hello <- "Printing a new message"
                            }, 
                            step_name = "R_hello", 
                            dependency = "new_echo")
codeLine(sal)
yamlinput(sal, step=1, paramName="results_path") <- list(results_path=list(
                                            class="Directory", path="./data"))
cmdlist(sal, step = 1, targets = 1)
appendCodeLine(sal, step=2, after = 0) <- "log <- log(10)"
codeLine(sal, 2)
replaceCodeLine(sal, step=2, line=1) <- LineWise(code = {
                            log <- log(50)
                            })
codeLine(sal, 2)
}
\keyword{ classes }

