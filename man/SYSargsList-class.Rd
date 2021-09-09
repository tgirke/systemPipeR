\name{SYSargsList}
\Rdversion{1.1}
\docType{class}
\alias{SYSargsList-class}
% Constructor
\alias{SYSargsList}
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
%\alias{[[,SYSargsList,ANY,ANY-method}
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
\alias{dependency<-}
\alias{dependency<-,SYSargsList-method}
\alias{appendStep<-}
\alias{appendStep<-,SYSargsList-method}
\alias{replaceStep<-}
\alias{replaceStep<-,SYSargsList-method}
\alias{renameStep<-}
\alias{renameStep<-,SYSargsList-method}
\alias{yamlinput<-}
\alias{yamlinput<-,SYSargsList-method}
\alias{appendCodeLine<-}
\alias{appendCodeLine<-,SYSargsList-method}
\alias{replaceCodeLine<-}
\alias{replaceCodeLine<-,SYSargsList-method}
\alias{updateColumn<-}
\alias{updateColumn<-,SYSargsList-method}
\alias{updateStatus<-}
\alias{updateStatus<-,SYSargsList-method}

\title{Class \code{"SYSargsList"}}
\description{
S4 class container for storing parameters and information of command-line- or 
R-based workflow.

\code{SYSargsList} instances are constructed by the \code{SYSargsList} function.
}
\usage{
## Constructor

SYSargsList(sysargs = NULL, step_name = "default", 
            targets = NULL, wf_file = NULL, input_file = NULL, dir_path = ".", 
            id = "SampleName",
            inputvars = NULL, rm_targets_col = NULL, 
            dir = TRUE, dependency = "", 
            run_step = "mandatory",
            run_session = "rsession",
            silent = FALSE, 
            projPath = getOption("projPath", getwd()))

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

baseCommand(x, ...)
SampleName(x, ...)
stepName(x)
targetsheader(x, ...)
yamlinput(x, ...)
cmdlist(x, ...)
codeLine(x, ...)
copyEnvir(x, list = character(), new.env = globalenv(), silent = FALSE)
viewEnvir(x)

## Subset Methods

subset(x, ...)
getColumn(x, step, position = c("outfiles", "targetsWF"), column = 1, 
          names = SampleName(x, step))



## Replacement

dependency(x, step, ...) <- value
appendStep(x, after = length(x), ...) <- value
replaceStep(x, step, step_name = "default") <- value
renameStep(x, step, ...) <- value
yamlinput(x, paramName, ...) <- value
appendCodeLine(x, after = length(x), ...) <- value
replaceCodeLine(x, line, ...) <- value
updateColumn(x, step, position = c("outfiles", "targetsWF")) <- value
}

\arguments{
  \item{sysargs}{
  \code{SYSargs2} object.
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
  variables list defined in the \code{input} file that matches the column names 
  defined in the \code{targets} file.
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
  Default is empty string \code{""}.
}
  \item{run_step}{
  character. If the step has "mandatory" or "optional" flag for the execution.
}
  \item{run_session}{
  character. If the step has "rsession" or "cluster" flag for the execution.
}
  \item{silent}{
  If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{projPath}{
  a character vector of a full project path name. Default is the current path. 
}
  \item{x}{
  An instance of class \code{SYSargsList}.
}
  \item{step}{
  character or numeric. Step name or index.
}
  \item{list}{
  a character vector naming objects to be copyied from the enviroment.
}
  \item{new.env}{
  An environment to copy to. Default is \code{globalenv()}.
}
  \item{position}{
  character. Options are \code{"outfiles" or "targetsWF"} slots.
}
  \item{column}{
  character or numeric. Which column will be replaces from the position argument.
}
  \item{names}{
  character vector. Names of the new column. 
}
  \item{after}{
  A subscript, after which the values are to be appended.
}
  \item{paramName}{
  character. Input name from \code{yamlinput(x)}, which value should be replaced. 
}
  \item{line}{
  An instance of class \code{LineWise}.
}
  \item{value}{
  object containing the values to be replaced to x. Values may be of the same 
  class as the orginal values. For \code{updateColumn<-}a dataframe that must 
  have the same rows as the modifying dataframe.
  However, if there is no column in the original dataframe. This new df will 
  replace the empty one. If there is a non-empty dataframe, any existing 
  columns with the same name as this new df will be replaced. Any columns 
  that do not exist in the original column but this new df will be added to 
  the original \code{data.frame}.
}
  \item{...}{
  object containing the values to be replaced to x. Values may be of the same 
  class as the orginal values. 
}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SYSargsList", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{stepsWF}:}{Object of class \code{"list"} storing all the 
    \code{SYSargs2} objects of the workflow }
    \item{\code{statusWF}:}{Object of class \code{"list"} storing all the 
    success and failure of each step in the workflow }
    \item{\code{targetsWF}:}{Object of class \code{"list"} storing all the 
    targets files for each step in the workflow}
    \item{\code{outfiles}:}{Object of class \code{"list"} storing all the 
    output files in the workflow }
    \item{\code{SE}:}{Object of class \code{"list"} storing all the 
    \code{SummarizedExperiment} objects in the workflow }
    \item{\code{dependency}:}{Object of class \code{"list"} storing all the 
    dependency steps in the workflow }
    \item{\code{targets_connection}:}{Object of class \code{"list"} storing all 
    targets files connection in the workflow }
    \item{\code{projectInfo}:}{Object of class \code{"list"} storing all the 
    \code{projectInfo} information of the workflow }
    \item{\code{runInfo}:}{Object of class \code{"list"} storing all the
    \code{runInfo} information of the workflow }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY", drop = "ANY")}: 
              subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY")}:
              subsetting of class with bracket operator }
    \item{[[<-}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY", value = "ANY")}: 
              replacement method for \code{SYSargsList} class }
    \item{$}{\code{signature(x = "SYSargsList")}: 
              extracting slots elements by name }
    \item{coerce}{\code{signature(from = "list", to = "SYSargsList")}:
              \code{as(list, "SYSargsList")} }
    \item{coerce}{\code{signature(from = "SYSargsList", to = "list")}: 
              \code{as(SYSargsList, "list")} }
    \item{length}{\code{signature(x = "SYSargsList")}: extracts number of 
              \code{SYSargs2} objects }
    \item{names}{\code{signature(x = "SYSargsList")}: extracts slot names }
    \item{show}{\code{signature(object = "SYSargsList")}: summary view of 
              \code{SYSargsList} steps }
    \item{stepsWF}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{stepsWF} slot }
    \item{statusWF}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{statusWF} slot }
    \item{dependency}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{dependency} slot }
    \item{projectInfo}{\code{signature(x = "SYSargsList")}: extract data from
              \code{projectInfo} slot }
    \item{runInfo}{\code{signature(x = "SYSargsList")}: extract data from
              \code{runInfo} slot }
    \item{targetsWF}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{targetsWF} slot }
    \item{sysargslist}{\code{signature(x = "SYSargsList")}: Coerce back to 
              list \code{as(SYSargsList, "list")} }
    \item{SE}{\code{signature(x = "SYSargsList")}: extract data from 
              \code{SE} slot }
    \item{cmdlist}{\code{signature(x = "SYSargsList", sample=NULL)}: 
              extracts data from \code{cmdlist} slot for each step }
    \item{appendStep<-}{\code{signature(x = "SYSargsList")}:
              replacement method for append a step }
    \item{appendStep<-}{\code{signature(x = "SYSargsList")}: 
              replacement method for append a step }
    \item{replaceStep<-}{\code{signature(x = "SYSargsList")}:
              replacement method for replace a step }
    \item{renameStep<-}{\code{signature(x = "SYSargsList")}: 
              replacement method for rename a step}
    \item{yamlinput<-}{\code{signature(x = "SYSargsList")}: 
              replacement method for \code{yamlinput} slot for each step }
    \item{dependency<-}{\code{signature(x = "SYSargsList")}: 
              replacement method for \code{dependency} slot for each step }
    \item{updateStatus<-}{\code{signature(x = "SYSargsList")}:
              updateStatus method for update a step status }
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

baseCommand(sal, 1) ## SYSargs2 step
SampleName(sal, step="echo")
stepName(sal)
targetsheader(sal, step=1) ## SYSargs2 step
yamlinput(sal, step=1) ## SYSargs2 step
cmdlist(sal, step=1, targets=1:2) ## SYSargs2 step
codeLine(sal, step=2) ## SYSargs2 step
copyEnvir(sal, list = character(), new.env = globalenv())
viewEnvir(sal)

## Subset Methods

sal_sub <- subset(sal, subset_steps=1, input_targets=1:2, keep_steps = TRUE)
targetsIn <- getColumn(sal, step=1, position = c("outfiles"))

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
#yamlinput(sal, step=1, paramName="results_path") <- list(results_path=list(
#                                            class="Directory", path="./data"))
appendCodeLine(sal, step=2, after = 0) <- "log <- log(10)"
codeLine(sal, 2)
replaceCodeLine(sal, step=2, line=1) <- LineWise(code = {
                            log <- log(50)
                            })
codeLine(sal, 2)
}
\keyword{ classes }

