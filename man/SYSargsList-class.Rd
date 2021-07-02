\name{SYSargsList-class}
\Rdversion{1.1}
\docType{class}
\alias{SYSargsList-class}
\alias{[,SYSargsList,ANY,ANY,ANY-method}
\alias{[[,SYSargsList,ANY,missing-method}
\alias{[[,SYSargsList,ANY,ANY-method}
\alias{[[<-,SYSargsList,ANY,ANY,ANY-method}
\alias{$,SYSargsList-method}
\alias{codeSteps,SYSargsList-method}
\alias{coerce,list,SYSargsList-method}
\alias{coerce,SYSargsList,list-method}
\alias{length,SYSargsList-method}
\alias{names,SYSargsList-method}
\alias{show,SYSargsList-method}
\alias{sysargslist}
\alias{sysargslist,SYSargsList-method}
\alias{stepsWF}
\alias{stepsWF,SYSargsList-method}
\alias{statusWF}
\alias{statusWF,SYSargsList-method}
\alias{dependency}
\alias{dependency,SYSargsList-method}
\alias{projectInfo}
\alias{projectInfo,SYSargsList-method}
\alias{runInfo}
\alias{runInfo,SYSargsList-method}
\alias{targetsWF}
\alias{targetsWF,SYSargsList-method}
\alias{SEobj}
\alias{SEobj,SYSargsList-method}
\alias{outfiles}
\alias{outfiles,SYSargsList-method}
\alias{targets_connection}
\alias{targets_connection,SYSargsList-method}
\alias{subsetTargets}
\alias{subsetTargets,SYSargsList-method}
\alias{getColumn}
\alias{getColumn,SYSargsList-method}
\alias{copyEnvir}
\alias{copyEnvir,SYSargsList-method}
\alias{viewEnvir}
\alias{viewEnvir,SYSargsList-method}
\alias{cmdlist,SYSargsList-method}

\alias{SampleName,SYSargsList-method}
\alias{stepName,SYSargsList-method}
\alias{targetsheader,SYSargsList-method}

\alias{appendStep<-}
\alias{appendStep<-,SYSargsList-method}
\alias{replaceStep<-}
\alias{replaceStep<-,SYSargsList-method}
\alias{renameStep<-}
\alias{renameStep<-,SYSargsList-method}
\alias{yamlinput,SYSargsList-method}
\alias{yamlinput<-,SYSargsList-method}
\alias{yamlinput<-,SYSargsList,ANY,ANY-method}
\alias{appendCodeLine<-,SYSargsList-method}
\alias{replaceCodeLine<-,SYSargsList-method}
\alias{codeLine,SYSargsList-method}

\title{Class \code{"SYSargsList"}}
\description{
S4 class container for storing parameters and information of command-line- or R-based workflow.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SYSargsList", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{stepsWF}:}{Object of class \code{"list"} storing all the \code{SYSargs2} objects of the workflow }
    \item{\code{statusWF}:}{Object of class \code{"list"} storing all the success and failure of each step in the workflow }
    \item{\code{targetsWF}:}{Object of class \code{"list"} storing all the targets files for each step in the workflow}
    \item{\code{outfiles}:}{Object of class \code{"list"} storing all the output files in the workflow }
    \item{\code{SEobj}:}{Object of class \code{"list"} storing all the SummarizedExperiment objects in the workflow }
    \item{\code{dependency}:}{Object of class \code{"list"} storing all the dependency steps in the workflow }
    \item{\code{targets_connection}:}{Object of class \code{"list"} storing all targets files connection in the workflow }
    \item{\code{projectInfo}:}{Object of class \code{"list"} storing all the \code{projectInfo} information of the workflow }
    \item{\code{runInfo}:}{Object of class \code{"list"} storing all the \code{runInfo} information of the workflow }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY", drop = "ANY")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY")}: subsetting of class with bracket operator }
    \item{[[<-}{\code{signature(x = "SYSargsList", i = "ANY", j = "ANY", value = "ANY")}: replacement method for \code{SYSargsList} class }
    \item{$}{\code{signature(x = "SYSargsList")}: extracting slots elements by name }
    \item{coerce}{\code{signature(from = "list", to = "SYSargsList")}: \code{as(list, "SYSargsList")} }
    \item{coerce}{\code{signature(from = "SYSargsList", to = "list")}: \code{as(SYSargsList, "list")} }
    \item{length}{\code{signature(x = "SYSargsList")}: extracts number of \code{SYSargs2} objects }
    \item{names}{\code{signature(x = "SYSargsList")}: extracts slot names }
    \item{show}{\code{signature(object = "SYSargsList")}: summary view of \code{SYSargsList} steps }
    \item{stepsWF}{\code{signature(x = "SYSargsList")}: extract data from \code{stepsWF} slot }
    \item{statusWF}{\code{signature(x = "SYSargsList")}: extract data from \code{statusWF} slot }
    \item{dependency}{\code{signature(x = "SYSargsList")}: extract data from \code{dependency} slot }
    \item{projectInfo}{\code{signature(x = "SYSargsList")}: extract data from \code{projectInfo} slot }
    \item{runInfo}{\code{signature(x = "SYSargsList")}: extract data from \code{runInfo} slot }
    \item{targetsWF}{\code{signature(x = "SYSargsList")}: extract data from \code{targetsWF} slot }
    \item{sysargslist}{\code{signature(x = "SYSargsList")}: Coerce back to list \code{as(SYSargsList, "list")} }
    \item{SEobj}{\code{signature(x = "SYSargsList")}: extract data from \code{SEobj} slot }
    \item{cmdlist}{\code{signature(x = "SYSargsList", sample=NULL)}: extracts data from \code{cmdlist} slot for each step }
    \item{appendStep<-}{\code{signature(x = "SYSargsList")}: replacement method for append a step }
    \item{appendStep<-}{\code{signature(x = "SYSargsList")}: replacement method for append a step }
    \item{replaceStep<-}{\code{signature(x = "SYSargsList")}: replacement method for replace a step }
    \item{renameStep<-}{\code{signature(x = "SYSargsList")}: replacement method for rename a step}
    \item{yamlinput<-}{\code{signature(x = "SYSargs2")}: replacement method for \code{yamlinput} slot for each step }
  }
}
\author{
Daniela Cassol and Thomas Girke
}
\seealso{
\code{SYSargs2} and \code{initWF}
}
\examples{
showClass("SYSargsList")
sal <- SPRproject(projPath = tempdir())
targetspath <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
appendStep(sal) <- SYSargsList(targets=targetspath, 
                   step_name="Quality",
                   wf_file = "trimmomatic/trimmomatic-pe.cwl", input_file = "trimmomatic/trimmomatic-pe.yml",
                   dir_path=system.file("extdata/cwl", 
                                        package = "systemPipeR"),
                   inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
}
\keyword{ classes }

