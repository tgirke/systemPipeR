\name{LineWise-class}
\Rdversion{1.1}
\docType{class}
% Class
\alias{LineWise-class}
% Constructor 
\alias{LineWise}
% Subsetting and Quick access
\alias{[,LineWise,ANY,ANY,ANY-method}
\alias{[[,LineWise,ANY,missing-method}
%\alias{[[,LineWise,ANY,ANY-method}
\alias{[[<-,LineWise,ANY,ANY,ANY-method}
\alias{$,LineWise-method}
\alias{length,LineWise-method}
\alias{names,LineWise-method}
\alias{show,LineWise-method}
% Coercion
\alias{coerce,LineWise,list-method}
\alias{coerce,list,LineWise-method}
\alias{linewise}
\alias{linewise,LineWise-method}
%Accessors
\alias{codeLine}
\alias{codeLine,LineWise-method}
\alias{codeChunkStart}
\alias{codeChunkStart,LineWise-method}
\alias{rmdPath}
\alias{rmdPath,LineWise-method}
\alias{stepName}
\alias{stepName,LineWise-method}
\alias{dependency,LineWise-method}
\alias{status,LineWise-method}
\alias{files,LineWise-method}
%Replacement
\alias{appendCodeLine,LineWise-method}
\alias{appendCodeLine<-}
\alias{appendCodeLine<-,LineWise-method}
%\alias{appendCodeLine<-,LineWise,ANY,ANY,ANY-method}
\alias{replaceCodeLine<-}
\alias{replaceCodeLine,LineWise-method}
\alias{replaceCodeLine<-,LineWise-method}
%\alias{replaceCodeLine<-,LineWise,ANY,ANY,ANY-method}

\title{Class \code{"LineWise"}}
\description{
S4 class container for storing R-based code for a workflow step.
\code{LineWise class} instances are constructed by the \code{LineWise} function, 
based on the R-based code, step name, and dependency tree. 
When the container is built from the R Markdown, using \code{\link{importWF}}
function, two other slots are populated: \code{codeChunkStart} and \code{rmdPath}.
\code{codeChunkStart} will store the first line of each R chunk, and \code{rmdPath}
will store the R Markdown file path. 
}
\usage{
## Constructor

LineWise(code, step_name="default", codeChunkStart=integer(), rmdPath=character(), 
        dependency=NA)

## Accessors

codeLine(x)
codeChunkStart(x)
rmdPath(x)
stepName(x)

## Replacement

appendCodeLine(x, after = length(x), ...) <- value
replaceCodeLine(x, line, ...) <- value

}
\arguments{
  \item{code}{
  R code separated either by a semi-colon (\code{\;}), or by a newline, and 
  enclosed by braces (\code{\{ \}}).
}
  \item{codeChunkStart}{
  integer. R Markdown code chunk line start. This element will be populated when
  the object is built by \code{\link{importWF}}.
}
  \item{rmdPath}{
  character. Path of R Markdown file used by \code{\link{importWF}}.
}
  \item{step_name}{
  character. Step name needs to be unique and is required when appending this
  step to the workflow.
}
  \item{dependency}{
  logical or character. Dependency tree, required when appending this step to 
  the workflow. Character name of a previous step in the workflow. 
  Default is \code{NA}.
}
  \item{x}{
  An instance of class \code{LineWise}.
}
  \item{after}{
  A subscript, after which the values are to be appended.
}
  \item{line}{
  An instance of class \code{LineWise}.
}
  \item{value}{
  character vector object containing the R code to be appended or replaced to x.
}
  \item{...}{
  Options from \code{\link{SYSargsList}} object.
}
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("LineWise", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{codeLine}:}{Object of class \code{"expression"} storing R-based code.}
    \item{\code{codeChunkStart}:}{Object of class \code{"interger"} storing start 
                                line from the \code{rmdPath} file, when the LineWise is 
                                built from R Markdown.}
    \item{\code{rmdPath}:}{Object of class \code{"character"} storing path from the R Markdown file.}
    \item{\code{stepName}:}{Object of class \code{"character"} storing step name.}
    \item{\code{dependency}:}{Object of class \code{"list"} storing dependency tree.}
    \item{\code{status}:}{Object of class \code{"list"} storing status steps.}
    \item{\code{files}:}{Object of class \code{"list"} storing file path containing
                        stdout and stderr after running the R-based code.}
  }
}
\section{Functions and Methods}{
   See 'Usage' for details on invocation.
  
  Constructor:
  \describe{
    \item{LineWise:}{Returns a \code{LineWise} object.}
  }

  Accessors:  
  \describe{
    % Accessors
    \item{codeLine}{Printing method for the \code{CodeLine} slot.}
    \item{codeChunkStart}{Extract start line of the R Markdown R chunk.}
    \item{rmdPath}{Extract Rmarkdown file path.}
    \item{stepName}{Extract the step name.}
    \item{dependency}{Extract the dependency tree.}
    \item{status}{Extract status of the step.}
    \item{files}{Extract log file path storing stdout and stderr after running step.}
    % Replacement
    \item{appendCodeLine<-}{Replacement method for append a R code line.}
    \item{replaceCodeLine<-}{Replacement method for replace a R code line.}
 }
 
   Methods:  
  \describe{
    % General
    \item{[}{Return a new \code{LineWise} object made of the selected R code lines.}
    \item{[[}{Extract the slot information from \code{LineWine} object.}
    \item{[[<-}{Replacement method for \code{LineWine} slots.}
    \item{$}{Extract slots elements by name.}
    \item{length}{Extract number of R-based code lines.}
    \item{names}{Extract slot names }
    \item{show}{summary view of \code{LineWise} elements.}
    % Coerce
    \item{coerce}{\code{signature(from = "LineWise", to = "list")}\code{as(LineWise, "list")} }
    \item{coerce}{\code{signature(from = "list", to = "LineWise")}\code{as(list, "LineWise")} }
    \item{linewise}{Coerce back to list \code{as(LineWise, "list")}.}
 }
}
\author{
Daniela Cassol
}
\seealso{
 \code{\linkS4class{SYSargsList}}
}
\examples{
showClass("LineWise")
lw <- LineWise(code = {
              log_out <- log(10)
              },
              step_name = "R_log")

codeLine(lw)

## ImportWF option
file_path <- system.file("extdata/spr_simple_lw.Rmd", package="systemPipeR")
sal <- SPRproject(overwrite = TRUE)
# file_path <- "../inst/extdata/spr_simple_lw.Rmd"
sal <- importWF(sal, file_path)
sal <- runWF(sal)
lw2 <- sal$stepsWF[[2]]
lw2
names(lw2)
length(lw2)

## Accessors
codeLine(lw2)
codeChunkStart(lw2)
rmdPath(lw2)
stepName(lw2)
dependency(lw2)
status(lw2)
files(lw2)

## Replacement
appendCodeLine(lw2, after = 0) <- "log <- log(10)"
codeLine(lw2)
replaceCodeLine(lw2, 1) <- "plot(iris)"
codeLine(lw2)

## Coerce
lw2 <- linewise(lw2) ## OR lw2 <- as(lw2, "list") 
lw2 <- as(lw2, "LineWise")
}
\keyword{classes}
