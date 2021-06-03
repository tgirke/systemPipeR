\name{LineWise-class}
\Rdversion{1.1}
\docType{class}
% Class
\alias{LineWise-class}
% Subsetting and Quick access
\alias{[,LineWise,ANY,ANY,ANY-method}
\alias{[[,LineWise,ANY,missing-method}
\alias{[[,LineWise,ANY,ANY-method}
\alias{[[<-,LineWise,ANY,ANY,ANY-method}
\alias{$,LineWise-method}
% Coercion
\alias{coerce,LineWise,list-method}
\alias{coerce,list,LineWise-method}
%Accessors
\alias{length,LineWise-method}
\alias{names,LineWise-method}
\alias{show,LineWise-method}
\alias{codeLine}
\alias{codeLine,LineWise-method}
\alias{rmdStart}
\alias{rmdStart,LineWise-method}
\alias{linewise}
\alias{linewise,LineWise-method}
\alias{script}
\alias{script,LineWise-method}
%Replacement
\alias{appendCodeLine,LineWise-method}
\alias{appendCodeLine<-}
\alias{appendCodeLine<-,LineWise-method}
\alias{appendCodeLine<-,LineWise,ANY,ANY,ANY-method}
\alias{replaceCodeLine<-}
\alias{replaceCodeLine,LineWise-method}
\alias{replaceCodeLine<-,LineWise-method}
\alias{replaceCodeLine<-,LineWise,ANY,ANY,ANY-method}

\title{Class \code{"LineWise"}}
\description{
\code{LineWise class} stores all the R chunk from a RMarkdown file. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LineWise", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{codeLine}:}{Object of class \code{"character"} storing data from Rc chunk. }
    \item{\code{rmdStart}:}{Object of class \code{"interger"} storing Start line from the script file. }
    \item{\code{script}:}{Object of class \code{"character"} storing path from the script file. }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "LineWise")}: subsetting of class with bracket operator }
    \item{[[}{\code{signature(x = "LineWise", i = "ANY", j = "missing")}: subsetting of class with bracket operator }
    \item{[[<-}{\code{signature(x = "LineWise")}: replacement method for LineWise class }
    \item{$}{\code{signature(x = "LineWise")}: extracting slots elements by name }
    \item{codeLine}{\code{signature(x = "LineWise")}: ... }
    
    \item{coerce}{\code{signature(from = "LineWise", to = "list")}: \code{as(LineWise, "list")} }
    \item{coerce}{\code{signature(from = "list", to = "LineWise")}: \code{as(list, "LineWise")} }
    \item{length}{\code{signature(x = "LineWise")}: extracts number of R code lines }
    \item{linewise}{\code{signature(x = "LineWise")}: Coerce back to list \code{as(LineWise, "list")} }
    \item{names}{\code{signature(x = "LineWise")}: extracts slot names }
    \item{script}{\code{signature(x = "LineWise")}: extracts path of the Rmarkdown files }
    \item{show}{\code{signature(object = "LineWise")}: summary view of \code{LineWise} objects }
    \item{appendCodeLine<-}{\code{signature(x = "LineWise")}: replacement method for append a code line }
    \item{replaceCodeLine<-}{\code{signature(x = "LineWise")}: replacement method for replace a code line }
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
}
\keyword{classes}
