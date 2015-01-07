\name{VENNset-class}
\Rdversion{1.1}
\docType{class}
\alias{VENNset-class}
\alias{as.list,VENNset-method}
\alias{coerce,list,VENNset-method}
\alias{complexitylevels,VENNset-method}
\alias{intersectmatrix,VENNset-method}
\alias{length,VENNset-method}
\alias{names,VENNset-method}
\alias{setlist,VENNset-method}
\alias{show,VENNset-method}
\alias{vennlist,VENNset-method}

\title{Class \code{"VENNset"}}
\description{
Container for storing Venn intersect results created by the \code{overLapper} function.
The \code{setlist} slot stores the original label sets as \code{vectors} in a \code{list};
\code{intersectmatrix} organizes the label sets in a present-absent matrix; \code{complexitylevels}
represents the number of comparisons considered for each comparison set as vector of integers;
and \code{vennlist} contains the Venn intersect vectors. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("VENNset", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{setlist}:}{Object of class \code{"list"}: \code{list} of \code{vectors} }
    \item{\code{intersectmatrix}:}{Object of class \code{"matrix"}: binary \code{matrix}}
    \item{\code{complexitylevels}:}{Object of class \code{"integer"}: \code{vector} of \code{integers}}
    \item{\code{vennlist}:}{Object of class \code{"list"}: \code{list} of \code{vectors} }
  }
}
\section{Methods}{
  \describe{
    \item{as.list}{\code{signature(x = "VENNset")}: coerces \code{VENNset} to \code{list} }
    \item{coerce}{\code{signature(from = "list", to = "VENNset")}: \code{as(list, "VENNset")} }
    \item{complexitylevels}{\code{signature(x = "VENNset")}: extracts data from \code{complexitylevels} slot }
    \item{intersectmatrix}{\code{signature(x = "VENNset")}: extracts data from \code{intersectmatrix} slot }
    \item{length}{\code{signature(x = "VENNset")}: returns number of original label sets }
    \item{names}{\code{signature(x = "VENNset")}: extracts slot names }
    \item{setlist}{\code{signature(x = "VENNset")}: extracts data from \code{setlist} slot }
    \item{show}{\code{signature(object = "VENNset")}: summary view of \code{VENNset} objects }
    \item{vennlist}{\code{signature(x = "VENNset")}: extracts data from \code{vennset} slot }
	 }
}
\author{
Thomas Girke
}
\seealso{
\code{overLapper}, \code{vennPlot}, \code{olBarplot}, \code{INTERSECTset-class}
}
\examples{
showClass("VENNset")

## Sample data
setlist <- list(A=sample(letters, 18), B=sample(letters, 16),
                C=sample(letters, 20), D=sample(letters, 22),
                E=sample(letters, 18), F=sample(letters, 22))

## Create VENNset
vennset <- overLapper(setlist[1:5], type="vennsets")
class(vennset)

## Accessor methods for VENNset/INTERSECTset objects
names(vennset)
setlist(vennset)
intersectmatrix(vennset)
complexitylevels(vennset)
vennlist(vennset)

## Coerce VENNset/INTERSECTset object to list
as.list(vennset)
}
\keyword{classes}
