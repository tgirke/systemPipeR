\name{INTERSECTset-class}
\Rdversion{1.1}
\docType{class}
\alias{INTERSECTset-class}
\alias{as.list,INTERSECTset-method}
\alias{coerce,list,INTERSECTset-method}
\alias{complexitylevels,INTERSECTset-method}
\alias{intersectlist,INTERSECTset-method}
\alias{intersectmatrix,INTERSECTset-method}
\alias{length,INTERSECTset-method}
\alias{names,INTERSECTset-method}
\alias{setlist,INTERSECTset-method}
\alias{show,INTERSECTset-method}

\title{Class \code{"INTERSECTset"}}
\description{
Container for storing standard intersect results created by the \code{overLapper} function.
The \code{setlist} slot stores the original label sets as \code{vectors} in a \code{list};
\code{intersectmatrix} organizes the label sets in a present-absent matrix; \code{complexitylevels}
represents the number of comparisons considered for each comparison set as vector of integers;
and \code{intersectlist} contains the standard intersect vectors. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("INTERSECTset", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{setlist}:}{Object of class \code{"list"}: \code{list} of \code{vectors} }
    \item{\code{intersectmatrix}:}{Object of class \code{"matrix"}: binary \code{matrix} }
    \item{\code{complexitylevels}:}{Object of class \code{"integer"}: \code{vector} of \code{integers} }
    \item{\code{intersectlist}:}{Object of class \code{"list"}: \code{list} of \code{vectors} }
  }
}
\section{Methods}{
  \describe{
    \item{as.list}{\code{signature(x = "INTERSECTset")}: coerces \code{INTERSECTset} to \code{list} }
    \item{coerce}{\code{signature(from = "list", to = "INTERSECTset")}: \code{as(list, "INTERSECTset")} }
    \item{complexitylevels}{\code{signature(x = "INTERSECTset")}: extracts data from \code{complexitylevels} slot }
    \item{intersectlist}{\code{signature(x = "INTERSECTset")}: extracts data from \code{intersectlist} slot }
    \item{intersectmatrix}{\code{signature(x = "INTERSECTset")}: extracts data from \code{intersectmatrix} slot }
    \item{length}{\code{signature(x = "INTERSECTset")}: returns number of original label sets }
    \item{names}{\code{signature(x = "INTERSECTset")}: extracts slot names }
    \item{setlist}{\code{signature(x = "INTERSECTset")}: extracts data from \code{setlist} slot }
    \item{show}{\code{signature(object = "INTERSECTset")}: summary view of \code{INTERSECTset} objects }
	 }
}
\author{
Thomas Girke
}
\seealso{
\code{overLapper}, \code{vennPlot}, \code{olBarplot}, \code{VENNset-class}
}
\examples{
showClass("INTERSECTset")

## Sample data
setlist <- list(A=sample(letters, 18), B=sample(letters, 16),
                C=sample(letters, 20), D=sample(letters, 22),
                E=sample(letters, 18), F=sample(letters, 22))

## Create VENNset
interset <- overLapper(setlist[1:5], type="intersects")
class(interset)

## Accessor methods for VENNset/INTERSECTset objects
names(interset)
setlist(interset)
intersectmatrix(interset)
complexitylevels(interset)
intersectlist(interset)

## Coerce VENNset/INTERSECTset object to list
as.list(interset)
}
\keyword{classes}
