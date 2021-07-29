\name{showDF}
\alias{showDF}
\title{Create an HTML table using DT package with fixed columns}
\usage{
showDF(data, ...)
}
\arguments{
\item{data}{data object (either a matrix or a data frame).}

\item{...}{Additional arguments used by dDT::atatable() function.}
}
\value{
returns an object of \code{datatables} and \code{htmlwidget}.
}
\description{
Create an HTML table using DT package with fixed columns
}
\examples{
showDF(iris)
}
