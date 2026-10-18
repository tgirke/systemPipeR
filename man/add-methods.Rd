\name{addAssay-methods}
\alias{addAssay-methods}
\alias{addAssay,SummarizedExperiment-method}
\alias{addMetadata-methods}
\alias{addMetadata,SummarizedExperiment-method}
\alias{addAssay}
\alias{addMetadata}

\title{Extension accessor methods for SummarizedExperiment object}
\description{
Accessors for adding new data to the 'assay' and 'metadata' slot of a SummarizedExperiment object
}
\usage{
addAssay(x, ...)
addMetadata(x, ...)
}
\arguments{
  \item{x}{
	Object of class \code{SummarizedExperiment}. 
}
  \item{...}{
	dots, name of the object. 
}
}
\section{Methods}{
\describe{
   \item{addAssay}{\code{signature(x = "SummarizedExperiment")}: add new dataset to \code{assays} slot }
    \item{addMetadata}{\code{signature(x = "SummarizedExperiment")}: add new dataset to \code{metadata} slot }
}
}
\value{
The \code{SummarizedExperiment} object \code{x}, with \code{value}
appended as a new assay (\code{addAssay}) or as a new \code{metadata}
element (\code{addMetadata}), named \code{xName}.
}
\author{
Daniela Cassol
}
\keyword{methods}

