\name{cwlFilesUpdate}
\alias{cwlFilesUpdate}
\title{Update CWL description files}
\description{
Function to sync and download the most updated CWL description files from
the systemPipeR package.
}
\usage{
cwlFilesUpdate(destdir, force = FALSE, verbose = TRUE)
}
\arguments{
\item{destdir}{a character string with the directory path where the param are
stored.}

\item{force}{logical. Force the download of the CWL files.}

\item{verbose}{logical. The setting \code{FALSE} suppresses all print messages.}
}
\author{
Daniela Cassol
}

\examples{
\dontrun{
destdir <- "param/"
cwlFilesUpdate(destdir)
}
}
\keyword{ utilities }
