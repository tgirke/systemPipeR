\name{cwlFilesUpdate}
\alias{cwlFilesUpdate}
\title{Update CWL description files}
\description{
Function to sync the CWL description files (and the matching docopt.R
scripts) in \code{destdir} with the copies bundled in the installed
systemPipeR package.
}
\usage{
cwlFilesUpdate(destdir, force = FALSE, verbose = TRUE)
}
\arguments{
\item{destdir}{a character string with the directory path where the param are
stored.}

\item{force}{logical. Force copying the CWL files, even if \code{destdir} is already up to date.}

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
