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
\value{
Called for its side effect: copies the \code{cwl} and \code{docopt.R}
directories bundled with the installed package into \code{destdir} when
\code{force = TRUE} or when the bundled version is newer. Returns
\code{NULL} invisibly.
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
