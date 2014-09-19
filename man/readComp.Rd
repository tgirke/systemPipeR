\name{readComp}
\alias{readComp}
\title{
Import sample comparisons from targets file
}
\description{
Parses sample comparisons specified in \code{<CMP>} line(s) of targets file. 
All possible comparisons can be specified with 'CMPset: ALL'.
}
\usage{
readComp(file, format = "vector", delim = "-")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
	Path to targets file.
}
  \item{format}{
	Object type to return: \code{vector} or  \code{matrix}.
}
  \item{delim}{
	Delimiter to use when sample comparisons are returned as \code{vector}.
}
}
\value{
	\code{list} where each component is named according to the name(s) used in the \code{<CMP>} line(s) of the targets file. The list will contain as many sample comparisons sets (list components) as there are sample comparisons lines in the corresponding targets file.
}
\author{
Thomas Girke
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")
readComp(file=targetspath, format="vector", delim="-")
}
\keyword{ utilities }
