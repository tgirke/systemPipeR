\name{writeTargetsout}
\alias{writeTargetsout}
\title{
Write updated targets out to file
}
\description{
Convenience write function for generating targets files with updated
\code{FileName} columns containing the output paths to files generated
by input/output processes. These processes can be commandline- or R-based software.
Typically, the paths to the inputs are stored in the targets infile 
(\code{targetsin(args)}) and the outputs are stored in the targets outfile 
(\code{targetsout(args)}). Note: the function cannot overwrite 
any existing files. If a file exists then the user has to explicitly 
remove it first. 
}
\usage{
writeTargetsout(x, file = "default", silent = FALSE, ...)
}
\arguments{
  \item{x}{
Object of class \code{SYSargs}.
}
  \item{file}{
Name and path of output file. If set to "default" then the name of
the output file will have the pattern 'targets_<software>.txt', where
<software> will be what \code{software(x)} returns.
}
  \item{silent}{
If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{\dots}{
To pass on additional arguments.
}
}
\value{
Writes tabular targes files containing the header/comment lines from \code{targetsheader(x)}
and the columns from \code{targetsout(x)}. 
}
\author{
Thomas Girke
}
\examples{
## Create SYSargs object
param <- system.file("extdata", "tophat.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args <- systemArgs(sysma=param, mytargets=targets)

\dontrun{
## Write targets out file
writeTargetsout(x=args, file="default")
}
}
\keyword{ utilities }
