\name{writeTargetsRef}
\alias{writeTargetsRef}
\title{
Generate targets file with reference
}
\description{
Generates targets file with sample-wise reference as required for some NGS
applications, such as ChIP-Seq containing input samples. The reference sample
information needs to be provided in the input file in a column called
\code{SampleReference} where the values reference the labels used in the
\code{SampleName} column. Sample rows without reference assignments will be
removed automatically.
}
\usage{
writeTargetsRef(infile, outfile, silent = FALSE, overwrite = FALSE, ...)
}
\arguments{
  \item{infile}{
Path to input targets file.
}
  \item{outfile}{
Path to output targets file.
}
  \item{silent}{
If set to \code{TRUE}, all messages returned by the function will be suppressed.
}
  \item{overwrite}{
If set to \code{TRUE}, existing files of same name will be overwritten.
}
  \item{\dots}{
To pass on additional arguments.
}
}
\value{
Generates modified targets file with the paths to the reference samples in the second column named \code{FileName2}. 
Note, sample rows not assigned reference samples are removed automatically.
}
\author{
Thomas Girke
}
\seealso{
\code{writeTargetsout}
}
\examples{
## Path to input targets file
targets <- system.file("extdata", "targets_chip.txt", package="systemPipeR")

\dontrun{
## Write modified targets file with reference (e.g. input) sample
writeTargetsRef(infile=targets, outfile="~/targets_refsample.txt", silent=FALSE, overwrite=FALSE)
}
}
\keyword{ utilities }
