\name{countRangeset}
\alias{countRangeset}
\title{
Read counting for several range sets
}
\description{
Convenience function to perform read counting over serveral different range
sets, e.g. peak range sets or feature types. Internally, the read counting
is performed with the \code{summarizeOverlaps} function from the \code{GenomicAlignments}
package. The resulting count tables are directly saved to files.
}
\usage{
countRangeset(bfl, rangefiles, extension = "_countDF.xls", ...)
}
\arguments{
  \item{bfl}{
\code{BamFileList} object containing paths to one or more BAM files.
}
  \item{rangefiles}{
Character vector containing paths to tabular range data files (e.g. peak
ranges). The files need to contain proper columns titles to coerce with
\code{as(..., "GRanges")} to \code{GRanges} objects after importing them
with \code{read.delim}.
}
  \item{extension}{
Chracter string specifying the file extensions of the result files.
}
  \item{\dots}{
Arguments to be passed on to internally used \code{summarizeOverlaps} function.
}
}
\value{
Named character vector containing the paths to the resulting count table files. 
The base name of the count files is the same as the corresponding input files 
specified under \code{rangefiles} and the value of \code{extension} appended. 
}
\author{
Thomas Girke
}
\seealso{
\code{summarizeOverlaps}
}
\examples{
## Paths to BAM files
param <- system.file("extdata", "bowtieSE.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args_bam <- systemArgs(sysma=param, mytargets=targets)
bfl <- BamFileList(outpaths(args_bam), yieldSize=50000, index=character())

\dontrun{
## Paths to range files
args_peak <- systemArgs(sysma="param/macs2_noinput.param", mytargets="targets_mergeBamByFactor.txt")
rangefiles <- outpaths(args_peak)

## Iterative read counting
countDFnames <- countRangeset(bfl, rangefiles, extension="_countDF.xls", mode="Union", ignore.strand=TRUE)
}
}
\keyword{ utilities }
