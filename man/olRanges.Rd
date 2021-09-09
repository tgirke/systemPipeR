\name{olRanges}
\alias{olRanges}
\title{
Identify Range Overlaps for IRanges and GRanges Object
}
\description{
Function for identifying consensus peak among two peaks sets sharing a minimum relative overlap.
}
\usage{
olRanges(query, subject, output = "gr")
}
\arguments{
  \item{query}{
  Object of class \code{GRanges}, which is a vector of genomic locations and associated annotations.
}
  \item{subject}{
  Object of class \code{GRanges}.

}
  \item{output}{
  By default \code{"gr"} returns any overlap with OL length information in an object of class \code{GRanges}. Also, can returns an object of class \code{data.frame} with \code{"df"}.
}
}
\author{
Thomas Girke
}
\examples{
## Sample Data Sets
grq <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), ranges = IRanges(seq(1, 100, by=10), end = seq(30, 120, by=10)), strand = Rle(strand(c("-", "+", "-")), c(1, 7, 2)))
grs <- shift(grq[c(2,5,6)], 5)
## Run olRanges function
olRanges(query=grq, subject=grs, output="df") 
olRanges(query=grq, subject=grs, output="gr") 
}
\keyword{ utilities }