\name{genFeatures}
\alias{genFeatures}
\title{
Generate feature ranges from TxDb
}
\description{
Function to generate various feature types from \code{TxDb} objects including
\code{intergenic}, \code{promoter}, \code{intron}, \code{exon}, \code{cds},
\code{5'/3'UTR} and different \code{transcript} types. The latter contains 
as many transcript types as available in the \code{tx_type} column when extracting
transcripts from \code{TxDb} objects as follows: 
\code{transcripts(txdb, c("tx_name", "gene_id", "tx_type"))}
}
\usage{
genFeatures(txdb, featuretype = "all", upstream = 1000, downstream = 0, verbose = TRUE)
}
\arguments{
  \item{txdb}{
\code{TxDb} object
}
  \item{featuretype}{
Feature types can be selected by assigning a \code{character} vector containing 
any of the following: c("tx_type", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic")
The default \code{all} is a shorthand to select all supported features.
}
  \item{upstream}{
Defines for promoter features the number of bases upstream from the transcription
start site. 
}
  \item{downstream}{
Defines for promoter features the number of bases downstream from the transcription
start site.
}
  \item{verbose}{
\code{verbose=FALSE} turns off all print messages.
}
}
\value{
The results are returned as a \code{GRangesList} where each component is a
\code{GRanges} object containing the range set of each feature type. Note,
overlapping feature components of the same type and from the same gene are
merged/reduced, e.g. two overlapping exons from the same gene are merged to
one.
}
\author{
Thomas Girke
}
\seealso{
\code{transcripts} and associated \code{TxDb} accessor functions from 
\code{GenomicFeatures} package.
}
\examples{
## Sample from GenomicFeatures package
library(GenomicFeatures)
gffFile <- system.file("extdata", "GFF3_files", "a.gff3", package="GenomicFeatures")
txdb <- makeTxDbFromGFF(file=gffFile, format="gff3", organism="Solanum lycopersicum")
feat <- genFeatures(txdb, featuretype="all", upstream=1000, downstream=0, verbose=TRUE)

## List extracted feature types
names(feat)

## Return all features in single GRanges object
unlist(feat)

\dontrun{
## Sample from systemPipeRdata package
file <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
txdb <- makeTxDbFromGFF(file=file, format="gff3", organism="Arabidopsis")
feat <- genFeatures(txdb, featuretype="all", upstream=1000, downstream=0, verbose=TRUE)
}
}
\keyword{ utilities }
