\name{genFeatures}
\alias{genFeatures}
\title{
Generate feature ranges from TxDb
}
\description{
Function to generate a variety of feature types from \code{TxDb} objects using 
utilities provided by the \code{GenomicFeatures} package. The feature types are 
organized per gene and can be returned on that level in their non-reduced or 
reduced form.

Currently, supported features include \code{intergenic}, \code{promoter},
\code{intron}, \code{exon}, \code{cds}, \code{5'/3'UTR} and different
\code{transcript} types. The latter contains as many transcript types as
available in the \code{tx_type} column when extracting transcripts from
\code{TxDb} objects as follows: 
\code{transcripts(txdb, c("tx_name", "gene_id", "tx_type"))}
}
\usage{
genFeatures(txdb, featuretype = "all", reduce_ranges, upstream = 1000, downstream = 0, verbose = TRUE)
}
\arguments{
  \item{txdb}{
\code{TxDb} object
}
  \item{featuretype}{
Feature types can be specified by assigning a \code{character} vector containing 
any of the following: c("tx_type", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic").
The default \code{all} is a shorthand to select all supported features.
}
  \item{reduce_ranges}{
If set to \code{TRUE} the feature ranges will be reduced on the gene level. As a result 
overlapping feature components of the same type and from the same gene will be
merged to a single range, e.g. two overlapping exons from the same gene are merged 
to one. Intergenic ranges are not affected by this setting. Note, all reduced feature
types are labeled with the suffix '_red'.
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
\code{GRanges} object containing the range set of each feature type. Intergenic
ranges are assigned unique identifiers and recorded in the \code{featuretype_id} 
column of the metadata block. For this the ids of their adjacent genes are concatenated with two 
underscores as separator. If the adjacent genes overlap with other genes then 
their identifiers are included in the id string as well and separated by a single
underscore.  
}
\author{
Thomas Girke
}
\seealso{
\code{transcripts} and associated \code{TxDb} accessor functions from 
the \code{GenomicFeatures} package.
}
\examples{
## Sample from txdbmaker package
library(txdbmaker)
gffFile <- system.file("extdata", "GFF3_files", "a.gff3", package="txdbmaker")
txdb <- makeTxDbFromGFF(file=gffFile, format="gff3", organism="Solanum lycopersicum")
feat <- genFeatures(txdb, featuretype="all", reduce_ranges=FALSE, upstream=1000, downstream=0)

## List extracted feature types
names(feat)

## Obtain feature lists by genes, here for promoter
split(feat$promoter, unlist(mcols(feat$promoter)$feature_by))

## Return all features in single GRanges object
unlist(feat)

\dontrun{
## Sample from systemPipeRdata package
file <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
txdb <- makeTxDbFromGFF(file=file, format="gff3", organism="Arabidopsis")
feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE, upstream=1000, downstream=0)
}
}
\keyword{ utilities }
