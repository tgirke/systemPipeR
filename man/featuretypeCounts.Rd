\name{featuretypeCounts}
\alias{featuretypeCounts}
\title{
Read counts for feature categories
}
\description{
Counts how many reads in short read alignment files (BAM format) overlap with
entire annotation categories or feature types (e.g. combined mRNA or itergenic
regions). To minimize memory consumption, the BAM files are processed in a stream 
using utilities from the \code{Rsamtools} and  \code{GenomicAlignment} packages. 
The counts can be reported for each read length separately or as a single value 
for reads of any length. Subsequently, the counting results can be plotted with 
the associated \code{plotfeaturetypeCounts} function.
}
\usage{
featuretypeCounts(bfl, grl, singleEnd = TRUE, readlength = NULL, type = "data.frame")
}
\arguments{
  \item{bfl}{
\code{BamFileList} object containing paths to BAM files stored on disk. Note,
memory-efficient processing is achieved by streaming through BAM files rather
than reading entire files into memory at once. The maximum number of alignments 
to process in each iteration is determined by \code{yieldSize}. For 
details see \code{?BamFileList}. 
}
  \item{grl}{
\code{GRangesList} object containing in each list component the range set of
a feature type. Typically, this object is generated with the \code{genFeatures}
function (see below for an example). 
}
  \item{singleEnd}{
Specifies whether BAM files contain alignments for single-end (SE) or paired-end
read data. \code{TRUE} is for SE and \code{FALSE} for PE data. 
}
  \item{readlength}{
Integer vector specifying the read length values for which to report counts
separately. If \code{readlength=NULL} the length of the reads will be ignored
resulting in a single value for each feature type. 

}
  \item{type}{
Determines whether results are returend as \code{data.frame} (\code{type="data.frame"}) 
or as \code{list} (\code{type="list"}). Each list component contains the counting
results for one BAM file and is named after the correponding sample. The 
\code{data.frame} result contains the sample assignemt information in a separate 
column.  
}
}
\value{
The results are returned as \code{data.frame} or \code{list} of \code{data.frames}.
For details see above under \code{types} argument. The result \code{data.frames} contain
the following columns:
\item{SampleName}{Sample names from \code{BamFileList} object.}
\item{Strand}{Sense or antisense strand of read mappings.}  
\item{Featuretype}{Name of feature type provided by \code{GRangesList} object. Note, the total number
of aligned reads is reported under the special feature type 'N_total_aligned'. This value 
is useful for scaling/normalization purposes in plots, e.g. counts per million reads.}
\item{Featuretypelength}{Total length of each reduced feature type in bases. This value is 
useful to normalize the read counts by genomic length units, e.g. in plots.}
\item{Subsequent columns}{Counts for reads of any length or for invidual read lengths.}
}
\author{
Thomas Girke
}
\seealso{
\code{plotfeaturetypeCounts}, \code{genFeatures}
}
\examples{

## Construct SYSargs object from param and targets files 
param <- system.file("extdata", "tophat.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args <- systemArgs(sysma=param, mytargets=targets)

\dontrun{
## Features from sample data of systemPipeRdata package
library(GenomicFeatures)
file <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
txdb <- makeTxDbFromGFF(file=file, format="gff3", organism="Arabidopsis")
feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE, upstream=1000, downstream=0, verbose=TRUE)

## Generate and plot feature counts for specific read lengths
fc <- featuretypeCounts(bfl=BamFileList(outpaths(args), yieldSize=50000), grl=feat, singleEnd=TRUE, readlength=c(74:76,99:102), type="data.frame")
p <- plotfeaturetypeCounts(x=fc, graphicsfile="featureCounts.pdf", graphicsformat="pdf", scales="fixed", anyreadlength=FALSE)

## Generate and plot feature counts for specific read lengths  
fc2 <- featuretypeCounts(bfl=BamFileList(outpaths(args), yieldSize=50000), grl=feat, singleEnd=TRUE, readlength=NULL, type="data.frame")
p2 <- plotfeaturetypeCounts(x=featureCounts2, graphicsfile="featureCounts2.pdf", graphicsformat="pdf", scales="fixed", anyreadlength=TRUE)
}
}
\keyword{ utilities }
