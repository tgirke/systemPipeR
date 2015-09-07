\name{featureCoverage}
\alias{featureCoverage}
\title{
Genome read coverage by transcript models
}
\description{
Computes read coverage along single and multi component features based on
genomic alignments. The coverage segments of component features are spliced to
continuous ranges, such as exons to transcripts or CDSs to ORFs.  The results
can be obtained with single nucleotide resolution (e.g. around start and stop
codons) or as mean coverage of relative bin sizes, such as 100 bins for each
feature. The latter allows comparisons of coverage trends among transcripts of
variable length. The results can be obtained for single or many features (e.g.
any number of transcritpts) at once.  Visualization of the coverage results is
facilitated by a downstream \code{plotfeatureCoverage} function. 
}
\usage{
featureCoverage(bfl, grl, resizereads = NULL, readlengthrange = NULL, Nbins = 20, 
                method = mean, fixedmatrix, resizefeatures, upstream, downstream, 
                outfile, overwrite = FALSE)

}
\arguments{
  \item{bfl}{
Paths to BAM files provided as \code{BamFileList} object. The name slot of the
BAM files will be used for naming samples in the results.
}
  \item{grl}{
Genomic ranges provided as \code{GRangesList} typically generated form
\code{txdb} instances with operations like: \code{cdsBy(txdb, "tx")} or
\code{exonsBy(txdb, "tx")}. Single component features will be processed the
same way as multi component features. 
}
  \item{resizereads}{
Positive integer defining the length alignments should be resized to prior to the 
coverage calculation. \code{NULL} will omit the resizing step. 
}
  \item{readlengthrange}{
Positive integer of length 2 determining the read length range to use for
the coverage calculation. Reads falling outside of the specified length range
will be excluded from the coverage calculation. For instance,
\code{readlengthrange=c(30:40)} will base the coverage calculation on reads
between 30 to 40 bps. Assigning \code{NULL} will skip this filtering step.
}
  \item{Nbins}{
Single positive integer defining the number of segments the coverage of each
feature should be binned into in order to obtain coverage summaries of constant
length, e.g. for plotting purposes. 
}
  \item{method}{
Defines the summary statistics to use for binning. The default is \code{method=mean}.
}
  \item{fixedmatrix}{
If set to \code{TRUE}, a coverage matrix with single nucleotide resolution will be
returned for any number of transcripts centered around precise anchor points in
a genome annotation, such a stop/start codons or transcription start sites. For instance, a 
matrix with coverage information 20bps upstream and downstream of the stop/start codons can be 
obtained with \code{fixedmatrix=TRUE, upstream=20, downstream=20} along with a \code{grl} 
instance containing the CDS exon ranges required for this operation, e.g. generated 
with \code{cdsBy(txdb, "tx")}.
}
  \item{resizefeatures}{
Needs to be set to \code{TRUE} when \code{fixedmatrix=TRUE}. Internally, this will use
the \code{systemPipeR::.resizeFeature} function to extend single and multi component
features at their most left and most right end coordinates. The corresponding extension 
values are specified under the \code{upstream} and \code{downstream} arguments.
}
  \item{upstream}{
Single positive integer specifying the upstream extension length relative to the orientation
of each feature in the genome. More details are given above.
}
  \item{downstream}{
Single positive integer specifying the downstream extension length relative to the orientation
of each feature in the genome. More details are given above.
}
  \item{outfile}{
Default \code{NULL} omits writing of the results to a file. If a file name is specified then
the results are written to a tabular file. If \code{bfl} contains the paths to several BAM files
then the results will be appended to the same file where the first column specifies the
sample labels. Redirecting the results to file is particularly useful when processing large 
files of many sample where computation times can be significant.
}
  \item{overwrite}{
If set to \code{TRUE} any existing file assigned to \code{outfile} will be overwritten.
}
}
\value{
The function allows to return the following four distinct outputs. The settings to return
these instances are illustrated below in the example section.
\item{(A)}{\code{data.frame} containing binned coverage where rows are features and columns 
coverage bins. The first four columns contain (i) the sample names, (ii) the number of 
total aligned reads in the corresponding BAM files (useful for normalization), (iii) the
feature IDs, (iv) strand of the coverage. All following columns are numeric and contain
the actual coverage data for the sense and antisense strand of each feature. 
}
\item{(B)}{\code{data.frame} containing coverage with single nucleotide
resolution around anchor points such as start and stop codons. The two matrix
components are appended column-wise. To clearly distinguish the two data components, 
they are separated by a specialty column containing pipe characters. The first four 
columns are the same as described under (A). The column title for the
anchor point is 0. For instance, if the features are CDSs then the first 0 corresponds to the
first nucleotide of the start codon and the second 0 to the last nucleotide of
the stop codon. Upstream and downstream positions are indicated by negative and
positive column numbers, respectively.
} 
\item{(C)}{\code{data.frame} containing combined results of (A) and (B) where the
first set of columns contains to the coverage around the start codons, the second one
the binned coverage of the CDSs and the third one the coverage around the stop codons 
separated by the same pipe columns mentioned under (B).
}
\item{(D)}{\code{Rle} list containing the nucleotide level coverage of each feature
}
}
\author{
Thomas Girke
}
\seealso{
\code{plotfeatureCoverage}
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

## (A) Generate binned coverage for two BAM files and 4 transcripts
grl <- cdsBy(txdb, "tx", use.names=TRUE)
fcov <- featureCoverage(bfl=BamFileList(outpaths(args)[1:2]), grl=grl[1:4], resizereads=NULL, 
                    readlengthrange=NULL, Nbins=20, method=mean, fixedmatrix=FALSE, 
                    resizefeatures=TRUE, upstream=20, downstream=20)

## (B) Coverage matrix upstream and downstream of start/stop codons
grl <- cdsBy(txdb, "tx", use.names=TRUE)
fcov <- featureCoverage(bfl=BamFileList(outpaths(args)[1:2]), grl=grl[1:4], resizereads=NULL, 
                    readlengthrange=NULL, Nbins=20, method=mean, fixedmatrix=TRUE, 
                    resizefeatures=TRUE, upstream=20, downstream=20)

## (C) Combined matrix for both binned and start/stop codon 
grl <- cdsBy(txdb, "tx", use.names=TRUE)
fcov <- featureCoverage(bfl=BamFileList(outpaths(args)[1:2]), grl=grl[1:4], resizereads=NULL, 
                    readlengthrange=NULL, Nbins=20, method=mean, fixedmatrix=TRUE, 
                    resizefeatures=TRUE, upstream=20, downstream=20)

## (D) Rle coverage objects one for each query feature
grl <- cdsBy(txdb, "tx", use.names=TRUE)
fcov <- featureCoverage(bfl=BamFileList(outpaths(args)[1:2]), grl=grl[1:4], resizereads=NULL, 
                    readlengthrange=NULL, Nbins=NULL, method=mean, fixedmatrix=FALSE, 
                    resizefeatures=TRUE, upstream=20, downstream=20)

}
}
\keyword{ utilities }
