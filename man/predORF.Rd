\name{predORF}
\alias{predORF}
\title{
Predict ORFs 
}
\description{
Predicts open reading frames (ORFs) and coding sequences (CDSs) in DNA sequences provided as \code{DNAString} or \code{DNAStringSet} objects. 
}
\usage{
predORF(x, n = 1, type = "grl", mode = "orf", strand = "sense", longest_disjoint=FALSE, startcodon = "ATG", stopcodon = c("TAA", "TAG", "TGA"))
}
\arguments{
  \item{x}{
DNA query sequence(s) provided as \code{DNAString} or \code{DNAStringSet} object.
}
  \item{n}{
Defines the maximum number of ORFs to return for each input sequence. The ORFs identified are sorted decreasingly 
by their length. For instance, \code{n=1} (default) returns the longest ORF, \code{n=2} the two longest ones, and so on. 
}
  \item{type}{
One of three options provided as character values: \code{'df'} returns results as \code{data.frame}, while \code{'gr'} and \code{'grl'} (default) return them as \code{GRanges} or \code{GRangesList} objects, respectively.
}
  \item{mode}{
The setting \code{mode='ORF'} returns a continuous reading frame that begins with a start codon and ends with a stop codon. The setting \code{mode='CDS'} 
return continuous reading frames that do not need to begin or end with start or stop codons, respectively. 
}
  \item{strand}{
One of three options passed on as character vector of length one: \code{'sense'} performs the predictions only for the sense strand of the query sequence(s), 
\code{'antisense'} does it only for the antisense strand and \code{'both'} does it for both strands.
}
  \item{longest_disjoint}{
    If set to \code{TRUE} and \code{n='all'}, the results will be subsetted to non-overlapping ORF set containing longest ORF. 
}
  \item{startcodon}{
Defines the start codon(s) for ORF predictions. The default is set to the standard start codon 'ATG'. Any custom set of triplet DNA sequences can be assigned here.
}
  \item{stopcodon}{
Defines the stop codon(s) for ORF predictions. The default is set to the three standard stop codons 'TAA', 'TAG' and 'TGA'. Any custom set of triplet DNA sequences can be assigned here.
}
}
\value{
Returns ORF/CDS ranges identified in query sequences as \code{GRanges} or
\code{data.frame} object. The \code{type} argument defines which one of them
will be returned. The objects contain the following columns:
    \itemize{
        \item{\code{seqnames}: names of query sequences}
        \item{\code{subject_id}: identified ORF/CDS ranges numbered by query} 
        \item{\code{start/end}: start and end positions of ORF/CDS ranges}
        \item{\code{strand}: strand of query sequence used for prediction}
        \item{\code{width}: length of subject range in bases}
        \item{\code{inframe2end}: frame of identified ORF/CDS relative to 3'
end of query sequence. This can be important if the query sequence was
extracted directly upstream of an ORF (e.g. 5' UTR upstream of main ORF). 
The value 1 stands for in-frame with downstream ORF, while 2 or 3 indicates 
a shift of one or two bases, respectively. } }
}
\author{
Thomas Girke
}

\seealso{
\code{scaleRanges}
}
\examples{
## Load DNA sample data set from Biostrings package
file <- system.file("extdata", "someORF.fa", package="Biostrings")
dna <- readDNAStringSet(file)

## Predict longest ORF for sense strand in each query sequence
(orf <- predORF(dna[1:4], n=1, type="gr", mode="orf", strand="sense"))

\dontrun{
## Usage for more complex example
library(txdbmaker); library(systemPipeRdata)
gff <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
txdb <- makeTxDbFromGFF(file=gff, format="gff3", organism="Arabidopsis")
futr <- fiveUTRsByTranscript(txdb, use.names=TRUE)
genome <- system.file("extdata/annotation", "tair10.fasta", package="systemPipeRdata")
dna <- extractTranscriptSeqs(FaFile(genome), futr)
uorf <- predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="sense")
grl_scaled <- scaleRanges(subject=futr, query=uorf, type="uORF", verbose=TRUE)
export.gff3(unlist(grl_scaled), "uorf.gff")
}
}
\keyword{ utilities }
