\name{scaleRanges}
\alias{scaleRanges}
\title{
Scale spliced ranges to genome coordinates
}
\description{
Function to scale mappings of spliced features (query ranges) to their
corresponding genome coordinates (subject ranges). The method accounts for
introns in the subject ranges that are absent in the query ranges. A use case
example are uORFs predicted in the 5' UTRs sequences using \code{predORF}.
These query ranges are given relative to the 5' UTR sequence. The
\code{scaleRanges} function will scale them to the corresponding genome
coordinates. This way they can be used in RNA-Seq expression experiments
like other gene ranges. 
}
\usage{
scaleRanges(subject, query, type = "custom", verbose = TRUE)
}
\arguments{
  \item{subject}{
Genomic ranges provided as \code{GRangesList} object. Their name and length
requirements are described under \code{query}.
}
  \item{query}{
Feature level ranges provided as \code{GRangesList} object. The names of the 
query ranges need to match the names of the \code{GRangesList} object 
assigened to the \code{subject} argument. In addition, the length of 
each \code{query} range cannot exceed the total length of the corresponding 
\code{subject} range set. 
}
  \item{type}{
Feature name to use in \code{type} column of \code{GRangesList} result.  
}
  \item{verbose}{
The setting \code{verbose=FALSE} suppresses all print messages.
}
}
\value{
Object of class \code{GRangesList}
}
\author{
Thomas Girke
}

\seealso{
\code{predORF}
}
\examples{
library(IRanges)
## Usage for simple example
subject <- GRanges(seqnames="Chr1", IRanges::IRanges(c(5,15,30),c(10,25,40)), strand="+") 
query <- GRanges(seqnames="myseq", IRanges::IRanges(1, 9), strand="+")
scaleRanges(GRangesList(myid1=subject), GRangesList(myid1=query), type="test")

\dontrun{
## Usage for more complex example
library(GenomicFeatures); library(systemPipeRdata)
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
