---
title: 6. Adding custom features to workflow
last_updated: Sun Oct 15 13:21:42 2017
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_06.html
---

## Predicting uORFs in 5' UTR regions

The function `predORF` can be used to identify open reading frames
(ORFs) and coding sequences (CDSs) in DNA sequences provided as
`DNAString` or `DNAStringSet` objects. The setting
`mode='ORF'` returns continuous reading frames that begin with a start
codon and end with a stop codon, while `mode='CDS'` returns continuous
reading frames that do not need to begin or end with start or stop codons,
respectively. Non-canonical start and stop condons are supported by allowing the
user to provide any custom set of triplets under the `startcodon` and `stopcodon`
arguments (`i.e.` non-ATG start codons). The argument `n` defines the maximum number of ORFs to return for each
input sequence (_e.g._ `n=1` returns only the longest ORF). It also
supports the identification of overlapping and nested ORFs. Alternatively, one 
can return all non-overlapping ORFs including the longest ORF for each input
sequence with `n="all"` and `longest_disjoint=TRUE`.


```r
library(systemPipeRdata); library(GenomicFeatures); library(rtracklayer)
txdb <- makeTxDbFromGFF(file="data/tair10.gff", format="gff3", organism="Arabidopsis")
futr <- fiveUTRsByTranscript(txdb, use.names=TRUE)
dna <- extractTranscriptSeqs(FaFile("data/tair10.fasta"), futr)
uorf <- predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="sense")
```

To use the predicted ORF ranges for expression analysis given genome alignments
as input, it is necessary to scale them to the corresponding genome
coordinates. The function `scaleRanges` does this by transforming the
mappings of spliced features (query ranges) to their corresponding genome
coordinates (subject ranges). The method accounts for introns in the subject
ranges that are absent in the query ranges. The above uORFs predicted in the
provided 5' UTRs sequences using `predORF` are a typical use case
for this application. These query ranges are given relative to the 5' UTR 
sequences and `scaleRanges` will convert them to the corresponding 
genome coordinates. The resulting `GRangesList` object (here `grl_scaled`) 
can be directly used for read counting.


```r
grl_scaled <- scaleRanges(subject=futr, query=uorf, type="uORF", verbose=TRUE)
export.gff3(unlist(grl_scaled), "results/uorf.gff")
```

To confirm the correctness of the obtained uORF ranges, one can parse their
corresponding DNA sequences from the reference genome with the `getSeq`
function and then translate them with the `translate` function into
proteins. Typically, the returned protein sequences should start with a
`M` (corresponding to start codon) and end with `*`
(corresponding to stop codon). The following example does this for a single uORF 
containing three exons.


```r
translate(unlist(getSeq(FaFile("data/tair10.fasta"), grl_scaled[[7]])))
```

## Adding custom features to other feature types

If required custom feature ranges can be added to the standard features
generated with the `genFeatures` function above. The following does this for the uORF ranges
predicted with `predORF`.


```r
feat <- genFeatures(txdb, featuretype="all", reduce_ranges=FALSE)
feat <- c(feat, GRangesList("uORF"=unlist(grl_scaled)))
```

## Predicting sORFs in intergenic regions
The following identifies continuous ORFs in intergenic regions. Note,
`predORF` can only identify continuous ORFs in query sequences. The
function does not identify and remove introns prior to the ORF prediction.  


```r
feat <- genFeatures(txdb, featuretype="intergenic", reduce_ranges=TRUE)
intergenic <- feat$intergenic
strand(intergenic) <- "+"
dna <- getSeq(FaFile("data/tair10.fasta"), intergenic)
names(dna) <- mcols(intergenic)$feature_by
sorf <- predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="both")
sorf <- sorf[width(sorf) > 60] # Remove sORFs below length cutoff, here 60bp
intergenic <- split(intergenic, mcols(intergenic)$feature_by)
grl_scaled_intergenic <- scaleRanges(subject=intergenic, query=sorf, type="sORF", verbose=TRUE)
export.gff3(unlist(grl_scaled_intergenic), "sorf.gff")
translate(getSeq(FaFile("data/tair10.fasta"), unlist(grl_scaled_intergenic)))
```

<br><br><center><a href="mydoc_systemPipeRIBOseq_05.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_07.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
