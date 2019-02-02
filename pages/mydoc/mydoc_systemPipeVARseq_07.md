---
title: 7. Annotate filtered variants
last_updated: Sat Feb  2 12:30:52 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_07.html
---

The function `variantReport` generates a variant report using
utilities provided by the `VariantAnnotation` package. The report for
each sample is written to a tabular file containing genomic context annotations
(_e.g._ coding or non-coding SNPs, amino acid changes, IDs of affected
genes, etc.) along with confidence statistics for each variant. The parameter
file `annotate_vars.param` defines the paths to the input and output
files which are stored in a new `SYSargs` instance. 

## Basics of annotating variants

Variants overlapping with common annotation features can be identified with `locateVariants`.

```r
library("GenomicFeatures")
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_gatk_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
vcf <- readVcf(infile1(args)[1], "A. thaliana")
locateVariants(vcf, txdb, CodingVariants())
```

Synonymous/non-synonymous variants of coding sequences are computed by the predictCoding function for variants overlapping with coding regions.


```r
fa <- FaFile(systemPipeR::reference(args))
predictCoding(vcf, txdb, seqSource = fa)
```

## Annotate filtered variants called by `GATK`


```r
library("GenomicFeatures")
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_gatk_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args = args, txdb = txdb, fa = fa, 
    organism = "A. thaliana"))
```

## Annotate filtered variants called by `BCFtools`


```r
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_sambcf_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args = args, txdb = txdb, fa = fa, 
    organism = "A. thaliana"))
```

## Annotate filtered variants called by `VariantTools`


```r
args <- systemArgs(sysma = "param/annotate_vars.param", mytargets = "targets_vartools_filtered.txt")
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(systemPipeR::reference(args))
suppressAll(variantReport(args = args, txdb = txdb, fa = fa, 
    organism = "A. thaliana"))
```

View annotation result for single sample

```r
read.delim(outpaths(args)[1])[38:40, ]
```

<br><br><center><a href="mydoc_systemPipeVARseq_06.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_08.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
