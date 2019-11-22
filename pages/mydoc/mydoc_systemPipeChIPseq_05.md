---
title: 5. Utilities for coverage data
last_updated: Thu Nov 21 15:49:32 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeChIPseq_05.html
---

The following introduces several utilities useful for ChIP-Seq data. They are not part of the actual workflow.

## Rle object stores coverage information


```r
library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
aligns <- readGAlignments(outpaths[1])
cov <- coverage(aligns)
cov
```

## Resizing aligned reads


```r
trim(resize(as(aligns, "GRanges"), width = 200))
```

## Naive peak calling


```r
islands <- slice(cov, lower = 15)
islands[[1]]
```

## Plot coverage for defined region


```r
library(ggbio)
myloc <- c("Chr1", 1, 1e+05)
ga <- readGAlignments(outpaths[1], use.names = TRUE, param = ScanBamParam(which = GRanges(myloc[1], 
    IRanges(as.numeric(myloc[2]), as.numeric(myloc[3])))))
autoplot(ga, aes(color = strand, fill = strand), facets = strand ~ 
    seqnames, stat = "coverage")
```

<br><br><center><a href="mydoc_systemPipeChIPseq_04.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeChIPseq_06.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
