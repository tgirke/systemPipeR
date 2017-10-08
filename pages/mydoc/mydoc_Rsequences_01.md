---
title: NGS Analysis Basics <br> <br> 1. Overview
last_updated: Wed Jun  7 19:39:39 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rsequences_01.html
---
Author: Thomas Girke

Last update: 07 June, 2017 

Alternative formats of this vignette:
[ [HTML](http://girke.bioinformatics.ucr.edu/GEN242/pages/mydoc/Rsequences.html){:target="_blank"} ],
[ [PDF](http://girke.bioinformatics.ucr.edu/GEN242/pages/mydoc/Rsequences.pdf){:target="_blank"} ],
[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242/gh-pages/_vignettes/09_Rsequences/Rsequences.Rmd){:target="_blank"} ],
[ [.R](https://raw.githubusercontent.com/tgirke/GEN242/gh-pages/_vignettes/09_Rsequences/Rsequences.R){:target="_blank"} ],
[ [Old PDF Slides](https://drive.google.com/open?id=0B-lLYVUOliJFWEZ0Vkdwckt5LTQ){:target="_blank"} ]



## Sequence Analysis in R and Bioconductor

__R Base__

* Some basic string handling utilities. Wide spectrum of numeric data analysis tools.

__Bioconductor__

Bioconductor packages provide much more sophisticated string handling utilities for sequence analysis (Lawrence et al., 2013, Huber et al., 2015).

* [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html): general sequence analysis environment
* [ShortRead](http://bioconductor.org/packages/release/bioc/html/ShortRead.html): pipeline for short read data
* [IRanges](http://bioconductor.org/packages/release/bioc/html/IRanges.html): low-level infrastructure for range data
* [GenomicRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html): high-level infrastructure for range data
* [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html): managing transcript centric annotations
* [GenomicAlignments](http://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html): handling short genomic alignments
* [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html): interface to  `samtools`, `bcftools` and `tabix` 
* [BSgenome](http://bioconductor.org/packages/release/bioc/html/BSgenome.html): genome annotation data
* [biomaRt](http://bioconductor.org/packages/release/bioc/html/biomaRt.html): interface to BioMart annotations
* [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html): Annotation imports, interface to online genome browsers
* [HelloRanges](http://bioconductor.org/packages/release/bioc/html/HelloRanges.html): Bedtools semantics in Bioc's Ranges infrastructure


<br><br><center><a href="mydoc_Rsequences_01.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rsequences_02.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
