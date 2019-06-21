---
title: 3. Read preprocessing
last_updated: Fri Jun 21 16:34:15 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_03.html
---

## Quality filtering and adaptor trimming

The following custom function trims adaptors hierarchically from the longest to
the shortest match of the right end of the reads. If `internalmatch=TRUE` then internal matches will trigger the same behavior.  The argument `minpatternlength` defines the shortest adaptor match to consider in this iterative process. In addition, the function
removes reads containing Ns or homopolymer regions. More detailed information
on read preprocessing is provided in `systemPipeR's` main vignette.


```r
args <- systemArgs(sysma = "param/trim.param", mytargets = "targets.txt")
fctpath <- system.file("extdata", "custom_Fct.R", package = "systemPipeR")
source(fctpath)
iterTrim <- ".iterTrimbatch1(fq, pattern='ACACGTCT', internalmatch=FALSE, minpatternlength=6, Nnumber=1, polyhomo=50, minreadlength=16, maxreadlength=101)"
preprocessReads(args = args, Fct = iterTrim, batchsize = 1e+05, 
    overwrite = TRUE, compress = TRUE)
writeTargetsout(x = args, file = "targets_trim.txt", overwrite = TRUE)
```

## FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of
useful quality statistics for a set of FASTQ files including per cycle quality
box plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads above
quality cutoffs and mean quality distribution. The results are written to a PDF file named `fastqReport.png`.


```r
args <- systemArgs(sysma = NULL, mytargets = "targets_trim.txt")
fqlist <- seeFastq(fastq = infile1(args), batchsize = 1e+05, 
    klength = 8)
png("./results/fastqReport.png", height = 18, width = 4 * length(fqlist), 
    units = "in", res = 72)
seeFastqPlot(fqlist)
dev.off()
```

![](./pages/mydoc/systemPipeRIBOseq_files/fastqReport.png)
<div align="center"><b>Figure 1:</b> FASTQ quality report. To zoom in, rigth click image and open it in a separate browser tab. </div>

<br><br><center><a href="mydoc_systemPipeRIBOseq_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
