---
title: 3. Read preprocessing
last_updated: Fri Jun 21 16:30:25 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRNAseq_03.html
---

## Read quality filtering and trimming

The function `preprocessReads` allows to apply predefined or custom
read preprocessing functions to all FASTQ files referenced in a
`SYSargs` container, such as quality filtering or adaptor trimming
routines.  The following example performs adaptor trimming with
the `trimLRPatterns` function from the `Biostrings` package.
After the trimming step a new targets file is generated (here
`targets_trim.txt`) containing the paths to the trimmed FASTQ files.
The new targets file can be used for the next workflow step with an updated
`SYSargs` instance, _e.g._ running the NGS alignments using the
trimmed FASTQ files.


```r
args <- systemArgs(sysma = "param/trim.param", mytargets = "targets.txt")
preprocessReads(args = args, Fct = "trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", 
    batchsize = 1e+05, overwrite = TRUE, compress = TRUE)
writeTargetsout(x = args, file = "targets_trim.txt", overwrite = TRUE)
```

## FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of useful 
quality statistics for a set of FASTQ files including per cycle quality box
plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads
above quality cutoffs and mean quality distribution. The results are
written to a PDF file named `fastqReport.pdf`.


```r
args <- systemArgs(sysma = "param/tophat.param", mytargets = "targets.txt")
fqlist <- seeFastq(fastq = infile1(args), batchsize = 1e+05, 
    klength = 8)
pdf("./results/fastqReport.pdf", height = 18, width = 4 * length(fqlist))
seeFastqPlot(fqlist)
dev.off()
```

![](./pages/mydoc/systemPipeRNAseq_files/fastqReport.png)
<div align="center">Figure 1: FASTQ quality report for 18 samples</div>

<br><br><center><a href="mydoc_systemPipeRNAseq_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRNAseq_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
