---
title: 3. Read preprocessing
last_updated: Fri Jun 21 16:33:06 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_03.html
---

## Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample comparisons of the analysis workflow.


```r
targetspath <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4, 1:4]
```

```
##                  FileName1                FileName2
## 1 ./data/SRR446027_1.fastq ./data/SRR446027_2.fastq
## 2 ./data/SRR446028_1.fastq ./data/SRR446028_2.fastq
## 3 ./data/SRR446029_1.fastq ./data/SRR446029_2.fastq
## 4 ./data/SRR446030_1.fastq ./data/SRR446030_2.fastq
##   SampleName Factor
## 1        M1A     M1
## 2        M1B     M1
## 3        A1A     A1
## 4        A1B     A1
```

## Read quality filtering and trimming

The following removes reads with low quality base calls (here Phred
scores below 20) from all FASTQ files.


```r
args <- systemArgs(sysma = "param/trimPE.param", mytargets = "targetsPE.txt")[1:4]
# Note: subsetting!
filterFct <- function(fq, cutoff = 20, Nexceptions = 0) {
    qcount <- rowSums(as(quality(fq), "matrix") <= cutoff, na.rm = TRUE)
    fq[qcount <= Nexceptions]
    # Retains reads where Phred scores are >= cutoff with N
    # exceptions
}
preprocessReads(args = args, Fct = "filterFct(fq, cutoff=20, Nexceptions=0)", 
    batchsize = 1e+05)
writeTargetsout(x = args, file = "targets_PEtrim.txt", overwrite = TRUE)
```

## FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of 
useful quality statistics for a set of FASTQ files including per cycle quality box
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

![](./pages/mydoc/systemPipeVARseq_files/fastqReport.png)
<div align="center">Figure 1: FASTQ quality report for 18 samples</div>

<br><br><center><a href="mydoc_systemPipeVARseq_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
