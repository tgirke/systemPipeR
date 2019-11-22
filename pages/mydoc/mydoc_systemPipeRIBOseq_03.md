---
title: 3. Read preprocessing
last_updated: Thu Nov 21 16:49:12 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_03.html
---

## Quality filtering and adaptor trimming

The following custom function trims adaptors hierarchically from the longest to
the shortest match of the right end of the reads. If `internalmatch=TRUE` then internal matches will trigger the same behavior. The argument `minpatternlength` defines the shortest adaptor match to consider in this iterative process. In addition, the function removes reads containing Ns or homopolymer regions. More detailed information on read preprocessing is provided in `systemPipeR's` main vignette.

First, we construct _`SYSargs2`_ object from _`cwl`_ and _`yml`_ param and _`targets`_ files.


```r
dir_path <- system.file("extdata/cwl/preprocessReads/trim-se", 
    package = "systemPipeR")
trim <- loadWorkflow(targets = targetspath, wf_file = "trim-se.cwl", 
    input_file = "trim-se.yml", dir_path = dir_path)
trim <- renderWF(trim, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
trim
output(trim)[1:2]
```

Next, we execute the code for trimming all the raw data. 


```r
fctpath <- system.file("extdata", "custom_Fct.R", package = "systemPipeR")
source(fctpath)
iterTrim <- ".iterTrimbatch1(fq, pattern='ACACGTCT', internalmatch=FALSE, minpatternlength=6, Nnumber=1, polyhomo=50, minreadlength=16, maxreadlength=101)"
preprocessReads(args = trim, Fct = iterTrim, batchsize = 1e+05, 
    overwrite = TRUE, compress = TRUE)
writeTargetsout(x = trim, file = "targets_trim.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```

## FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of
useful quality statistics for a set of FASTQ files including per cycle quality
box plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads above
quality cutoffs and mean quality distribution. The results are written to a PDF file named `fastqReport.png`.


```r
fqlist <- seeFastq(fastq = infile1(trim), batchsize = 10000, 
    klength = 8)
png("./results/fastqReport.png", height = 18, width = 4 * length(fqlist), 
    units = "in", res = 72)
seeFastqPlot(fqlist)
dev.off()
```

![](./pages/mydoc/systemPipeRIBOseq_files/fastqReport.png)
<div align="center"><b>Figure 1:</b> FASTQ quality report. To zoom in, right click image and open it in a separate browser tab. </div>

<br><br><center><a href="mydoc_systemPipeRIBOseq_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
