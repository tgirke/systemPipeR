---
title: 3. Read preprocessing
last_updated: Sat Apr 18 12:43:59 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRNAseq_03.html
---

## Read quality filtering and trimming

The function `preprocessReads` allows to apply predefined or custom
read preprocessing functions to all FASTQ files referenced in a
`SYSargs2` container, such as quality filtering or adapter trimming
routines. The paths to the resulting output FASTQ files are stored in the
`output` slot of the `SYSargs2` object. The following example performs adapter trimming with
the `trimLRPatterns` function from the `Biostrings` package.
After the trimming step a new targets file is generated (here
`targets_trim.txt`) containing the paths to the trimmed FASTQ files.
The new targets file can be used for the next workflow step with an updated
`SYSargs2` instance, _e.g._ running the NGS alignments using the
trimmed FASTQ files.

Construct _`SYSargs2`_ object from _`cwl`_ and _`yml`_ param and _`targets`_ files.


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


```r
preprocessReads(args = trim, Fct = "trimLRPatterns(Rpattern='GCCCGGGTAA', 
                subject=fq)", 
    batchsize = 1e+05, overwrite = TRUE, compress = TRUE)
writeTargetsout(x = trim, file = "targets_trimPE.txt", step = 1, 
    new_col = "FileName1", new_col_output_index = 1, overwrite = TRUE)
```

## FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of useful 
quality statistics for a set of FASTQ files including per cycle quality box
plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads
above quality cutoffs and mean quality distribution. The results are
written to a PDF file named `fastqReport.pdf`.


```r
fqlist <- seeFastq(fastq = infile1(trim), batchsize = 10000, 
    klength = 8)
pdf("./results/fastqReport.pdf", height = 18, width = 4 * length(fqlist))
seeFastqPlot(fqlist)
dev.off()
```

![](./pages/mydoc/systemPipeRNAseq_files/fastqReport.png)
<div align="center">Figure 1: FASTQ quality report for 18 samples</div>

<br><br><center><a href="mydoc_systemPipeRNAseq_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRNAseq_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
