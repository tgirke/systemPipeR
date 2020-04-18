---
title: 4. Alignments
last_updated: Sat Apr 18 12:33:39 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_04.html
---

## Read mapping with `HISAT2`

The following steps will demonstrate how to use the short read aligner `Hisat2`
(Kim et al., 2015) in both interactive job submissions and batch submissions to
queuing systems of clusters using the _`systemPipeR's`_ new CWL command-line interface.

Build `Hisat2` index.


```r
dir_path <- system.file("extdata/cwl/hisat2/hisat2-idx", package = "systemPipeR")
idx <- loadWorkflow(targets = NULL, wf_file = "hisat2-index.cwl", 
    input_file = "hisat2-index.yml", dir_path = dir_path)
idx <- renderWF(idx)
idx
cmdlist(idx)

## Run
runCommandline(idx, make_bam = FALSE)
```

The parameter settings of the aligner are defined in the `hisat2-mapping-se.cwl` 
and `hisat2-mapping-se.yml` files. The following shows how to construct the 
corresponding *SYSargs2* object, here *args*.


```r
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package = "systemPipeR")
args <- loadWorkflow(targets = targetspath, wf_file = "hisat2-mapping-se.cwl", 
    input_file = "hisat2-mapping-se.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
args
cmdlist(args)[1:2]
output(args)[1:2]

## Run
args <- runCommandline(args)
```


```r
library(batchtools)
resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, 
    make_bam = TRUE, dir = FALSE), conffile = ".batchtools.conf.R", 
    template = "batchtools.slurm.tmpl", Njobs = 18, runid = "01", 
    resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
```

Check whether all BAM files have been created.


```r
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
file.exists(outpaths)
```

## Read and alignment stats

The following provides an overview of the number of reads in each sample and how many of them aligned to the reference.


```r
read_statsDF <- alignStats(args = args)
write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, 
    quote = FALSE, sep = "\t")
```


```r
read.table(system.file("extdata", "alignStats.xls", package = "systemPipeR"), 
    header = TRUE)[1:4, ]
```

```
##   FileName Nreads2x Nalign Perc_Aligned Nalign_Primary
## 1      M1A   192918 177961     92.24697         177961
## 2      M1B   197484 159378     80.70426         159378
## 3      A1A   189870 176055     92.72397         176055
## 4      A1B   188854 147768     78.24457         147768
##   Perc_Aligned_Primary
## 1             92.24697
## 2             80.70426
## 3             92.72397
## 4             78.24457
```

## Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV. The corresponding URLs are written to a file
with a path specified under `urlfile` in the `results` directory.


```r
symLink2bam(sysargs = args, htmldir = c("~/.html/", "projects/tests/"), 
    urlbase = "http://biocluster.ucr.edu/~tgirke/", urlfile = "./results/IGVurl.txt")
```

<br><br><center><a href="mydoc_systemPipeRIBOseq_03.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_05.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
