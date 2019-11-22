---
title: 4. Alignments
last_updated: Thu Nov 21 15:49:32 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeChIPseq_04.html
---

## Read mapping with `Bowtie2` 

The NGS reads of this project will be aligned with `Bowtie2` against the
reference genome sequence (Langmead et al., 2012). The parameter settings of the
aligner are defined in the `bowtie2-index.cwl` and `bowtie2-index.yml` files. 
In ChIP-Seq experiments it is usually more appropriate to eliminate reads mapping 
to multiple locations. To achieve this, users want to remove the argument setting 
`-k 50 non-deterministic` in the configuration files.

Building the index:


```r
dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-idx", package = "systemPipeR")
idx <- loadWorkflow(targets = NULL, wf_file = "bowtie2-index.cwl", 
    input_file = "bowtie2-index.yml", dir_path = dir_path)
idx <- renderWF(idx)
idx
cmdlist(idx)

## Run in single machine
runCommandline(idx, make_bam = FALSE)
```

The following submits 18 alignment jobs via a scheduler to a computer cluster.


```r
targets <- system.file("extdata", "targets_chip.txt", package = "systemPipeR")
dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-se", package = "systemPipeR")
args <- loadWF(targets = targets, wf_file = "bowtie2-mapping-se.cwl", 
    input_file = "bowtie2-mapping-se.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
args
cmdlist(args)[1:2]
output(args)[1:2]
```


```r
moduleload(modules(args))  # Skip if a module system is not used
resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, 
    dir = FALSE), conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", 
    Njobs = 18, runid = "01", resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
args <- output_update(args, dir = FALSE, replace = TRUE, extension = c(".sam", 
    ".bam"))  ## Updates the output(args) to the right location in the subfolders
output(args)
```

Alternatively, one can run the alignments sequentially on a single system. 


```r
args <- runCommandline(args, force = F)
```

Check whether all BAM files have been created and write out the new targets file. 


```r
writeTargetsout(x = args, file = "targets_bam.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
file.exists(outpaths)
```

## Read and alignment stats

The following provides an overview of the number of reads in each sample
and how many of them aligned to the reference.


```r
read_statsDF <- alignStats(args = args)
write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, 
    quote = FALSE, sep = "\t")
read.delim("results/alignStats.xls")
```

## Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV without moving these large files to a local
system. The corresponding URLs are written to a file with a path
specified under `urlfile`, here `IGVurl.txt`. Please replace the directory and the user name. 


```r
symLink2bam(sysargs = args, htmldir = c("~/.html/", "somedir/"), 
    urlbase = "http://cluster.hpcc.ucr.edu/~tgirke/", urlfile = "./results/IGVurl.txt")
```

<br><br><center><a href="mydoc_systemPipeChIPseq_03.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeChIPseq_05.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
