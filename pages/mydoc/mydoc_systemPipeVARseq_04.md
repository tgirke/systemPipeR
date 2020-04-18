---
title: 4. Alignments
last_updated: Sat Apr 18 12:30:50 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_04.html
---

## Read mapping with `BWA-MEM` 

The NGS reads of this project are aligned against the reference genome
sequence using the highly variant tolerant short read aligner `BWA-MEM`
(Li , 2013; Li et al., 2009). The parameter settings of the aligner are
defined in the `gatk/bwa-pe.cwl`

DNA sequencing nowadays are usually solid for base quality and therefore, 
trimming is usually not needed in most cases. Also, variant calling tools like 
`GATK` will automatically not consider low quality bases. Therefore, this test 
code uses untrimmed fastqs. However, it is best to test with `FASTQ quality report` 
function provided above to verify on your real data first.

### Build index and dictionary files for BWA and GATK

The following object `dir_path` is the folder where all `BWA` and `GATK` param files are located.


```r
dir_path <- system.file("extdata/cwl/gatk", package = "systemPipeR")
```

Build the index and dictionary files for BWA and GATK to run. 


```r
## Index for BWA
args <- loadWorkflow(targets = NULL, wf_file = "bwa-index.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args)
cmdlist(args)  # shows the command
output(args)  # shows the expected output files
# Run single Machine
runCommandline(args, make_bam = FALSE)

## Index needed for gatk tools
args <- loadWorkflow(wf_file = "fasta_dict.cwl", input_file = "gatk.yaml", 
    dir_path = dir_path)
args <- renderWF(args)
args <- runCommandline(args, make_bam = FALSE)

## Index
args <- loadWorkflow(wf_file = "fasta_faidx.cwl", input_file = "gatk.yaml", 
    dir_path = dir_path)
args <- renderWF(args)
args <- runCommandline(args, make_bam = FALSE)
```

### Run the read mapping


```r
targetsPE <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
args <- loadWorkflow(targets = targetsPE, wf_file = "bwa-pe.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
cmdlist(args)[1:2]
output(args)[1:2]
```

Runs the alignments sequentially (_e.g._ on a single machine) by `runCommandline` function.


```r
args <- runCommandline(args = args, make_bam = FALSE)
writeTargetsout(x = args[1:2], file = "./results/targetsPE.txt", 
    step = 1, new_col = "BWA_SAM", new_col_output_index = 1, 
    overwrite = TRUE)
```

Alternatively, the alignment jobs can be submitted to a compute cluster. Here is the 
example to run cluster jobs by `clusterRun` on a `slurm` based system. 4 cpus for 
each task for 18 samples, totally 72 cpus.


```r
library(batchtools)
resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, 
    dir = FALSE, make_bam = FALSE), conffile = ".batchtools.conf.R", 
    template = "batchtools.slurm.tmpl", Njobs = 18, runid = "01", 
    resourceList = resources)
getStatus(reg = reg)
writeTargetsout(x = args, file = "./results/targetsPE.txt", step = 1, 
    new_col = "BWA_SAM", new_col_output_index = 1, overwrite = TRUE)
```

Check whether all BAM files have been created.


```r
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
file.exists(outpaths)
```

## Read and alignment stats

The following generates a summary table of the number of reads in each
sample and how many of them aligned to the reference.


```r
read_statsDF <- alignStats(args = args)
write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, 
    quote = FALSE, sep = "\t")
```

## Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV. The corresponding URLs are written to a file
with a path specified under `urlfile`, here `IGVurl.txt`.


```r
symLink2bam(sysargs = args, htmldir = c("~/.html/", "somedir/"), 
    urlbase = "http://cluster.hpcc.ucr.edu/~tgirke/", urlfile = "IGVurl.txt")
```

<br><br><center><a href="mydoc_systemPipeVARseq_03.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_05.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
