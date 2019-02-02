---
title: 4. Alignments
last_updated: Sat Feb  2 11:43:49 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeChIPseq_04.html
---

## Read mapping with `Bowtie2` 

The NGS reads of this project will be aligned with `Bowtie2` against the
reference genome sequence (Langmead et al., 2012). The parameter settings of the
aligner are defined in the `bowtieSE.param` file. In ChIP-Seq experiments it is
usually more appropriate to eliminate reads mapping to multiple locations. To
achieve this, users want to remove the argument setting `-k 50 non-deterministic` 
in the `bowtieSE.param` file.

The following submits 18 alignment jobs via a scheduler to a computer cluster.


```r
args <- systemArgs(sysma = "param/bowtieSE.param", mytargets = "targets_chip_trim.txt")
sysargs(args)[1]  # Command-line parameters for first FASTQ file
moduleload(modules(args))  # Skip if a module system is not used
system("bowtie2-build ./data/tair10.fasta ./data/tair10.fasta")
# Indexes reference genome
resources <- list(walltime = 120, ntasks = 1, ncpus = cores(args), 
    memory = 1024)
reg <- clusterRun(args, conffile = ".batchtools.conf.R", Njobs = 18, 
    template = "batchtools.slurm.tmpl", runid = "01", resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
writeTargetsout(x = args, file = "targets_bam.txt", overwrite = TRUE)
```

Alternatively, one can run the alignments sequentially on a single system. 


```r
runCommandline(args)
```

Check whether all BAM files have been created

```r
file.exists(outpaths(args))
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
specified under `urlfile`, here `IGVurl.txt`.


```r
symLink2bam(sysargs = args, htmldir = c("~/.html/", "somedir/"), 
    urlbase = "http://biocluster.ucr.edu/~tgirke/", urlfile = "./results/IGVurl.txt")
```

<br><br><center><a href="mydoc_systemPipeChIPseq_03.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeChIPseq_05.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
