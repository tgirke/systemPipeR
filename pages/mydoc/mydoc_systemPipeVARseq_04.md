---
title: 4. Alignments
last_updated: Sat Feb  2 12:30:52 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_04.html
---

## Read mapping with `BWA-MEM` 

The NGS reads of this project are aligned against the reference genome
sequence using the highly variant tolerant short read aligner `BWA-MEM`
(Li , 2013; Li et al., 2009). The parameter settings of the aligner are
defined in the `bwa.param` file.


```r
args <- systemArgs(sysma = "param/bwa.param", mytargets = "targets.txt")
sysargs(args)[1]  # Command-line parameters for first FASTQ file
```

Runs the alignments sequentially (_e.g._ on a single machine)


```r
moduleload(modules(args))
system("bwa index -a bwtsw ./data/tair10.fasta")
bampaths <- runCommandline(args = args)
writeTargetsout(x = args, file = "targets_bam.txt", overwrite = TRUE)
```

Alternatively, the alignment jobs can be submitted to a compute cluster,
here using 72 CPU cores (18 `qsub` processes each with 4 CPU cores).


```r
moduleload(modules(args))
system("bwa index -a bwtsw ./data/tair10.fasta")
resources <- list(walltime = 120, ntasks = 1, ncpus = cores(args), 
    memory = 1024)
reg <- clusterRun(args, conffile = ".batchtools.conf.R", Njobs = 18, 
    template = "batchtools.slurm.tmpl", runid = "01", resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
writeTargetsout(x = args, file = "targets_bam.txt", overwrite = TRUE)
```

Check whether all BAM files have been created


```r
file.exists(outpaths(args))
```

## Read mapping with `gsnap` 

An alternative variant tolerant aligner is `gsnap` from the `gmapR` package
(Wu et al., 2010). The following code shows how to run this aligner on
multiple nodes of a computer cluster that uses Torque as scheduler.


```r
library(gmapR)
library(BiocParallel)
library(batchtools)
args <- systemArgs(sysma = "param/gsnap.param", mytargets = "targetsPE.txt")
gmapGenome <- GmapGenome(systemPipeR::reference(args), directory = "data", 
    name = "gmap_tair10chr", create = TRUE)
f <- function(x) {
    library(gmapR)
    library(systemPipeR)
    args <- systemArgs(sysma = "param/gsnap.param", mytargets = "targetsPE.txt")
    gmapGenome <- GmapGenome(reference(args), directory = "data", 
        name = "gmap_tair10chr", create = FALSE)
    p <- GsnapParam(genome = gmapGenome, unique_only = TRUE, 
        molecule = "DNA", max_mismatches = 3)
    o <- gsnap(input_a = infile1(args)[x], input_b = infile2(args)[x], 
        params = p, output = outfile1(args)[x])
}
resources <- list(walltime = 120, ntasks = 1, ncpus = cores(args), 
    memory = 1024)
param <- BatchtoolsParam(workers = 4, cluster = "slurm", template = "batchtools.slurm.tmpl", 
    resources = resources)
d <- bplapply(seq(along = args), f, BPPARAM = param)
writeTargetsout(x = args, file = "targets_gsnap_bam.txt", overwrite = TRUE)
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
symLink2bam(sysargs = args, htmldir = c("~/.html/", "projects/gen242/"), 
    urlbase = "http://biocluster.ucr.edu/~tgirke/", urlfile = "./results/IGVurl.txt")
```

<br><br><center><a href="mydoc_systemPipeVARseq_03.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_05.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
