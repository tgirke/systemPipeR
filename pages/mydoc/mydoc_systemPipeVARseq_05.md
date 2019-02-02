---
title: 5. Variant calling
last_updated: Sat Feb  2 12:30:52 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_05.html
---

The following performs variant calling with `GATK`, `BCFtools` and `VariantTools` 
in parallel mode on a compute cluster (McKenna et al., 2010; Li , 2011). If a cluster is not
available, the `runCommandline` function can be used to run the variant calling with `GATK` 
and `BCFtools` for each sample sequentially on a single machine, or `callVariants` in case 
of `VariantTools`. Typically, the user would choose here only one variant caller rather
than running several ones.

## Variant calling with `GATK`

The following creates in the inital step a new `targets` file
(`targets_bam.txt`). The first column of this file gives the paths to
the BAM files created in the alignment step. The new targets file and the
parameter file `gatk.param` are used to create a new `SYSargs`
instance for running GATK. Since GATK involves many processing steps, it is
executed by a bash script `gatk_run.sh` where the user can specify the
detailed run parameters. All three files are expected to be located in the
current working directory. Samples files for `gatk.param` and
`gatk_run.sh` are available in the `param` subdirectory
provided by `systemPipeRdata`.


```r
moduleload("picard/1.130")
moduleload("samtools/1.3")
system("picard CreateSequenceDictionary R=./data/tair10.fasta O=./data/tair10.dict")
system("samtools faidx data/tair10.fasta")
args <- systemArgs(sysma = "param/gatk.param", mytargets = "targets_bam.txt")
resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
reg <- clusterRun(args, conffile = ".batchtools.conf.R", Njobs = 18, 
    template = "batchtools.slurm.tmpl", runid = "01", resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
# unlink(outfile1(args), recursive = TRUE, force = TRUE)
writeTargetsout(x = args, file = "targets_gatk.txt", overwrite = TRUE)
```

## Variant calling with `BCFtools`

The following runs the variant calling with `BCFtools`. This step requires
in the current working directory the parameter file `sambcf.param` and the bash script 
`sambcf_run.sh`.


```r
args <- systemArgs(sysma = "param/sambcf.param", mytargets = "targets_bam.txt")
resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
reg <- clusterRun(args, conffile = ".batchtools.conf.R", Njobs = 18, 
    template = "batchtools.slurm.tmpl", runid = "01", resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
# unlink(outfile1(args), recursive = TRUE, force = TRUE)
writeTargetsout(x = args, file = "targets_sambcf.txt", overwrite = TRUE)
```

## Variant calling with `VariantTools`  


```r
library(gmapR)
library(BiocParallel)
library(batchtools)
args <- systemArgs(sysma = "param/vartools.param", mytargets = "targets_gsnap_bam.txt")
f <- function(x) {
    library(VariantTools)
    library(gmapR)
    library(systemPipeR)
    args <- systemArgs(sysma = "param/vartools.param", mytargets = "targets_gsnap_bam.txt")
    gmapGenome <- GmapGenome(systemPipeR::reference(args), directory = "data", 
        name = "gmap_tair10chr", create = FALSE)
    tally.param <- TallyVariantsParam(gmapGenome, high_base_quality = 23L, 
        indels = TRUE)
    bfl <- BamFileList(infile1(args)[x], index = character())
    var <- callVariants(bfl[[1]], tally.param)
    sampleNames(var) <- names(bfl)
    writeVcf(asVCF(var), outfile1(args)[x], index = TRUE)
}
resources <- list(walltime = 120, ntasks = 1, ncpus = cores(args), 
    memory = 1024)
param <- BatchtoolsParam(workers = 4, cluster = "slurm", template = "batchtools.slurm.tmpl", 
    resources = resources)
d <- bplapply(seq(along = args), f, BPPARAM = param)
writeTargetsout(x = args, file = "targets_vartools.txt", overwrite = TRUE)
```

## Inspect VCF file 

VCF files can be imported into R with the `readVcf` function. Both `VCF` and `VRanges` objects provide
convenient data structure for working with variant data (_e.g._ SNP quality filtering). 


```r
library(VariantAnnotation)
args <- systemArgs(sysma = "param/filter_gatk.param", mytargets = "targets_gatk.txt")
vcf <- readVcf(infile1(args)[1], "A. thaliana")
vcf
vr <- as(vcf, "VRanges")
vr
```

<br><br><center><a href="mydoc_systemPipeVARseq_04.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_06.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
