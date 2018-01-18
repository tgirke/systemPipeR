---
title: 2. Samples and environment settings
last_updated: Mon Nov 13 18:37:51 2017
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRNAseq_02.html
---

## Environment settings and input data

Typically, the user wants to record here the sources and versions of the
reference genome sequence along with the corresponding annotations. In
the provided sample data set all data inputs are stored in a `data`
subdirectory and all results will be written to a separate `results` directory,
while the `systemPipeRNAseq.Rmd` script and the `targets` file are expected to be located in the parent
directory. The R session is expected to run from this parent directory.

To run this sample report, mini sample FASTQ and reference genome files
can be downloaded from
[here](http://biocluster.ucr.edu/~tgirke/projects/systemPipeR_test_data.zip).
The chosen data set [SRP010938](http://www.ncbi.nlm.nih.gov/sra/?term=SRP010938)
contains 18 paired-end (PE) read sets from *Arabidposis thaliana*
(Howard et al., 2013). To minimize processing time during testing, each FASTQ
file has been subsetted to 90,000-100,000 randomly sampled PE reads that
map to the first 100,000 nucleotides of each chromosome of the *A.
thalina* genome. The corresponding reference genome sequence (FASTA) and
its GFF annotion files (provided in the same download) have been
truncated accordingly. This way the entire test sample data set is less
than 200MB in storage space. A PE read set has been chosen for this test
data set for flexibility, because it can be used for testing both types
of analysis routines requiring either SE (single end) reads or PE reads.

The following loads one of the available NGS workflow templates (here RNA-Seq)
into the user's current working directory. At the moment, the package includes
workflow templates for RNA-Seq, ChIP-Seq, VAR-Seq and Ribo-Seq. Templates for
additional NGS applications will be provided in the future.


```r
library(systemPipeRdata)
genWorkenvir(workflow="rnaseq")
setwd("rnaseq")
```

Alternatively, this can be done from the command-line as follows:


```sh
Rscript -e "systemPipeRdata::genWorkenvir(workflow='rnaseq')"
```

Now open the R markdown script `systemPipeRNAseq.Rmd`in your R IDE (_e.g._
vim-r or RStudio) and run the workflow as outlined below. If you work under
Vim-R-Tmux, the following command sequence will connect the user in an
interactive session with a node on the cluster. The code of the `Rmd`
script can then be sent from Vim on the login (head) node to an open R session running
on the corresponding computer node. This is important since Tmux sessions
should not be run on the computer nodes. 


```r
# push `F2` on your keyboard to open interactive R session
q("no") # closes R session on head node
srun --x11 --partition=intel --mem=2gb --cpus-per-task 1 --ntasks 1 --time 2:00:00 --pty bash -l
R
```

Now check whether your R session is running on a computer node of the cluster and not on a head node.


```r
system("hostname") # should return name of a compute node starting with i or c 
getwd() # checks current working directory of R session
dir() # returns content of current working directory
```

## Required packages and resources

The `systemPipeR` package needs to be loaded to perform the analysis steps shown in
this report (H Backman et al., 2016).


```r
library(systemPipeR)
```

If applicable load custom functions not provided by


```r
source("systemPipeRNAseq_Fct.R")
```
## Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample
comparisons of the analysis workflow.


```r
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")[,1:4]
targets
```

```
##                    FileName SampleName Factor SampleLong
## 1  ./data/SRR446027_1.fastq        M1A     M1  Mock.1h.A
## 2  ./data/SRR446028_1.fastq        M1B     M1  Mock.1h.B
## 3  ./data/SRR446029_1.fastq        A1A     A1   Avr.1h.A
## 4  ./data/SRR446030_1.fastq        A1B     A1   Avr.1h.B
## 5  ./data/SRR446031_1.fastq        V1A     V1   Vir.1h.A
## 6  ./data/SRR446032_1.fastq        V1B     V1   Vir.1h.B
## 7  ./data/SRR446033_1.fastq        M6A     M6  Mock.6h.A
## 8  ./data/SRR446034_1.fastq        M6B     M6  Mock.6h.B
## 9  ./data/SRR446035_1.fastq        A6A     A6   Avr.6h.A
## 10 ./data/SRR446036_1.fastq        A6B     A6   Avr.6h.B
## 11 ./data/SRR446037_1.fastq        V6A     V6   Vir.6h.A
## 12 ./data/SRR446038_1.fastq        V6B     V6   Vir.6h.B
## 13 ./data/SRR446039_1.fastq       M12A    M12 Mock.12h.A
## 14 ./data/SRR446040_1.fastq       M12B    M12 Mock.12h.B
## 15 ./data/SRR446041_1.fastq       A12A    A12  Avr.12h.A
## 16 ./data/SRR446042_1.fastq       A12B    A12  Avr.12h.B
## 17 ./data/SRR446043_1.fastq       V12A    V12  Vir.12h.A
## 18 ./data/SRR446044_1.fastq       V12B    V12  Vir.12h.B
```

<br><br><center><a href="mydoc_systemPipeRNAseq_01.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRNAseq_03.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
