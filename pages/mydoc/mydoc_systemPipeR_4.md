---
title: 4. Workflow templates
last_updated: Sat Feb  2 11:29:18 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeR_4.html
---

## RNA-Seq sample

Load the RNA-Seq sample workflow into your current working directory.


```r
library(systemPipeRdata)
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
```

### Run workflow

Next, run the chosen sample workflow _`systemPipeRNAseq`_ ([PDF](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/rnaseq/systemPipeRNAseq.pdf?raw=true), [Rmd](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/rnaseq/systemPipeRNAseq_single.Rmd)) by executing from the command-line _`make -B`_ within the _`rnaseq`_ directory. Alternatively, one can run the code from the provided _`*.Rmd`_ template file from within R interactively. 

Workflow includes following steps:

1. Read preprocessing
    + Quality filtering (trimming)
    + FASTQ quality report
2. Alignments: _`Tophat2`_ (or any other RNA-Seq aligner)
3. Alignment stats 
4. Read counting 
5. Sample-wise correlation analysis
6. Analysis of differentially expressed genes (DEGs)
7. GO term enrichment analysis
8. Gene-wise clustering

## ChIP-Seq sample

Load the ChIP-Seq sample workflow into your current working directory.


```r
library(systemPipeRdata)
genWorkenvir(workflow = "chipseq")
setwd("chipseq")
```

### Run workflow

Next, run the chosen sample workflow _`systemPipeChIPseq_single`_ ([PDF](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/chipseq/systemPipeChIPseq.pdf?raw=true), [Rmd](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/chipseq/systemPipeChIPseq_single.Rmd)) by executing from the command-line _`make -B`_ within the _`chipseq`_ directory. Alternatively, one can run the code from the provided _`*.Rmd`_ template file from within R interactively. 

Workflow includes following steps:

1. Read preprocessing
    + Quality filtering (trimming)
    + FASTQ quality report
2. Alignments: _`Bowtie2`_ or _`rsubread`_
3. Alignment stats 
4. Peak calling: _`MACS2`_, _`BayesPeak`_ 
5. Peak annotation with genomic context
6. Differential binding analysis
7. GO term enrichment analysis
8. Motif analysis

## VAR-Seq sample 

### VAR-Seq workflow for single machine

Load the VAR-Seq sample workflow into your current working directory.


```r
library(systemPipeRdata)
genWorkenvir(workflow = "varseq")
setwd("varseq")
```

### Run workflow

Next, run the chosen sample workflow _`systemPipeVARseq_single`_ ([PDF](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/varseq/systemPipeVARseq_single.pdf?raw=true), [Rmd](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/varseq/systemPipeVARseq_single.Rmd)) by executing from the command-line _`make -B`_ within the _`varseq`_ directory. Alternatively, one can run the code from the provided _`*.Rmd`_ template file from within R interactively. 

Workflow includes following steps:

1. Read preprocessing
    + Quality filtering (trimming)
    + FASTQ quality report
2. Alignments: _`gsnap`_, _`bwa`_
3. Variant calling: _`VariantTools`_, _`GATK`_, _`BCFtools`_
4. Variant filtering: _`VariantTools`_ and _`VariantAnnotation`_
5. Variant annotation: _`VariantAnnotation`_
6. Combine results from many samples
7. Summary statistics of samples

### VAR-Seq workflow for computer cluster

The workflow template provided for this step is called _`systemPipeVARseq.Rmd`_ ([PDF](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/varseq/systemPipeVARseq.pdf?raw=true), [Rmd](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/varseq/systemPipeVARseq.Rmd)).
It runs the above VAR-Seq workflow in parallel on multiple computer nodes of an HPC system using Slurm as scheduler. 

## Ribo-Seq sample

Load the Ribo-Seq sample workflow into your current working directory.


```r
library(systemPipeRdata)
genWorkenvir(workflow = "riboseq")
setwd("riboseq")
```

### Run workflow

Next, run the chosen sample workflow _`systemPipeRIBOseq`_ ([PDF](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/riboseq/systemPipeRIBOseq.pdf?raw=true), [Rmd](https://github.com/tgirke/systemPipeRdata/blob/master/inst/extdata/workflows/ribseq/systemPipeRIBOseq_single.Rmd)) by executing from the command-line _`make -B`_ within the _`ribseq`_ directory. Alternatively, one can run the code from the provided _`*.Rmd`_ template file from within R interactively. 

Workflow includes following steps:

1. Read preprocessing
    + Adaptor trimming and quality filtering
    + FASTQ quality report
2. Alignments: _`Tophat2`_ (or any other RNA-Seq aligner)
3. Alignment stats
4. Compute read distribution across genomic features
5. Adding custom features to workflow (e.g. uORFs)
6. Genomic read coverage along transcripts
7. Read counting 
8. Sample-wise correlation analysis
9. Analysis of differentially expressed genes (DEGs)
10. GO term enrichment analysis
11. Gene-wise clustering
12. Differential ribosome binding (translational efficiency)

<br><br><center><a href="mydoc_systemPipeR_3.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeR_5.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
