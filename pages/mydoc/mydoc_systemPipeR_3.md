---
title: 3. How to run a Workflow
last_updated: Sat May  2 14:44:00 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeR_3.html
---

This tutorial introduces the basic ideas and tools needed to build a specific workflow from preconfigured templates.

## Load sample data and workflow templates


```r
library(systemPipeRdata)
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
```

## Setup and Requirements 

To go through this tutorial, you need the following software installed:

* R/>=3.6.2 
* systemPipeR R package (version 1.22)
* Hisat2/2.1.0

If you desire to build your pipeline with any different software, make sure to have the respective software installed and configured in your PATH. To make sure if the configuration is right, you always can test as follow:


```r
tryCL(command = "hisat2")  ## 'All set up, proceed!'
```

## Project Initialization 

The Project management structure is essential, especially for reproducibility and efficiency in the analysis. Here we show how to construct an instance of this S4 object class by the _`initWF`_ function. The object of class _`SYSarsgsList`_ storing all the configuration information for the project and allows management and control at a high level. 


```r
script <- "systemPipeRNAseq.Rmd"
targetspath <- "targets.txt"
sysargslist <- initWF(script = script, targets = targets)
```

## Project Initialization in a Temporary Directory


```r
library(systemPipeRdata)
script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package = "systemPipeRdata")
targets <- system.file("extdata", "targets.txt", package = "systemPipeR")
dir_path <- tempdir()
SYSconfig <- initProject(projPath = dir_path, targets = targets, script = script, 
    overwrite = TRUE)
sysargslist_temp <- initWF(sysconfig = "SYSconfig.yml")
```

## Configuration and run of the project


```r
sysargslist <- configWF(x = sysargslist, input_steps = "1:3")
sysargslist <- runWF(sysargslist = sysargslist, steps = "ALL")
sysargslist <- runWF(sysargslist = sysargslist, steps = "1:2")
```

<br><br><center><a href="mydoc_systemPipeR_2.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeR_4.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
