---
title: systemPipeR - NGS workflow and report generation environment  <br> <br> 1. Introduction
last_updated: Fri Jun 21 16:39:14 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeR_1.html
---
Author: Daniela Cassol (danielac@ucr.edu) and Thomas Girke (thomas.girke@ucr.edu)

Last update: 21 June, 2019 

Alternative formats of this tutorial:
[ [HTML](http://girke.bioinformatics.ucr.edu/systemPipeR/pages/mydoc/systemPipeR.html){:target="_blank"} ],
[ [PDF](http://girke.bioinformatics.ucr.edu/systemPipeR/pages/mydoc/systemPipeR.pdf){:target="_blank"} ],
[ [.Rmd](https://raw.githubusercontent.com/tgirke/systemPipeR/gh-pages/_vignettes/10_Rworkflows/systemPipeR.Rmd){:target="_blank"} ],
[ [.R](https://raw.githubusercontent.com/tgirke/systemPipeR/gh-pages/_vignettes/10_Rworkflows/systemPipeR.R){:target="_blank"} ],
[ [Slides](https://docs.google.com/presentation/d/175aup31LvnbIJUAvEEoSkpGsKgtBJ2RpQYd0Gs23dLo/embed?start=false&loop=false&delayms=60000){:target="_blank"} ] 


[_`systemPipeR`_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html) provides utilities for building and running automated end-to-end analysis workflows for a wide range of research applications, including next generation sequence (NGS) experiments, such as RNA-Seq, ChIP-Seq, VAR-Seq and Ribo-Seq (H Backman et al., 2016)]. Important features include a uniform workflow interface across different data analysis applications, automated report generation, and support for running both R and command-line software, such as NGS aligners or peak/variant callers, on local computers or compute clusters (Figure 1). The latter supports interactive job submissions and batch submissions to queuing systems of clusters. For instance, _`systemPipeR`_ can be used with most command-line aligners such as `BWA` (Li , 2013; Li et al., 2009), `HISAT2` (Kim et al., 2015), `TopHat2` (Kim et al., 2013) and `Bowtie2` (Langmead et al., 2012), as well as the R-based NGS aligners [_`Rsubread`_](http://www.bioconductor.org/packages/devel/bioc/html/Rsubread.html) (Liao et al., 2013) and [_`gsnap (gmapR)`_](http://www.bioconductor.org/packages/devel/bioc/html/gmapR.html) (Wu et al., 2010). Efficient handling of complex sample sets (_e.g._ FASTQ/BAM files) and experimental designs is facilitated by a well-defined sample annotation infrastructure which improves reproducibility and user-friendliness of many typical analysis workflows in the NGS area (Lawrence et al., 2013). 

The main motivation and advantages of using _`systemPipeR`_ for complex data analysis tasks are:

1. Facilitates design of complex NGS workflows involving multiple R/Bioconductor packages
2. Common workflow interface for different NGS applications
3. Makes NGS analysis with Bioconductor utilities more accessible to new users
4. Simplifies usage of command-line software from within R
5. Reduces complexity of using compute clusters for R and command-line software
6. Accelerates runtime of workflows via parallelzation on computer systems with mutiple CPU cores and/or multiple compute nodes
6. Automates generation of analysis reports to improve reproducibility

<center><img src="./pages/mydoc/systemPipeR_files/utilities.png"></center>

**Figure 1:** Relevant features in _`systemPipeR`_.
Workflow design concepts are illustrated under (A & B). Examples of
*systemPipeR’s* visualization functionalities are given under (C). 

A central concept for designing workflows within the _`sytemPipeR`_ environment 
is the use of workflow management containers. In previous versions, _`sytemPipeR`_ 
used a custom command-line interface called _`SYSargs`_ (see Figure 2) and for 
this purpose will continue to be supported for some time. With the latest [Bioconductor Release 3.9](http://www.bioconductor.org/packages/release/bioc/html/systemPipeR.html), 
we are adopting for this functionality the widely used community standard 
[Common Workflow Language](https://www.commonwl.org/) (CWL) for describing 
analysis workflows in a generic and reproducible manner, introducing _`SYSargs2`_
workflow control class (see Figure 3). Using this community standard in _`sytemPipeR`_
has many advantages. For instance, the integration of CWL allows running _`sytemPipeR`_
workflows from a single specification instance either entirely from within R, from various command-line
wrappers (e.g., *cwl-runner*) or from other languages (*, e.g.,* Bash or Python).
_`sytemPipeR`_ includes support for both command-line and R/Bioconductor software 
as well as resources for containerization, parallel evaluations on computer clusters 
along with the automated generation of interactive analysis reports.

An important feature of _`sytemPipeR's`_ CWL interface is that it provides two
options to run command-line tools and workflows based on CWL. First, one can
run CWL in its native way via an R-based wrapper utility for *cwl-runner* or
*cwl-tools* (CWL-based approach). Second, one can run workflows using CWL's
command-line and workflow instructions from within R (R-based approach). In the
latter case the same CWL workflow definition files (*e.g.* `*.cwl` and `*.yml`)
are used but rendered and exectuted entirely with R functions defined by
_`sytemPipeR`_, and thus use CWL mainly as a command-line and workflow
definition format rather than a software to run workflows. In this regard
_`sytemPipeR`_ also provides several convenience functions that are useful for
designing and debugging workflows, such as a command-line rendering function to
retrieve the exact command-line strings for each data set and processing step
prior to running a command-line.

This tutorial introduces the design of a new CWL S4 class in _`systemPipeR`_, 
as well as the custom command-line interface, combined with the overview of all
the common analysis steps of NGS experiments.


## Workflow design structure using _`SYSargs`_ 

Instances of this S4 object class are constructed by the _`systemArgs`_ function from two simple tabular files: a _`targets`_ file and a _`param`_ file. The latter is optional for workflow steps lacking command-line software. Typically, a _`SYSargs`_ instance stores all sample-level inputs as well as the paths to the corresponding outputs generated by command-line- or R-based software generating sample-level output files, such as read preprocessors (trimmed/filtered FASTQ files), aligners (SAM/BAM files), variant callers (VCF/BCF files) or peak callers (BED/WIG files). Each sample level input/outfile operation uses its own _`SYSargs`_ instance. The outpaths of _`SYSargs`_ usually define the sample inputs for the next _`SYSargs`_ instance. This connectivity is established by writing the outpaths with the _`writeTargetsout`_ function to a new _`targets`_ file that serves as input to the next _`systemArgs`_ call. Typically, the user has to provide only the initial _`targets`_ file. All downstream _`targets`_ files are generated automatically. By chaining several _`SYSargs`_ steps together one can construct complex workflows involving many sample-level input/output file operations with any combinaton of command-line or R-based software. 


<center><img src="./pages/mydoc/systemPipeR_files/SystemPipeR_Workflow.png"></center>

**Figure 2:** Workflow design structure of _`systemPipeR`_ using _`SYSargs`_. 

## Workflow design structure using _`SYSargs2`_ 

The flexibility of _`sytemPipeR's`_ new interface workflow control class is the driving factor behind 
the use of as many steps necessary for the analysis, as well as the connection 
between command-line- or R-based software. The connectivity among all
workflow steps is achieved by the _`SYSargs2`_ workflow control class (see Figure 3).
This S4 class is a list-like container where each instance stores all the
input/output paths and parameter components requried for a particular data
analysis step. _`SYSargs2`_ * instances are generated by two constructor
functions, *loadWorkflow* and *renderWF*, using as data input *targets* or
*yaml* files as well as two *cwl* parameter files (for details see below). When
running preconfigured workflows, the only input the user needs to provide is
the initial *targets* file containing the paths to the input files (*e.g.*
FASTQ) along with unique sample labels. Subsequent targets instances are
created automatically. The parameters required for running command-line
software are provided by the parameter (*.cwl*) files described below. 

We also introduce the *`SYSargs2Pipe`* class that organizes one or many
SYSargs2 containers in a single compound object capturing all information
required to run, control and monitor complex workflows from start to finish. This
design enhances the *`systemPipeR`* workflow framework with a generalized,
flexible, and robust design.

<center><img src="./pages/mydoc/systemPipeR_files/SYSargs2.png"></center>

**Figure 3:** Workflow steps with input/output file operations are controlled by 
_`SYSargs2`_ objects. Each _`SYSargs2`_ instance is constructed from one *targets* 
and two *param* files. The only input provided by the user is the initial *targets* 
file. Subsequent *targets* instances are created automatically, from the previous 
output files. Any number of predefined or custom workflow steps are supported. One
or many _`SYSargs2`_ objects are organized in a *`SYSargs2Pipe`* container.

<br><br><center><a href="mydoc_systemPipeR_1.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeR_2.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
