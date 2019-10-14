### systemPipeR: NGS workflow and report generation environment 

[![platforms](http://www.bioconductor.org/shields/availability/3.10/systemPipeR.svg)](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html#archives)
[![rank](http://www.bioconductor.org/shields/downloads/devel/systemPipeR.svg)](http://bioconductor.org/packages/stats/bioc/systemPipeR/)
[![posts](http://www.bioconductor.org/shields/posts/systemPipeR.svg)](https://support.bioconductor.org/t/systempiper/)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/systemPipeR.svg)](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html#since)
[![build](http://www.bioconductor.org/shields/build/devel/bioc/systemPipeR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/systemPipeR/)
[![updated](http://www.bioconductor.org/shields/lastcommit/devel/bioc/systemPipeR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/systemPipeR/)
[![dependencies](http://www.bioconductor.org/shields/dependencies/devel/systemPipeR.svg)](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html#since)

### Introduction

[_systemPipeR_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html)
is an R/Bioconductor package for building and running automated *end-to-end*
analysis workflows for a wide range of research applications, including next-generation 
sequencing (NGS) experiments, such as RNA-Seq, ChIP-Seq, VAR-Seq and Ribo-Seq.
Important features include a uniform workflow interface across different data analysis 
applications, automated report generation, and support for running both R and command-line software,
such as NGS aligners or peak/variant callers, on local computers or compute
clusters. The latter supports interactive job submissions and batch submissions
to queuing systems of clusters. Efficient handling of complex sample sets and
experimental designs is facilitated by a well-defined sample annotation
infrastructure which improves reproducibility and user-friendliness of many
typical analysis workflows in the NGS area.

#### Installation 
To install the package, please use the _`BiocManager::install`_ command:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("systemPipeR")
```

To obtain the most recent updates immediately, one can install it directly from
github as follow:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("tgirke/systemPipeR", build_vignettes=TRUE, dependencies=TRUE)
```

#### Usage

Instructions for running _systemPipeR_ are given in its main
[_vignette_](http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html) (manual).
The sample data set used in the vignette are provided by the data package [_systemPipeRdata_](http://www.bioconductor.org/packages/devel/data/experiment/html/systemPipeRdata.html).
The expected format to define NGS samples (_e.g._ FASTQ files) and their
labels are given in
[_targets.txt_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targets.txt)
and
[_targetsPE.txt_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targetsPE.txt)
(latter is for PE reads).
With the latest [Bioconductor Release 3.9](http://www.bioconductor.org/packages/release/bioc/html/systemPipeR.html), 
we are adopting for this functionality the widely used community standard 
[Common Workflow Language](https://www.commonwl.org/) (CWL) for describing 
analysis workflows in a generic and reproducible manner, introducing _`SYSargs2`_
workflow control class. Using this community standard in _`systemPipeR`_
has many advantages. For instance, the integration of CWL allows running _`sytemPipeR`_
workflows from a single specification instance either entirely from within R, from various command-line
wrappers (e.g., *cwl-runner*) or from other languages (*, e.g.,* Bash or Python).
The run parameters of command-line software are defined by param files that
have a simplified YAML name/value structure. Here is a sample param file
for _Hisat2_:
[_hisat2.cwl_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/cwl/hisat2/hisat2-pe/hisat2-mapping-pe.cwl).
Templates for setting up custom project reports are provided by [_systemPipeRdata_](https://github.com/tgirke/systemPipeRdata).
The corresponding PDFs of these report templates are linked here:
[systemPipeRNAseq](http://www.bioconductor.org/packages/devel/data/experiment/vignettes/systemPipeRdata/inst/doc/systemPipeRNAseq.html),
[systemPipeRIBOseq](http://www.bioconductor.org/packages/devel/data/experiment/vignettes/systemPipeRdata/inst/doc/systemPipeRIBOseq.html),
[systemPipeChIPseq](http://www.bioconductor.org/packages/devel/data/experiment/vignettes/systemPipeRdata/inst/doc/systemPipeChIPseq.html)
and
[systemPipeVARseq](http://www.bioconductor.org/packages/devel/data/experiment/vignettes/systemPipeRdata/inst/doc/systemPipeVARseq.html).

#### Slides

+ [Overview Slide Show](http://girke.bioinformatics.ucr.edu/systemPipeR/pages/mydoc/systemPipeRslides.html).