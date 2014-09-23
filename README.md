### systemPipeR: NGS workflow and report generation environment 

[_systemPipeR_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html)
is an R/Bioconductor package for building *end-to-end* analysis pipelines with
automated report generation for next generation sequence (NGS) applications
such as RNA-Seq, ChIP-Seq, BS-Seq, VAR-Seq and many others. An important feature is
support for running command-line software, such as NGS aligners, on both single
machines or compute clusters. This includes both interactive job submissions or
batch submissions to queuing systems of clusters (tested only with Torque).
Efficient handling of complex sample sets and experimental designs is
facilitated by a well-defined sample annotation infrastructure which improves
reproducibility and user-friendliness of many typical analysis workflows in the
NGS area.

#### Installation 
To install the package, please use the _biocLite_ method as instructed 
[_here_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html).


#### Usage
Instructions for running _systemPipeR_ are given in its
[_vignette_](https://github.com/tgirke/systemPipeR/blob/master/vignettes/systemPipeR.pdf?raw=true) (manual).
The sample data set used in the vignette can be downloaded [_here_](http://biocluster.ucr.edu/~tgirke/projects/systemPipeR_test_data.zip). 
The expected format to define NGS samples (_e.g._ FASTQ files) and their
labels are given in
[_targets.txt_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targets.txt)
and
[_targetsPE.txt_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targetsPE.txt)
(latter is for PE reads). 
The run parameters of command-line software are defined by param files that
have a simplified JSON-like name/value structure. Here is a sample param file
for _Tophat2_:
[_tophat.param_](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/tophat.param).
Templates for setting up custom project reports are provided as _*.Rnw_
files in the [vignettes](https://github.com/tgirke/systemPipeR/tree/master/vignettes) subdirectory of this package. The
corresponding PDFs of these report templates are linked here:
[systemPipeRNAseq](https://github.com/tgirke/systemPipeR/blob/master/vignettes/systemPipeRNAseq.pdf?raw=true),
[systemPipeChIPseq](https://github.com/tgirke/systemPipeR/blob/master/vignettes/systemPipeChIPseq.pdf?raw=true)
and
[systemPipeVARseq](https://github.com/tgirke/systemPipeR/blob/master/vignettes/systemPipeVARseq.pdf?raw=true).

