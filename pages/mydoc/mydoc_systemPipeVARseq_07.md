---
title: 7. Annotate filtered variants
last_updated: Sat Apr 18 12:30:50 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_07.html
---

The function `variantReport` generates a variant report using
utilities provided by the `VariantAnnotation` package. The report for
each sample is written to a tabular file containing genomic context annotations
(_e.g._ coding or non-coding SNPs, amino acid changes, IDs of affected
genes, etc.) along with confidence statistics for each variant. The CWL
file `param/cwl/varseq_downstream/annotate.cwl` defines the paths to the input 
and output files which are stored in a `SYSargs2` instance. 

## Basics of annotating variants

Variants overlapping with common annotation features can be identified with `locateVariants`.

```r
dir_path <- system.file("extdata/cwl/varseq", package = "systemPipeR")
library("GenomicFeatures")
args <- loadWorkflow(targets = "./results/targets_filter_gatk.txt", 
    wf_file = "annotate.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
txdb <- loadDb("./data/tair10.sqlite")
vcf <- readVcf(infile1(args)[1], "A. thaliana")
locateVariants(vcf, txdb, CodingVariants())
```

Synonymous/non-synonymous variants of coding sequences are computed by the predictCoding function for variants overlapping with coding regions.


```r
fa <- FaFile(normalizePath(file.path(args$yamlinput$data_path$path, 
    args$yamlinput$ref_name)))
predictCoding(vcf, txdb, seqSource = fa)
```

## Annotate filtered variants `GATK`


```r
library("GenomicFeatures")
args <- loadWorkflow(targets = "./results/targets_filter_gatk.txt", 
    wf_file = "annotate.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(normalizePath(file.path(args$yamlinput$data_path$path, 
    args$yamlinput$ref_name)))
suppressAll(variantReport(args = args, txdb = txdb, fa = fa, 
    organism = "A. thaliana"))
writeTargetsout(x = args, file = "./results/targets_report_gatk.txt", 
    step = 1, new_col = "FileName1", new_col_output_index = 1, 
    overwrite = TRUE)
```

## Annotate filtered variants `bcftools`


```r
library("GenomicFeatures")
args <- loadWorkflow(targets = "./results/targets_filter_bcf.txt", 
    wf_file = "annotate.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
txdb <- loadDb("./data/tair10.sqlite")
fa <- FaFile(normalizePath(file.path(args$yamlinput$data_path$path, 
    args$yamlinput$ref_name)))
suppressAll(variantReport(args = args, txdb = txdb, fa = fa, 
    organism = "A. thaliana"))
writeTargetsout(x = args, file = "./results/targets_report_bcf.txt", 
    step = 1, new_col = "FileName1", new_col_output_index = 1, 
    overwrite = TRUE)
```

View annotation result for single sample

```r
read.delim(output(args)[[1]][[1]])[38:40, ]
```

<br><br><center><a href="mydoc_systemPipeVARseq_06.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_08.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
