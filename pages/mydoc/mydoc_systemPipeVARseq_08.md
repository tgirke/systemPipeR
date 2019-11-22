---
title: 8. Combine annotation results among samples
last_updated: Thu Nov 21 15:58:58 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_08.html
---

To simplify comparisons among samples, the `combineVarReports`
function combines all variant annotation reports referenced in a
`SYSargs2` instance (here `args`). At the same time the function
allows to consider only certain feature types of interest. For instance, the
below setting `filtercol=c(Consequence="nonsynonymous")` will include
only nonsysynonymous variances listed in the `Consequence` column of
the annotation reports. To omit filtering, one can use the setting
`filtercol="All"`.

## Combine results `GATK`


```r
dir_path <- system.file("extdata/cwl/varseq", package = "systemPipeR")
args <- loadWorkflow(targets = "results/targets_report_gatk.txt", 
    wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
combineDF <- combineVarReports(args, filtercol = c(Consequence = "nonsynonymous"))
write.table(combineDF, "./results/combineDF_nonsyn_gatk.xls", 
    quote = FALSE, row.names = FALSE, sep = "\t")
```

## Combine results  `bcftools`


```r
args <- loadWorkflow(targets = "results/targets_report_bcf.txt", 
    wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
combineDF <- combineVarReports(args, filtercol = c(Consequence = "nonsynonymous"))
write.table(combineDF, "./results/combineDF_nonsyn_bcf.xls", 
    quote = FALSE, row.names = FALSE, sep = "\t")
```

<br><br><center><a href="mydoc_systemPipeVARseq_07.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_09.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
