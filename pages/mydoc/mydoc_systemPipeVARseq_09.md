---
title: 9. Summary statistics of variants
last_updated: Sat Apr 18 12:30:50 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_09.html
---

The `varSummary` function counts the number of variants for each feature type
included in the anntation reports.

## Summary of variants `GATK`


```r
args <- loadWorkflow(targets = "./results/targets_report_gatk.txt", 
    wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
varSummary(args)
write.table(varSummary(args), "./results/variantStats_gatk.xls", 
    quote = FALSE, col.names = NA, sep = "\t")
```

## Summary of variants `bcf`


```r
args <- loadWorkflow(targets = "./results/targets_report_bcf.txt", 
    wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
varSummary(args)
write.table(varSummary(args), "./results/variantStats_bcf.xls", 
    quote = FALSE, col.names = NA, sep = "\t")
```

<br><br><center><a href="mydoc_systemPipeVARseq_08.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_10.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
