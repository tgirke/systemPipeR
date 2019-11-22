---
title: 10. Venn diagram of variants
last_updated: Thu Nov 21 15:58:58 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_10.html
---

The venn diagram utilities defined by the `systemPipeR` package can be used to
identify common and unique variants reported for different samples
and/or variant callers. The below generates a 4-way venn diagram
comparing four sampes for each of the two variant callers.


```r
dir_path <- system.file("extdata/cwl/varseq", package = "systemPipeR")
## gatk
args <- loadWorkflow(targets = "results/targets_report_gatk.txt", 
    wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
varlist <- sapply(names(subsetWF(args[1:2], slot = "output", 
    subset = 1, index = 1)), function(x) as.character(read.delim(subsetWF(args[1:2], 
    slot = "output", subset = 1, index = 1)[x])$VARID))
vennset_gatk <- overLapper(varlist, type = "vennsets")

## bcf
args <- loadWorkflow(targets = "./results/targets_report_bcf.txt", 
    wf_file = "combine.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_", SampleName = "_SampleName_"))
varlist <- sapply(names(subsetWF(args[1:2], slot = "output", 
    subset = 1, index = 1)), function(x) as.character(read.delim(subsetWF(args[1:2], 
    slot = "output", subset = 1, index = 1)[x])$VARID))
vennset_bcf <- overLapper(varlist, type = "vennsets")

pdf("./results/vennplot_var.pdf")
vennPlot(list(vennset_gatk, vennset_bcf), mymain = "", mysub = "GATK: red; BCFtools: blue", 
    colmode = 2, ccol = c("red", "blue"))
dev.off()
```

![](./pages/mydoc/systemPipeVARseq_files/vennplot_var.png)
<div align="center">Figure 2: Venn Diagram for 4 samples from GATK and BCFtools</div>


<br><br><center><a href="mydoc_systemPipeVARseq_09.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_11.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
