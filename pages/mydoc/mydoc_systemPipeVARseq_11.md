---
title: 11. Plot variants programmatically 
last_updated: Thu Nov 21 15:58:58 2019
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_11.html
---

The following plots a selected variant with `ggbio`.

In this example, the input `BAM` file is from the `GATK` step 5, analysis ready bam. 
You can use other aligned `BAMs` as well, but make sure they are indexed. The `VCF` 
file is taken from `Inspect VCF file ` section or you can load your own vcf.


```r
library(ggbio)
mychr <- "ChrC"
mystart <- 11000
myend <- 13000
args <- loadWorkflow(targets = "results/targets_gatk.txt", wf_file = "combine.cwl", 
    input_file = "varseq.yml", dir_path = "param/cwl/varseq_downstream/")
args <- renderWF(args, inputvars = c(GATK_FIXED = "_FILE1_", 
    SampleName = "_SampleName_"))
ga <- readGAlignments(subsetWF(args, slot = "input", subset = 1)[1], 
    use.names = TRUE, param = ScanBamParam(which = GRanges(mychr, 
        IRanges(mystart, myend))))
p1 <- autoplot(ga, geom = "rect")
p2 <- autoplot(ga, geom = "line", stat = "coverage")
p3 <- autoplot(vcf[seqnames(vcf) == mychr], type = "fixed") + 
    xlim(mystart, myend) + theme(legend.position = "none", axis.text.y = element_blank(), 
    axis.ticks.y = element_blank())
p4 <- autoplot(loadDb("./data/tair10.sqlite"), which = GRanges(mychr, 
    IRanges(mystart, myend)), names.expr = "gene_id")
png("./results/plot_variant.png")
tracks(Reads = p1, Coverage = p2, Variant = p3, Transcripts = p4, 
    heights = c(0.3, 0.2, 0.1, 0.35)) + ylab("")
dev.off()
```

![](./pages/mydoc/systemPipeVARseq_files/plot_variant.png)
<div align="center">Figure 3: Plot variants with programmatically.</div>

<br><br><center><a href="mydoc_systemPipeVARseq_10.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_12.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
