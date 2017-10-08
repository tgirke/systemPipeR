---
title: RNA-Seq1 - RNA-Seq Aligners 
sidebar: mydoc_sidebar
permalink: mydoc_project_rnaseq1.html 
---

## RNA-Seq Workflow  

1. Read quality assessment, filtering and trimming 
2. Map reads against reference genome 
3. Perform read counting for required ranges (_e.g._ exonic gene ranges)
4. Normalization of read counts
5. Identification of differentially expressed genes (DEGs)
6. Clustering of gene expression profiles 
7. Gene set enrichment analysis

## Challenge Project: Comparison of RNA-Seq Aligners 

+ Run workflow from start to finish (steps 1-7) on RNA-Seq data set from Howard et al. (2013)
+ Challenge project tasks
    + Run at least 2-3 RNA-Seq alignment tools such as Bowtie2/Tophat2, HISAT and Kallisto, and evaluate the impact of the aligner on:
        + Read counts
        + Differentially expressed genes (DEGs)
        + Generate plots to compare the results efficiently

## References

+ Howard, B.E. et al., 2013. High-throughput RNA sequencing of pseudomonas-infected Arabidopsis reveals hidden transcriptome complexity and novel splice variants. PloS one, 8(10), p.e74183. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24098335)
+ Kim D, Pertea G, Trapnell C, Pimentel H, Kelley R, Salzberg SL (2013) TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. Genome Biol. doi: 10.1186/gb-2013-14-4-r36 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/23618408)
+ Kim D, Langmead B, Salzberg SL (2015) HISAT: a fast spliced aligner with low memory requirements. Nat Methods 12: 357–360 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/25751142)
+ Liao Y, Smyth GK, Shi W (2013) The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Res 41: e108 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/23558742)
+ Bray NL, Pimentel H, Melsted P, Pachter L (2016) Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol. doi: 10.1038/nbt.3519 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/27043002)


