---
title: RIBO-Seq Workflow  <br> <br> 1. Introduction
last_updated: Mon Nov 13 16:18:42 2017
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeRIBOseq_01.html
---
Author: Daniela Cassol (danielac@ucr.edu) and Thomas Girke (thomas.girke@ucr.edu)

Last update: 13 November, 2017 

Alternative formats of this tutorial:
[ [HTML](http://girke.bioinformatics.ucr.edu/systemPipeR/pages/mydoc/systemPipeRIBOseq.html){:target="_blank"} ],
[ [PDF](http://girke.bioinformatics.ucr.edu/systemPipeR/pages/mydoc/systemPipeRIBOseq.pdf){:target="_blank"} ],
[ [.Rmd](https://raw.githubusercontent.com/tgirke/systemPipeR/gh-pages/_vignettes/14_RIBOseqWorkflow/systemPipeRIBOseq.Rmd){:target="_blank"} ],
[ [.R](https://raw.githubusercontent.com/tgirke/systemPipeR/gh-pages/_vignettes/14_RIBOseqWorkflow/systemPipeRIBOseq.R){:target="_blank"} ]


Ribo-Seq and polyRibo-Seq are a specific form of RNA-Seq gene expression
experiments utilizing mRNA subpopulations directly bound to ribosomes.
Compared to standard RNA-Seq, their readout of gene expression provides a
better approximation of downstream protein abundance profiles due to their
close association with translational processes. The most important difference
among the two is that polyRibo-Seq utilizes polyribosomal RNA for sequencing,
whereas Ribo-Seq is a footprinting approach restricted to sequencing RNA
fragments protected by ribosomes (Ingolia et al., 2009; Aspden et al., 2014; @Juntawong2015-ru}. 

The workflow presented in this vignette contains most of the data analysis
steps described by (Juntawong et al., 2014) including functionalities useful for
processing both polyRibo-Seq and Ribo-Seq experiments. To improve re-usability
and adapt to recent changes of software versions (_e.g._ R, Bioconductor and
short read aligners), the code has been optimized accordingly. Thus, the
results obtained with the updated workflow are expected to be similar but not
necessarily identical with the published results described in the original
paper. 

Relevant analysis steps of this workflow include read preprocessing, read
alignments against a reference genome, counting of reads overlapping with a
wide range of genomic features (_e.g._ CDSs, UTRs, uORFs, rRNAs, etc.),
differential gene expression and differential ribosome binding analyses, as
well as a variety of genome-wide summary plots for visualizing RNA expression
trends. Functions are provided for evaluating the quality of Ribo-seq data,
for identifying novel expressed regions in the genomes, and for gaining
insights into gene regulation at the post-transcriptional and translational
levels. For example, the functions `genFeatures` and
`featuretypeCounts` can be used to quantify the expression output for
all feature types included in a genome annotation (`e.g.` genes,
introns, exons, miRNAs, intergenic regions, etc.). To determine the approximate
read length of ribosome footprints in Ribo-Seq experiments, these feature type
counts can be obtained and plotted for specific read lengths separately.
Typically, the most abundant read length obtained for translated features
corresponds to the approximate footprint length occupied by the ribosomes of a
given organism group. Based on the results from several Ribo-Seq studies, these
ribosome footprints are typically ~30 nucleotides long
(Ingolia et al., 2011; Ingolia et al., 2009; Juntawong et al., 2014).  However, their
length can vary by several nucleotides depending upon the optimization of the
RNA digestion step and various factors associated with translational
regulation.  For quality control purposes of Ribo-Seq experiments it is also
useful to monitor the abundance of reads mapping to rRNA genes due to the high
rRNA content of ribosomes. This information can be generated with the 
`featuretypeCounts` function described above.

Coverage trends along transcripts summarized for any number of transcripts can
be obtained and plotted with the functions `featureCoverage` and
`plotfeatureCoverage`, respectively. Their results allow monitoring
of the phasing of ribosome movements along triplets of coding sequences.
Commonly, high quality data will display here for the first nucleotide of each
codon the highest depth of coverage computed for the 5' ends of the aligned
reads. 
 
Ribo-seq data can also be used to evaluate various aspects of translational
control due to ribosome occupancy in upstream open reading frames (uORFs). The
latter are frequently present in (or near) 5' UTRs of transcripts. For this,
the function `predORFs` can be used to identify ORFs in the
nucleotide sequences of transcripts or their subcomponents such as UTR regions.
After scaling the resulting ORF coordinates back to the corresponding genome
locations using `scaleRanges`, one can use these novel features
(_e.g._ uORFs) for expression analysis routines similar to those
employed for pre-existing annotations, such as the exonic regions of genes. For
instance, in Ribo-Seq experiments one can use this approach to systematically identify all
transcripts occupied by ribosomes in their uORF regions. The binding of
ribosomes to uORF regions may indicate a regulatory role in the translation of
the downstream main ORFs and/or translation of the uORFs into functionally
relevant peptides. 

## Experimental design
Typically, users want to specify here all information relevant for the analysis
of their NGS study. This includes detailed descriptions of FASTQ files,
experimental design, reference genome, gene annotations, etc.  

<br><br><center><a href="mydoc_systemPipeRIBOseq_01.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeRIBOseq_02.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
