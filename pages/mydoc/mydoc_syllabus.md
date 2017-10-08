---
layout: page
title: Syllabus - GEN242
sidebar: mydoc_sidebar
permalink: mydoc_syllabus.html 
---

## Course title

Data Analysis in Genome Biology <br>
GEN242 - Spring 2017

## Instructor
Name: Thomas Girke <br> 
Email: thomas.girke@ucr.edu <br>
Office location: 1207F Genomics <br>
Office hours: Tue & Thu 4:30 - 5:30 PM 


## Description

Introduction to algorithms, statistical methods and data analysis programming
routines relevant for genome biology. The class consists of three
main components: lectures, hands-on practicals and student course projects. The
lecture topics cover databases, sequence (NGS) analysis, phylogenetics,
comparative genomics, genome-wide profiling methods, network biology and more.
The hands-on practicals include homework assignments and course projects
focusing on data analysis programming of next generation genome data using
command-line tools on a computer cluster and the programming environment R.
Credit: 4 units (4 hours lecture, 2 hours discussion)

## Objectives of course

- Acquire understanding of algorithms used in bioinformatics
- Obtain hands-on experience in large scale data analysis.

## Prerequisites
The main prerequisite for this course is a strong interest in acquiring the
skills required for mastering the computational aspects of modern genome
research.

## Structure of course
Two lectures per week (2 hours each) plus one discussion section (1 hour).
During the first weeks the discussion section will be used for data analysis
tutorials using Linux command-line tools and R. 

## Time
Lecture: Tue/Thu 2:10-3:30 PM <br>
Discussion: Thu 3:40-4:30 PM

## Location
[2130 CHASS Interdisciplinary Bldg-South (INST)](http://classrooms.cmsdev.ucr.edu/buildings/chassINTS/2130.html)

## Grading
1. Homework assignments: 40%
2. Scientific paper presentation: 20%
3. Course project presentations: 20%
4. Final project report: 20%

## Materials needed
Students are expected to bring to each class meeting a laptop with a functional wireless
connection and a recent internet browser version (e.g. Firefox, Chrome or
Safari) preinstalled. Tablet computers with mobile operating systems are not
suitable for running the required software. User accounts on a research
computer cluster will be provided at the beginning of the course. To log in to
the cluster, students also need to install a terminal application for their
operating system (_e.g._ [iTerm2](http://www.iterm2.com/#/section/home) on OS X,
and [PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) or 
[MobaXterm](http://mobaxterm.mobatek.net/) on Windows) as well as a file exchange software such as
[FileZilla](https://filezilla-project.org/download.php?type=client). In
addition, a recent version of [R](http://www.r-project.org) and
[RStudio](http://rstudio.org) should be installed on each laptop.

## Schedule

|Week     |Topic                                                    |
|---------|---------------------------------------------------------|
| Week 1  | Course Introduction                                     |
|         | Databases and Software for Genome Biology               |
|         | Tutorial: Introduction to Linux and Biocluster          |
| Week 2  | Sequencing Technologies                                 |
|         | Tutorial: Introduction to R                             | 
| Week 3  | Sequence Alignments and Searching                       |
|         | Tutorial: Programming in R                              |
| Week 4  | Multiple Sequence Alignments                            |
|         | Short Read Alignment Algorithms                         |
|         | Tutorial: Basics of NGS Analysis                        |
| Week 5  | Gene Expression Analysis using Microarrays and RNA-Seqs |
|         | Tutorial: NGS Workflow Overview                         |
|         | Tutorial: RNA-Seq Analysis                              |
| Week 6  | Analysis of ChIP-Seq and VAR-Seq Experiments            |
|         | Tutorial: ChIP-Seq Analysis                             |
|         | Tutorial: VAR-Seq Analysis                              |
| Week 7  | Student Paper Presentations                             |
| Week 8  | Clustering algorithms                                   |
|         | Annotation Systems and Gene Set Enrichment Analysis     |
|         | Tutorial: Gene Set Enrichment Analysis                  |
| Week 9  | Profile HMMs for Protein Family Modeling                |
|         | Phylogenetics                                           |
|         | Tutorial: Graphics and Data Visualization               |
| Week 10 | Student Project Presentations                           |
|         | Final Course Discussion                                 |
|---------|---------------------------------------------------------|

## Reading list

### Journal articles

Alkan, C, Sajjadian, S, Eichler, E E (2011) Limitations of next-generation genome sequence assembly. Nat Methods, 8: 61-65. http://www.hubmed.org/display.cgi?uids=21102452

Anders, S, Reyes, A, Huber, W (2012) Detecting differential usage of exons from RNA-seq data. Genome Res, 22: 2008-2017. http://www.hubmed.org/display.cgi?uids=22722343

DePristo, M A, Banks, E, Poplin, R, Garimella, K V, Maguire, J R, Hartl, C, Philippakis, A A, del Angel, G, Rivas, M A, Hanna, M, McKenna, A, Fennell, T J, Kernytsky, A M, Sivachenko, A Y, Cibulskis, K, Gabriel, S B, Altshuler, D, Daly, M J (2011) A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nat Genet, 43: 491-498.
http://www.hubmed.org/display.cgi?uids=21478889

Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., Gingeras, T.R., 2012. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21.
http://www.hubmed.org/display.cgi?uids=23104886

Grabherr, M G, Haas, B J, Yassour, M, Levin, J Z, Thompson, D A, Amit, I, Adiconis, X, Fan, L, Raychowdhury, R, Zeng, Q, Chen, Z, Mauceli, E, Hacohen, N, Gnirke, A, Rhind, N, di Palma, F, Birren, B W, Nusbaum, C, Lindblad-Toh, K, Friedman, N, Regev, A (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nat Biotechnol, 29: 644-652.
http://www.hubmed.org/display.cgi?uids=21572440

Langmead, B, Salzberg, S L (2012) Fast gapped-read alignment with Bowtie 2. Nat Methods, 9: 357-359. http://www.hubmed.org/display.cgi?uids=22388286

Landt et al. (2012) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Res, 22: 1813-1831. http://www.hubmed.org/display.cgi?uids=22955991

Li, H, Durbin, R (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25: 1754-1760. http://www.hubmed.org/display.cgi?uids=19451168

Machanick, P, Bailey, T L (2011) MEME-ChIP: motif analysis of large DNA datasets. Bioinformatics, 27: 1696-1697.http://www.hubmed.org/display.cgi?uids=21486936

Soneson, C, Delorenzi, M (2013) A comparison of methods for differential expression analysis of RNA-seq data. BMC Bioinformatics, 14: 91-91. http://www.hubmed.org/display.cgi?uids=23497356

Tompa, M, Li, N, Bailey, T L, Church, G M, De Moor, B, Eskin, E, Favorov, A V, Frith, M C, Fu, Y, Kent, W J, Makeev, V J, Mironov, A A, Noble, W S, Pavesi, G, Pesole, G, R{\'e}gnier, M, Simonis, N, Sinha, S, Thijs, G, van Helden, J, Vandenbogaert, M, Weng, Z, Workman, C, Ye, C, Zhu, Z (2005) Assessing computational tools for the discovery of transcription factor binding sites. Nat Biotechnol, 23: 137-144.
http://www.hubmed.org/display.cgi?uids=15637633

Trapnell, C, Hendrickson, D G, Sauvageau, M, Goff, L, Rinn, J L, Pachter, L (2013) Differential analysis of gene regulation at transcript resolution with RNA-seq. Nat Biotechnol, 31: 46-53. 
http://www.hubmed.org/display.cgi?uids=23222703

Wilbanks, E G, Facciotti, M T (2010) Evaluation of algorithm performance in ChIP-seq peak detection. PLoS One, 5. http://www.hubmed.org/display.cgi?uids=20628599

Zeitouni, B, Boeva, V, Janoueix-Lerosey, I, Loeillet, S, Legoix-n{\'e}, P, Nicolas, A, Delattre, O, Barillot, E (2010) SVDetect: a tool to identify genomic structural variations from paired-end and mate-pair sequencing data. Bioinformatics, 26: 1895-1896. http://www.hubmed.org/display.cgi?uids=20639544


### Books 

Note: there is no need to purchase any books for this course as most reading material will be based on journal articles!

General
Jonathan Pevsner (2009) Bioinformatics and Functional Genomics. Wiley-Blackwell; 2nd Edition, 992 pages.

Algorithms
Jones N and Pevzner P (2004) An Introduction to Bioinformatics Algorithms. MIT Press, Massachusetts, 435 pages.

Sequence Analysis
Durbin, R, Eddy, S, Krogh, A, Mitchison, G. (1998) Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids. Cambridge University Press, UK, 356 pages.

Parida L (2008) Pattern Discovery in Bioinformatics: Theory & Algorithms. CRC Press, London, 526 pages.

Profiling Bioinformatics
Gentleman, R, Carey, V, Dudoit, S, Irizarry, R, Huber, W (2005) Bioinformatics and Computational Biology Solutions Using R and Bioconductor. Springer, New York, 473 pages.

Phylogenetics 
Felsenstein, J (2004) Inferring Phylogenies. Sinauer, Massachusetts, 664 pages.

Paradis (2006) Analysis of Phylogenetics and Evolution with R. Springer, New York, 211 pages.

