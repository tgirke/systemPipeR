---
title: HW4 - Pairwise Alignments
sidebar: mydoc_sidebar
permalink: mydoc_homework_04.html 
---

## A. Choice of Sequence Type

- __Task 1__: Which sequence type - amino acid or nucleotide - is more appropriate to search databases for remotely related sequences? Provide at least three reasons for your decision. 
    
## B. Dynamic Programming for Pairwise Alignments
- __Task 2__: Create manually (or write an R script for it) one global and one local alignment for the following two protein sequences using the Needleman-Wusch and Smith-Waterman algorithms, respectively:
    
{% highlight sh %}
O15528: PFGFGKRSCMGRRLA
P98187: FIPFSAGPRNCIGQK
{% endhighlight %}

Use in each case BLOSUM50 as substitution matrix and 8 as gap opening and extension penalties. 
Note, [here](https://github.com/tgirke/GEN242/blob/gh-pages/_vignettes/06_Homework/mydoc_homework_04.R) is some R code to create the initial matrix programmatically for upload to a spreadsheet program. Alternatively, solve the entire homework by writing an R script.
Your answers should contain the following components: 

1. Manually populated dynamic programming matrices
2. The optimal pairwise alignments created by traceback 
3. The final scores of the alignments

	
## C. Alignments with Different Substitution Matrices

- __Task 1__: Load the `Biostrings` package in R, import the following two cytochrome P450 sequences `O15528` and `P98187` from [NCBI](http://www.ncbi.nlm.nih.gov/protein/O15528,P98187) (save as `myseq.fasta`), and create a global alignment with the `pairwiseAlignment` function from `Biostrings` as follows:

{% highlight r %}
library(Biostrings)
myseq <- readAAStringSet("myseq.fasta", "fasta")
(p <- pairwiseAlignment(myseq[[1]], myseq[[2]], type="global", substitutionMatrix="BLOSUM50"))
writePairwiseAlignments(p)
{% endhighlight %}

Your answers should address the following items: 
		
1. Record the scores for the scoring matrices BLOSUM50, BLOSUM62 and BLOSUM80.
2. How and why do the scores differ for the three scoring matrices?

## Homework submission

Assemble the results from this homework in one PDF file (`HW4.pdf`) and upload it to your private GitHub repository under `Homework/HW4/HW4.pdf`.

## Due date

This homework is due in two weeks on Tue, April 25th at 6:00 PM.
