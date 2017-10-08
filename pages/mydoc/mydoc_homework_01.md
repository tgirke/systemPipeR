---
title: HW1 - Online Databases and Software
sidebar: mydoc_sidebar
permalink: mydoc_homework_01.html 
---

## Topic: Databases and Software

This is an easy warm-up homework exposing students to a variety of online databases and software tools.

1. Go to [http://www.ncbi.nlm.nih.gov](http://www.ncbi.nlm.nih.gov), select protein DB, run query: `P450 & hydroxylase & human [organism]`, select under _Source_ databases UniProtKB/Swiss-Prot 

    1. Report final query syntax from _Details_ page. 

2. Save GIs from this final query to file (select GI List format under display) 
    1. Report the number of retrieved GIs.

3. Retrieve the corresponding sequences through [Batch-Entrez](http://www.ncbi.nlm.nih.gov/sites/batchentrez) using GI list file as query input -> save sequences in FASTA format

4. Generate multiple alignment and tree of these sequences using [MultAalin](http://bioinfo.genotoul.fr/multalin)

    1. Save multiple alignment and tree to file
    2. Identify putative heme binding cysteine in multiple alignment

5. Open corresponding [UniProt page](http://www.uniprot.org) and search for first P450 sequence in your list.

    1. Compare putative heme binding cysteine with consensus pattern from Prosite database ([Syntax](http://prosite.expasy.org/scanprosite/scanprosite_doc.html#mo_motifs))
	2. Report corresponding Pfam ID
	3. How many mouse (Mus musculus) sequences are in this family (use species tree from Pfam db)

6. [BLASTP](http://www.ncbi.nlm.nih.gov/blast/Blast.cgi) against nr database (use again first P450 in your list); on next two pages click on P450 and CypX domains, respectively (this runs RPS-BLAST). 
    1. Compare resulting alignment with result from MultAlin
	2. View 3D structure in Cn3D, save structure (screen shot) and highlight heme binding cysteine 

## Homework submission

Please assemble the results of this homework in one PDF file and upload it to your private course GitHub repository under `Homework/HW1/HW1.pdf`.

## Due date

Most homeworks will be due one week after they are assigned. This one is due on Thu, April 13th at 6:00 PM.
