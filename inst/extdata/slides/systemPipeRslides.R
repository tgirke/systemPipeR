## ----install, eval=FALSE-------------------------------------------------
## source("http://bioconductor.org/biocLite.R")
## biocLite("systemPipeR")
## biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE) # From github

## ----load_systemPipeR_hidden, eval=TRUE, include=FALSE, messages=FALSE, warnings=FALSE, cache=TRUE----
library("systemPipeR"); library("systemPipeRdata")

## ----load_systemPipeR_print, eval=FALSE, messages=FALSE, warnings=FALSE, cache=TRUE----
## library("systemPipeR"); library("systemPipeRdata")

## ----show_targets, eval=TRUE, messages=FALSE, warnings=FALSE, cache=TRUE----
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:4,1:5]

## ----generate_workenvir, eval=FALSE, cache=TRUE--------------------------
## ### <b>
## genWorkenvir(workflow="varseq", mydirname=NULL)
## ### </b>
## setwd("varseq")

## ----generate_workenvir_from_shell, eval=FALSE, cache=TRUE, engine="sh"----
## $ echo 'library(systemPipeRdata);
##         genWorkenvir(workflow="varseq", mydirname=NULL)' | R --slave

## ----workflow_template_structure, eval=FALSE-----------------------------
## workflow_name/            # *.Rnw/*.Rmd scripts and targets file
##                 param/    # parameter files for command-line software
##                 data/     # inputs e.g. FASTQ, reference, annotations
##                 results/  # analysis result files

## ----run_make, eval=FALSE, engine="sh"-----------------------------------
## $ make -B

