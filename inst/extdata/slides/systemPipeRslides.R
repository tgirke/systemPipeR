## ----install, eval=FALSE-------------------------------------------------
## source("http://bioconductor.org/biocLite.R")
## biocLite("systemPipeR")
## biocLite("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE) # From github

## ----load_systemPipeR_hidden, eval=TRUE, include=FALSE, messages=FALSE, warnings=FALSE, cache=TRUE----
library("systemPipeR"); library("systemPipeRdata")

## ----load_systemPipeR_print, eval=FALSE, messages=FALSE, warnings=FALSE, cache=TRUE----
## library("systemPipeR")
## library("systemPipeRdata")

## ----accessing_help, eval=FALSE, messages=FALSE, warnings=FALSE, cache=TRUE----
## library(help="systemPipeR")
## vignette("systemPipeR")

## ----show_targetsSE, eval=TRUE, messages=FALSE, warnings=FALSE, cache=TRUE----
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:3,1:5]

## ----show_targetsPE, eval=TRUE, messages=FALSE, warnings=FALSE, cache=TRUE----
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:3,1:4]

## ----sysargs_instance, eval=TRUE, messages=FALSE, warnings=FALSE, cache=TRUE----
parampath <- system.file("extdata", "tophat.param", package="systemPipeR")
args <- suppressWarnings(systemArgs(sysma=parampath, mytargets=targetspath))
args

## ----sysargs_names, eval=TRUE, messages=FALSE, warnings=FALSE, cache=TRUE----
names(args)[c(5,8,13)]

## ----sysargs_args, eval=FALSE, messages=FALSE, warnings=FALSE, cache=TRUE----
## sysargs(args)[1]

## ----run_args_single, eval=FALSE, messages=FALSE, warnings=FALSE, cache=TRUE----
## runCommandline(args)

## ----run_args_cluster, eval=FALSE, messages=FALSE, warnings=FALSE, cache=TRUE----
## clusterRun(args, ...)

## ----generate_workenvir, eval=FALSE, cache=TRUE--------------------------
## ### <b>
## genWorkenvir(workflow="varseq", mydirname=NULL)
## ### </b>
## setwd("varseq")

## ----generate_workenvir_from_shell, eval=FALSE, cache=TRUE, engine="sh"----
## $ echo 'library(systemPipeRdata);
##         genWorkenvir(workflow="varseq", mydirname=NULL)' | R --slave

## ----workflow_template_structure, eval=FALSE-----------------------------
## ### <b>
## workflow_name/            # *.Rnw/*.Rmd scripts, targets file, etc.
##                 param/    # parameter files for command-line software
##                 data/     # inputs e.g. FASTQ, reference, annotations
##                 results/  # analysis result files
## ### </b>

## ----run_make, eval=FALSE, engine="sh"-----------------------------------
## $ make -B

