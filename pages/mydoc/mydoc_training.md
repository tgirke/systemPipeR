---
title: Training
sidebar: mydoc_sidebar
permalink: mydoc_training.html
---

## BioC2019

+ [Overview Slide Show - 2019](./Bioc2019.html)
+ [Workshop Vignette](./)

## Slide Shows

+ [Overview Slide Show - 2016](https://docs.google.com/presentation/d/175aup31LvnbIJUAvEEoSkpGsKgtBJ2RpQYd0Gs23dLo/embed?start=false&loop=false&delayms=60000)
+ [Overview Slide Show - 2015](./systemPipeRslides.html)


## Material for Bioc2016 Tutorial 
+ [Intro Slide Show](https://docs.google.com/presentation/d/175aup31LvnbIJUAvEEoSkpGsKgtBJ2RpQYd0Gs23dLo/embed?start=false&loop=false&delayms=60000)
+ Tutorial Material 
    + Introduction to `systemPipeR` 
        + [HTML](https://htmlpreview.github.io/?https://raw.githubusercontent.com/tgirke/systemPipeRdata/master_github_deprecated/vignettes/systemPipeR_Intro.html)
        + [PDF](https://raw.githubusercontent.com/tgirke/systemPipeRdata/master_github_deprecated/vignettes/systemPipeR_Intro.pdf)
        + [Rmd](https://raw.githubusercontent.com/tgirke/systemPipeRdata/master_github_deprecated/vignettes/systemPipeR_Intro.Rmd)
    + Demo Workflow: RIBO-Seq 
        + [HTML](https://htmlpreview.github.io/?https://raw.githubusercontent.com/tgirke/systemPipeRdata/master_github_deprecated/inst/extdata/workflows/riboseq/systemPipeRIBOseq.html)
        + [PDF](https://raw.githubusercontent.com/tgirke/systemPipeRdata/master_github_deprecated/inst/extdata/workflows/riboseq/systemPipeRIBOseq.pdf)
        + [Rmd](https://raw.githubusercontent.com/tgirke/systemPipeRdata/master_github_deprecated/inst/extdata/workflows/riboseq/systemPipeRIBOseq.Rmd)

### Installation 
Please install _systemPipeRdata_ from this GitHub repository as shown below. This package provides the data and `Rmd` files for the tutorial. 
Its parent package _systemPipeR_ is a dependency and it should install along with all its own dependencies automatically. If it doesn't then please also install the package, using the BiocManager::install command given below.

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("tgirke/systemPipeRdata", build_vignettes=TRUE, dependencies=TRUE)
```
