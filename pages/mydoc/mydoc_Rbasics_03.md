---
title: 3. Installation of R Packages
last_updated: Sun Jun 25 17:41:27 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rbasics_03.html
---

1. Install R for your operating system from [CRAN](http://cran.at.r-project.org/).

2. Install RStudio from [RStudio](http://www.rstudio.com/ide/download).

3. Install CRAN Packages from R console like this:

    
    ```r
    install.packages(c("pkg1", "pkg2")) 
    install.packages("pkg.zip", repos=NULL)
    ```

4. Install Bioconductor packages as follows:

    
    ```r
    source("http://www.bioconductor.org/biocLite.R")
    library(BiocInstaller)
    BiocVersion()
    biocLite()
    biocLite(c("pkg1", "pkg2"))
    ```

5. For more details consult the [Bioc Install page](http://www.bioconductor.org/install/)
and [BiocInstaller](http://www.bioconductor.org/packages/release/bioc/html/BiocInstaller.html) package.

<br><br><center><a href="mydoc_Rbasics_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rbasics_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
