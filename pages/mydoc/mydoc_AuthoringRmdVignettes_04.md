---
title: 4. Style macros
last_updated: Sat Oct 21 17:03:44 2017
sidebar: mydoc_sidebar
permalink: mydoc_AuthoringRmdVignettes_04.html
---



*[BiocStyle](http://bioconductor.org/packages/BiocStyle)* introduces the following macros useful when
referring to _R_ packages:

* `Biocpkg("IRanges")` for _Bioconductor_ software, annotation and experiment data packages,
  including a link to the release landing page or if the package is only in devel, to the devel landing page, *[IRanges](http://bioconductor.org/packages/IRanges)*.

* `CRANpkg("data.table")` for _R_ packages available on
  CRAN, including a link to the FHCRC CRAN mirror landing page, *[data.table](https://CRAN.R-project.org/package=data.table)*.

* `Githubpkg("rstudio/rmarkdown")` for _R_ packages
  available on GitHub, including a link to the package repository, *[rmarkdown](https://github.com/rstudio/rmarkdown)*.

* `Rpackage("MyPkg")` for _R_ packages that are _not_
  available on _Bioconductor_, CRAN or GitHub; *MyPkg*.

These are meant to be called inline, e.g., `` `r Biocpkg("IRanges")` ``.


<br><br><center><a href="mydoc_AuthoringRmdVignettes_03.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_AuthoringRmdVignettes_05.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
