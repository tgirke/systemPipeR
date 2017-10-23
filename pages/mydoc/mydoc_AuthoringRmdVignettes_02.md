---
title: 2. Getting started
last_updated: Sat Oct 21 17:03:44 2017
sidebar: mydoc_sidebar
permalink: mydoc_AuthoringRmdVignettes_02.html
---

To enable the _Bioconductor_ style in your R Markdown vignette you need to (H Backman et al., 2016):

- Edit the `DESCRIPTION` file by adding

        VignetteBuilder: knitr
        Suggests: BiocStyle, knitr, rmarkdown

- Specify `BiocStyle::html_document` or `BiocStyle::pdf_document` as output format and add vignette metadata in the document header:

        ---
        title: "Vignette Title"
        author: "Vignette Author"
        package: PackageName
        output: 
          BiocStyle::html_document
        vignette: >
          %\VignetteIndexEntry{Vignette Title}
          %\VignetteEngine{knitr::rmarkdown}
          %\VignetteEncoding{UTF-8}  
        ---

The `vignette` section is required in order to instruct _R_ how to build the vignette.^[`\VignetteIndexEntry` should match the `title` of your vignette] The `package` field which should contain the package name is used to print the package version in the output document header. It is not necessary to specify `date` as by default the document compilation date will be automatically included.
See the following section for details on specifying author affiliations and abstract.

BiocStyle's `html_document` and `pdf_document` format functions extend the corresponding original _rmarkdown_ formats, so they accept the same arguments as `html_document` and `pdf_document`, respectively. For example, use `toc_float: true` to obtain a floating TOC as in this vignette.

## Use with R markdown v1

Apart from the default markdown engine implemented in the *[rmarkdown](https://CRAN.R-project.org/package=rmarkdown)* package,
it is also possible to compile _Bioconductor_ documents with the older markdown v1 engine 
from the package *[markdown](https://CRAN.R-project.org/package=markdown)*. There are some
differences in setup and the resulting output between these two
engines. 

To use the *[markdown](https://CRAN.R-project.org/package=markdown)* vignette builder engine:

- Edit the `DESCRIPTION` file to include

        VignetteBuilder: knitr
        Suggests: BiocStyle, knitr
        
- Specify the vignette engine in the `.Rmd` files (inside HTML
  comments)

        <!--
        %% \VignetteEngine{knitr::knitr}
        -->

- Add the following code chunk at the beginning of your `.Rmd`
  vignettes

        ```{r style, echo = FALSE, results = 'asis'}
        BiocStyle::markdown()
        ```

The way of attaching CSS files when using
*[markdown](https://CRAN.R-project.org/package=markdown)* differs from how this is done with *[rmarkdown](https://CRAN.R-project.org/package=rmarkdown)*.
In the former case additional style sheets can be
used by providing them to the `BiocStyle::markdown` function.
To include `custom.css` file use

    ```{r style, echo = FALSE, results = 'asis'}
    BiocStyle::markdown(css.files = c('custom.css'))
    ```

<br><br><center><a href="mydoc_AuthoringRmdVignettes_01.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_AuthoringRmdVignettes_03.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
