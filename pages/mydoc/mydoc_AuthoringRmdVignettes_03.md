---
title: 3. Document header
last_updated: Sat Oct 21 17:03:44 2017
sidebar: mydoc_sidebar
permalink: mydoc_AuthoringRmdVignettes_03.html
---

## Author affiliations

The `author` field allows for specifing author names along with affiliation and email information.

In the basic case, when no additional information apart from author names is provided, these can be entered as a single character string

    author: "Single Author"
    
or a list

    author:
    - First Author
    - Second Author
    - Last Author
    
which will print as "First Author, Second Author and Last Author".

Author affiliations and emails can be entered in named sublists of the author list. Multiple affiliations per author can be specified this way.

    author:
    - name: First Author
      affiliation: 
      - Shared affiliation
      - Additional affiliation
    - name: Second Author
      affiliation: Shared affiliation
      email: corresponding@author.com

A list of unique affiliations will be displayed below the authors, similar as in this document.

For clarity, compactness, and to avoid errors,
repeated nodes in YAML header can be initially denoted by an anchor entered with an ampersand &,
and later referenced with an asterisk *. For example, the above affiliation metadata is equivalent to the shorthand notation

    author:
    - name: First Author
      affiliation: 
      - &id Shared affiliation
      - Additional affiliation
    - name: Second Author
      affiliation: *id
      email: corresponding@author.com


## Abstract and running headers

Abstract can be entered in the corresponding field of the document front matter,
as in the example below.

    ---
    title: "Full title for title page"
    shorttitle: "Short title for headers"
    author: "Vignette Author"
    package: PackageName
    abstract: >
      Document summary
    output: 
      BiocStyle::pdf_document
    ---

The `shorttitle` option specifies the title used in running headers 
instead of the document title.^[only relevant to PDF output]


<br><br><center><a href="mydoc_AuthoringRmdVignettes_02.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_AuthoringRmdVignettes_04.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
