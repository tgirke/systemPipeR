---
title: 9. Cross-references
last_updated: Sat Oct 21 17:03:44 2017
sidebar: mydoc_sidebar
permalink: mydoc_AuthoringRmdVignettes_09.html
---

Apart from referencing figures (Section \@ref(figures)), tables (Section \@ref(tables)), and equations (Section \@ref(equations)), you can also use the same syntax `\@ref(label)` to reference sections, where `label` is the section ID. By default, Pandoc will generate IDs for all section headers, e.g., `# Hello World` will have an ID `hello-world`. In order to avoid forgetting to update the reference label after you change the section header, you may also manually assign an ID to a section header by appending `{#id}` to it.

When a referenced label cannot be found, you will see two question marks like \@ref(non-existing-label), as well as a warning message in the _R_ console when rendering the document.


<br><br><center><a href="mydoc_AuthoringRmdVignettes_08.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_AuthoringRmdVignettes_10.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
