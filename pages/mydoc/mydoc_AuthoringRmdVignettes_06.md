---
title: 6. Figures
last_updated: Sat Oct 21 17:03:44 2017
sidebar: mydoc_sidebar
permalink: mydoc_AuthoringRmdVignettes_06.html
---

_BiocStyle_ comes with three predefined figure sizes. Regular figures not otherwise specified appear indented with respect to the paragraph text, as in the example below.


```r
plot(cars)
```

<img src="./pages/mydoc/AuthoringRmdVignettes_files/no-cap-1.png" width="100%" />

Figures which have no captions are just placed wherever they were generated in the _R_ code.
If you assign a caption to a figure via the code chunk option `fig.cap`, the plot will be automatically labeled and numbered^[for PDF output it will be placed in a floating figure environment], and it will be also possible to reference it. These features are provided by *[bookdown](https://CRAN.R-project.org/package=bookdown)*, which defines a format-independent syntax for specifying cross-references, see Section \@ref(cross-references). The figure label is generated from the code chunk label^[for cross-references to work the chunk label may only contain alphanumeric characters (a-z, A-Z, 0-9), slashes (/), or dashes (-)] by prefixing it with `fig:`, e.g., the label of a figure originating from the chunk `foo` will be `fig:foo`. To reference a figure, use the syntax `\@ref(label)`^[do not forget the leading backslash!], where `label` is the figure label, e.g., `fig:foo`. For example, the following code chunk was used to produce Figure \@ref(fig:plot).

    ```{r plot, fig.cap = "Regular figure. The first sentence...", echo = FALSE}
    plot(cars)
    ```

<div class="figure">
<img src="./pages/mydoc/AuthoringRmdVignettes_files/plot-1.png" alt="Regular figure. The first sentence of the figure caption is automatically emphasized to serve as figure title." width="100%" />
<p class="caption">(\#fig:plot)Regular figure. The first sentence of the figure caption is automatically emphasized to serve as figure title.</p>
</div>

In addition to regular figures, *BiocStyle* provides small and wide figures which can be specified by `fig.small` and `fig.wide` code chunk options.
Wide figures are left-aligned with the paragraph and extend on the right margin, as Figure \@ref(fig:wide).
Small figures are meant for possibly rectangular plots which are centered with respect to the text column, see Figure \@ref(fig:small).

<div class="figure">
<img src="./pages/mydoc/AuthoringRmdVignettes_files/wide-1.png" alt="Wide figure. A plot produced by a code chunk with option `fig.wide = TRUE`." width="100%"  class="widefigure" />
<p class="caption">(\#fig:wide)Wide figure. A plot produced by a code chunk with option `fig.wide = TRUE`.</p>
</div>

<div class="figure">
<img src="./pages/mydoc/AuthoringRmdVignettes_files/small-1.png" alt="Small figure. A plot produced by a code chunk with option `fig.small = TRUE`." width="100%"  class="smallfigure" />
<p class="caption">(\#fig:small)Small figure. A plot produced by a code chunk with option `fig.small = TRUE`.</p>
</div>


<br><br><center><a href="mydoc_AuthoringRmdVignettes_05.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_AuthoringRmdVignettes_07.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
