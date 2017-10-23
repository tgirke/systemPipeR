---
title: 8. Equations
last_updated: Sat Oct 21 17:03:44 2017
sidebar: mydoc_sidebar
permalink: mydoc_AuthoringRmdVignettes_08.html
---

To number and reference equations, put them in equation environments and append labels to them following the syntax `(\#eq:label)`^[due to technical constraints equation labels must start with `eq:`], e.g.,

```latex
\begin{equation}
  f\left(k\right) = \binom{n}{k} p^k\left(1-p\right)^{n-k}
  (\#eq:binom)
\end{equation}
```
renders the equation below.

$$
  f\left(k\right) = \binom{n}{k} p^k\left(1-p\right)^{n-k}
  (\#eq:binom)
$$

You may then refer to Equation \@ref(eq:binom) by `\@ref(eq:binom)`. Note that in HTML output only labeled equations will appear numbered.


<br><br><center><a href="mydoc_AuthoringRmdVignettes_07.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_AuthoringRmdVignettes_09.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
