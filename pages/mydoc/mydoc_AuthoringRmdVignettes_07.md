---
title: 7. Tables
last_updated: Sat Oct 21 17:03:44 2017
sidebar: mydoc_sidebar
permalink: mydoc_AuthoringRmdVignettes_07.html
---

Like figures, tables with captions will also be numbered and can be referenced. The caption is entered as a paragraph starting with `Table:`^[or just `:`], which may appear either before or after the table. When adding labels, make sure that the label appears at the beginning of the table caption in the form `(\#tab:label)`, and use `\@ref(tab:label)` to refer to it. For example, Table \@ref(tab:table) has been produced with the following code.

```markdown
## Fruit   | Price
bananas | 1.2
apples  | 1.0
oranges | 2.5
  
: (\#tab:table) A simple table. With caption.
```
        
## Fruit   | Price
bananas | 1.2
apples  | 1.0
oranges | 2.5
  
: (\#tab:table) A simple table. With caption.
    
The function `knitr::kable()` will automatically generate a label for a table environment, which is the chunk label prefixed by `tab:`, see Table \@ref(tab:mtcars).


```r
knitr::kable(
  head(mtcars[, 1:8], 10), caption = 'A table of the first 10 rows of `mtcars`.'
)
```



Table: (\#tab:mtcars)A table of the first 10 rows of `mtcars`.

##                       mpg   cyl    disp    hp   drat      wt    qsec   vs
Mazda RX4            21.0     6   160.0   110   3.90   2.620   16.46    0
Mazda RX4 Wag        21.0     6   160.0   110   3.90   2.875   17.02    0
Datsun 710           22.8     4   108.0    93   3.85   2.320   18.61    1
Hornet 4 Drive       21.4     6   258.0   110   3.08   3.215   19.44    1
Hornet Sportabout    18.7     8   360.0   175   3.15   3.440   17.02    0
Valiant              18.1     6   225.0   105   2.76   3.460   20.22    1
Duster 360           14.3     8   360.0   245   3.21   3.570   15.84    0
Merc 240D            24.4     4   146.7    62   3.69   3.190   20.00    1
Merc 230             22.8     4   140.8    95   3.92   3.150   22.90    1
Merc 280             19.2     6   167.6   123   3.92   3.440   18.30    1


<br><br><center><a href="mydoc_AuthoringRmdVignettes_06.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_AuthoringRmdVignettes_08.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
