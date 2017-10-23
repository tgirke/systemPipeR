## ----macro, echo=FALSE---------------------------------------------------
macro <- function(name, pkg, description)
    sprintf('`%s("%s")` %s %s', name, pkg,
            description, do.call(name, list(pkg)))
inline <- function(cmd) sprintf('`` `r %s` ``', cmd)

## ----no-cap--------------------------------------------------------------
plot(cars)

## ----plot, fig.cap="Regular figure. The first sentence of the figure caption is automatically emphasized to serve as figure title.", echo=FALSE----
plot(cars)

## ----wide, fig.cap="Wide figure. A plot produced by a code chunk with option `fig.wide = TRUE`.", fig.wide=TRUE, echo=FALSE----
plot(cars)

## ----small, fig.cap="Small figure. A plot produced by a code chunk with option `fig.small = TRUE`.", fig.small=TRUE, echo=FALSE----
plot(cars)

## ----mtcars--------------------------------------------------------------
knitr::kable(
  head(mtcars[, 1:8], 10), caption = 'A table of the first 10 rows of `mtcars`.'
)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

