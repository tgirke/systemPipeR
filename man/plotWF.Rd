\name{plotWF}
\alias{plotWF}
\title{
Graphviz Plot Workflow
}
\description{
Graphviz Plot Workflow
}
\usage{
plotWF(sysargslist, plot_style = "detect", out_type = "html", out_path = "default", height = NULL, width = NULL)
}
\arguments{
  \item{sysargslist}{
data frame, the standard df generated from subsetRmd function
}
  \item{plot_style}{
  one of 'detect', 'none', 'linear'
}
  \item{out_type}{
  choose from 'html', 'svg', 'png'
}
  \item{out_path}{
string, path of output plot, only apply from svg or png
}
  \item{height}{
int, height of svg or png in pixels, default NULL, automatic
}
  \item{width}{
int, width of svg or png in pixels, default NULL, automatic
}
}
\author{
Le Zhang and Daniela Cassol
}
\examples{
\dontrun{
## Construct SYSargsList object from Rmd file
library(systemPipeRdata)
script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- tempdir()
SYSconfig <- initProject(projPath=dir_path, targets=targets, script=script, overwrite = TRUE)
sysargslist <- initWF(sysconfig ="SYSconfig.yml")
plotWF(sysargslist, plot_style = "none", out_type = "png")
}
}
\keyword{ utilities }
