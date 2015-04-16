\name{run_DESeq2}
\alias{run_DESeq2}
\title{
Runs DESeq2 
}
\description{
Convenience wrapper function to identify differentially expressed genes
(DEGs) in batch mode with \code{DESeq2} for any number of pairwise sample
comparisons specified under the \code{cmp} argument.  Users are strongly
encouraged to consult the \code{DESeq2} vignette for more detailed information
on this topic and how to properly run \code{DESeq2} on data sets with more
complex experimental designs.
}
\usage{
run_DESeq2(countDF, targets, cmp, independent = FALSE)
}
\arguments{
  \item{countDF}{
\code{date.frame} containing raw read counts
}
  \item{targets}{
	targets \code{data.frame}
}
  \item{cmp}{
	\code{character matrix} where comparisons are defined in two columns. This matrix should be generated with the \code{readComp()} function from the targets file. Values used for comparisons need to match those in the \code{Factor} column of the targets file.
}
  \item{independent}{
	If \code{independent=TRUE} then the \code{countDF} will be subsetted for each comparison. This behavior can be useful when working with samples from unrelated studies. For samples from the same or comparable studies, the setting \code{independent=FALSE} is usually preferred.
}
}
\value{
\code{data.frame} containing \code{DESeq2} results from all comparisons. Comparison labels are appended to column titles for tracking.
}
\references{
Please properly cite the \code{DESeq2} papers when using this function:
http://www.bioconductor.org/packages/devel/bioc/html/DESeq2.html
}
\author{
Thomas Girke
}
\seealso{
\code{run_edgeR}, \code{readComp} and \code{DESeq2} vignette
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
countfile <- system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDF <- read.delim(countfile, row.names=1)
degseqDF <- run_DESeq2(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE)
pval <- degseqDF[, grep("_FDR$", colnames(degseqDF)), drop=FALSE]
fold <- degseqDF[, grep("_logFC$", colnames(degseqDF)), drop=FALSE]
DEG_list <- filterDEGs(degDF=degseqDF, filter=c(Fold=2, FDR=10))
names(DEG_list)
DEG_list$Summary
}
\keyword{ utilities }
