\name{runDiff}
\alias{runDiff}
\title{
Differential abundance analysis for many range sets
}
\description{
Convenience wrapper function for \code{run_edgeR} and \code{run_DESeq2} to
perform differential expression or abundance analysis iteratively for several
count tables. The latter can be peak calling results for several samples or
counts generated for different genomic feature types. The function also returns 
the filtering results and plots from \code{filterDEGs}. 
}
\usage{
runDiff(args, outfiles=NULL, diffFct, targets, cmp, dbrfilter, ...)
}
\arguments{
  \item{args}{
    An instance of \code{SYSargs} or \code{SYSargs2} constructed from a
    \code{targets} file where the first column (\code{targetsin(args)} or 
    \code{targets.as.df(targets(args))}) contains the paths to the  tabular 
    read count data files. Another possibily is \code{named character vector}
    with the paths to the tabular range data files and the elements names should
    be the sampleID.
}
  \item{outfiles}{
    Default is \code{NULL}. When \code{args} is an object of \code{named 
    character vector} class, \code{outfile} argument is required. 
    \code{Named character vector} with the paths to the resulting count tables
    and the elements names should be the sampleID.
}
  \item{diffFct}{
    Defines which function should be used for the differential abundance analysis.
    Can be \code{diffFct=run_edgeR} or \code{diffFct=run_DESeq2}.
}
  \item{targets}{
    targets \code{data.frame}
}
  \item{cmp}{
    code{character matrix} where comparisons are defined in two columns. This matrix should be generated with \code{readComp()} from the targets file. Values used for comparisons need to match those in the \code{Factor} column of the targets file.
}
  \item{dbrfilter}{
    Named vector with filter cutoffs of format \code{c(Fold=2, FDR=1)} where \code{Fold} refers to the fold change cutoff (unlogged) and \code{FDR} to the p-value cutoff. Those values are passed on to the \code{filterDEGs} function.
}
  \item{\dots}{
    Arguments to be passed on to the internally used \code{run_edgeR} or 
    \code{run_DESeq2} function.
}
}
\value{
Returns \code{list} containing the \code{filterDEGs} results for each
count table. Each result set is a \code{list} with four components 
which are described under \code{?filterDEGs}. The result files
contain the \code{edgeR} or \code{DESeq2} results from the comparisons
specified under \code{cmp}. The base names of the result files are the 
same as the corresponding input files specified under \code{countfiles} 
and the value of \code{extension} appended. 
}
\author{
Thomas Girke
}
\seealso{
\code{run_edgeR}, \code{run_DESeq2}, \code{filterDEGs}
}
\examples{
## Paths to BAM files
param <- system.file("extdata", "bowtieSE.param", package="systemPipeR")
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
args_bam <- systemArgs(sysma=param, mytargets=targets)
bfl <- BamFileList(outpaths(args_bam), yieldSize=50000, index=character())

\dontrun{
## SYSargs with paths to range data and count files
args <- systemArgs(sysma="param/count_rangesets.param", mytargets="targets_macs.txt")

## Iterative read counting
countDFnames <- countRangeset(bfl, args, mode="Union", ignore.strand=TRUE)
writeTargetsout(x=args, file="targets_countDF.txt", overwrite=TRUE)

## Run differential abundance analysis
cmp <- readComp(file=args_bam, format="matrix")
args_diff <- systemArgs(sysma="param/rundiff.param", mytargets="targets_countDF.txt")
dbrlist <- runDiff(args, diffFct=run_edgeR, targets=targetsin(args_bam), cmp=cmp[[1]], independent=TRUE, dbrfilter=c(Fold=2, FDR=1))
writeTargetsout(x=args_diff, file="targets_rundiff.txt", overwrite=TRUE)
}
}
\keyword{ utilities }
