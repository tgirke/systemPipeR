\name{runDiff}
\alias{runDiff}
\title{
Iterative differential abundance analysis for many range sets
}
\description{
Convenience wrapper function for \code{run_edgeR} and \code{run_DESeq2} 
to iterate over several count tables.
}
\usage{
runDiff(countfiles, diffFct, targets, cmp, extension, dbrfilter, ...)
}
\arguments{
  \item{countfiles}{
Character vector containing paths to tabular count files. 
}
  \item{diffFct}{
Defines which function should be used for the differential abundance analysis.
Can be \code{diffFct=run_edgeR} or \code{diffFct=run_DESeq2}.
}
  \item{targets}{
	targets \code{data.frame}
}
  \item{cmp
	\code{character matrix} where comparisons are defined in two columns. This matrix should be generated with \code{readComp()} from the targets file. Values used for comparisons need to match those in the \code{Factor} column of the targets file.
}
  \item{extension}{
Chracter string specifying the file extensions of the result files.
}
  \item{dbrfilter}{
	Named vector with filter cutoffs of format \code{c(Fold=2, FDR=1)} where \code{Fold} refers to the fold change cutoff (unlogged) and \code{FDR} to the p-value cutoff.
}
  \item{\dots}{
Arguments to be passed on to the internally used \code{run_edgeR} or \code{run_DESeq2} function.
}
}
\value{
Returns \code{list} containing the \code{filterDEG} results for each
count table. Each result set is a \code{list} with four components 
which are described under \code{?filterDEG}. The result result files
contain the \code{edgeR} or \code{DESeq2} results from the comparisons
specified under \code{cmp}. The base names of the result files are the 
same as the corresponding input files specified under \code{countfiles} 
and the value of \code{extension} appended. 
}
\author{
Thomas Girke
}
\seealso{
\code{run_edgeR}, \code{run_DESeq2}, \code{filterDEG}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (countfiles, diffFct, targets, cmp, extension, dbrfilter, 
    ...) 
{
    absent_count_files <- countfiles[!file.exists(countfiles)]
    if (length(absent_count_files) != 0) 
        stop("The following files assigned to 'countfiles' do not exist: ", 
            paste(basename(absent_count_files), collapse = ", "))
    dbrDFnames <- paste0(countfiles, extension)
    names(dbrDFnames) <- names(countfiles)
    dbrlists <- list()
    for (i in seq(along = dbrDFnames)) {
        countDF <- read.delim(countfiles[i], row.names = 1)
        edgeDF <- diffFct(countDF = countDF, targets, cmp, ...)
        write.table(edgeDF, dbrDFnames[i], quote = FALSE, sep = "\t", 
            col.names = NA)
        pdf(paste0(dbrDFnames[i], ".pdf"))
        DBR_list <- filterDEGs(degDF = edgeDF, filter = dbrfilter)
        dev.off()
        dbrlists <- c(dbrlists, list(DBR_list))
        names(dbrlists)[i] <- names(dbrDFnames[i])
        cat("Wrote count result", i, "to", basename(dbrDFnames[i]), 
            "\n")
        cat("Saved plot", i, "to", basename(paste0(dbrDFnames[i], 
            ".pdf")), "\n")
    }
    return(dbrlists)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
