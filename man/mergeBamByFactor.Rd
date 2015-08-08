\name{mergeBamByFactor}
\alias{mergeBamByFactor}
\title{
Merge BAM files based on factor
}
\description{
Merges BAM files based on sample groupings provided by a factor using internally
the \code{mergeBam} function from the \code{Rsamtools} package. The function also returns
an updated \code{SYSargs} object containing the paths to the merged BAM files as well as 
to the unmerged BAM files if any. All rows of merged samples are removed. 

The functionality provided by \code{mergeBamByFactor} is useful for experiments where 
pooling of replicates is advantageous to maximize the depth of coverage, such as prior to 
peak calling in ChIP-Seq or miRNA gene prediction experiments.
}
\usage{
mergeBamByFactor(args, mergefactor = targetsin(args)$Factor, overwrite = FALSE, silent = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{args}{
    An instance of \code{SYSargs} constructed from a \code{targets} file where the first 
column (\code{targetsin(args)}) contains the paths to the BAM files along with the column title 
\code{FileName}.
}
  \item{mergefactor}{
    \code{factor} containing the grouping information required for merging the BAM files 
referenced in the first column of \code{targetsin(args)}. The default uses 
\code{targetsin(args)$Factor}. Typically, the latter will merge BAM files
for which replicates are available.  
}
  \item{overwrite}{
    If \code{overwrite=FALSE} existing BAM files of the same name will not be overwritten. 
}
  \item{silent}{
    If \code{silent=TRUE} print statements will be suppressed.
}
  \item{\dots}{
    To pass on additional arguments to \code{mergeBam} from \code{Rsamtools}.
}
}
\value{
The merged BAM files will be written to output files with the following naming
convention: <first_BAM_file_name>_<grouping_label_of_factor>.<bam>. In
addition, the function returns an updated \code{SYSargs} object where all
output file paths contain the paths to the merged BAM files. The rows of the
parent/merged samples are removed and the rows of the unmerged samples remain unchanged.
}
\author{
Thomas Girke
}
\seealso{
\code{writeTargetsout}, \code{writeTargetsRef}
}
\examples{
## Construct initial SYSargs object
targetspath <- system.file("extdata", "targets_chip.txt", package="systemPipeR")
parampath <- system.file("extdata", "bowtieSE.param", package="systemPipeR")
args <- systemArgs(sysma=parampath, mytargets=targetspath)

\dontrun{
## After running alignmets (e.g. with Bowtie2) generate targets file 
## for the corresponding BAM files. The alignment step is skipped here.
writeTargetsout(x=args, file="targets_bam.txt", overwrite=TRUE) 
args <- systemArgs(sysma=NULL, mytargets="targets_bam.txt")

## Merge BAM files and return updated SYSargs object
args_merge <- mergeBamByFactor(args, overwrite=TRUE, silent=FALSE)

## Export modified targets file 
writeTargetsout(x=args_merge, file="targets_mergeBamByFactor.txt", overwrite=TRUE) 
}
}
\keyword{ utilities }
