\name{mergeBamByFactor}
\alias{mergeBamByFactor}
\title{
Merge BAM files based on factor
}
\description{
Merges BAM files based on sample groupings provided by a factor using internally
the \code{mergeBam} function from the \code{Rsamtools} package. The function 
also returns an updated \code{SYSargs} or \code{SYSargs2} object containing
the paths to the merged BAM files as well as to the unmerged BAM files if there
are any. All rows of merged parent samples are removed. When a \code{named 
character vector} is provided as input, a \code{data.frame} with a target 
containing the paths to the merged BAM files as output. 

The functionality provided by \code{mergeBamByFactor} is useful for experiments 
where pooling of replicates is advantageous to maximize the depth of read 
coverage, such as prior to peak calling in ChIP-Seq or miRNA gene prediction 
experiments.
}
\usage{
mergeBamByFactor(args, targetsDF = NULL, mergefactor = "Factor", 
                overwrite = FALSE, silent = FALSE, ...)
}
\arguments{
  \item{args}{
    An instance of \code{SYSargs} or \code{SYSargs2} constructed from a
    \code{targets} file where the first column (\code{targetsin(args)} or 
    \code{targets.as.df(targets(args))}) contains the paths to the BAM files 
    along with the column title \code{FileName}. Another possibily is 
    \code{named character vector} with BAM files PATH and the elements names 
    should be the sampleID.
}
  \item{targetsDF}{
    This argument is required when \code{named character vector} is provided as 
    input. Default is \code{NULL}. 
    Object of class \code{DFrame}, and it can be obtained with
    \code{targetsWF(<SYSargsList>)}. 
}
  \item{mergefactor}{
    \code{factor} containing the grouping information required for merging the
    BAM files referenced in the first column of \code{targetsin(args)} or 
    \code{targets.as.df(targets(args))}. The default uses \code{Factor} column
    from the \code{targets} files as factor. The latter merges BAM files for 
    which replicates are specified in the \code{Factor} column.  
}
  \item{overwrite}{
    If \code{overwrite=FALSE} existing BAM files of the same name will not be
    overwritten. 
}
  \item{silent}{
    If \code{silent=TRUE} print statements will be suppressed.
}
  \item{\dots}{
    To pass on additional arguments to the internally used \code{mergeBam} 
    function from \code{Rsamtools}.
}
}
\value{
The merged BAM files will be written to output files with the following naming
convention: <first_BAM_file_name>_<grouping_label_of_factor>.<bam>. In
addition, the function returns an updated \code{SYSargs} or \code{SYSargs2} 
object where all output file paths contain the paths to the merged BAM files. 
When a \code{named character vector} is provided as input, a \code{data.frame} 
with a target containing the paths to the merged BAM files as output. 
The rows of the merged parent samples are removed and the rows of the unmerged
samples remain unchanged.
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
