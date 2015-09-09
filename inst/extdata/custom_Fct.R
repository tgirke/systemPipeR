#################################################
## Collection of Non-Exported Custom Functions ##
#################################################

######################
## Adaptor Trimming ##
######################
## Trims adaptors hierarchically from longest to shortest match from right end of read. 
## If 'internalmatch=TRUE' then internal matches will trigger the same behavior. The 
## argument minpatternlength defines shortest adaptor match to consider for reads 
## containing only partial adaptors at the right end.
.iterTrimbatch1 <- function(fq, pattern, internalmatch=FALSE, minpatternlength=8, Nnumber=1, polyhomo=100, minreadlength=18, maxreadlength) {
    mylength <- max(width(fq))
    if(internalmatch==TRUE) {
        regpart <- ".*" 
    } else {
        regpart <- "$" 
    }
    processedreads <- ShortReadQ()
    ## Iterative trimming based on regular expression
    for(j in seq(along=minpatternlength:(nchar(pattern)))-1) {
        mypattern <- paste(substring(pattern, 1, nchar(pattern)-j), regpart, sep="")
        d <- DNAStringSet(gsub(mypattern, "", sread(fq)))
        q <- BStringSet(quality(quality(fq)), 1, width(d))
        fq <- ShortReadQ(sread=d, quality=q, id=id(fq))
        processedreads <- append(processedreads, fq[width(fq) < mylength])
        fq <- fq[width(fq) >= mylength]
    }
    processedreads <- append(processedreads, fq)
    processedreads <- processedreads[width(processedreads)>=minreadlength & width(processedreads)<=maxreadlength] # Read length filter
    filter1 <- nFilter(threshold=Nnumber) # keep only reads with <= x Ns
    filter2 <- polynFilter(threshold=polyhomo, nuc=c("A", "C", "T", "G")) # remove reads with x or more of the same nucleotide
    filter <- compose(filter1, filter2) # Combine filters into one
    processedreads[filter(processedreads)]
}
## Usage:
# iterTrim <- ".iterTrimbatch(fq, pattern="ACACGTCT", minpatternlength=6, Nnumber=1, polyhomo=20, minreadlength=16) 
# preprocessReads(args=args, Fct="iterTrim", batchsize=100000, overwrite=TRUE, compress=TRUE)
