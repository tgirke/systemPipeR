'preprocessReads

Usage:
      preprocessReads_se.doc.R [--FileName1 <file>] [--FileName2 <file>] [--outfile1 <file>] [--outfile2 <file>] [--Fct <string>] [--batchsize <int>] [--overwrite <logical>] [--compress <logical>]

Options:
    --FileName1 <file> FileName1 file path.
    --FileName2 <file> FileName2 file path.
    --outfile1 <file> outfile1 file path.
    --outfile2 <file> outfile2 file path.
    --Fct <string> character string of custom read preprocessing function call where both the input and output needs to be an object of class ShortReadQ. The name of the input ShortReadQ object needs to be fq.
    --batchsize <int> Number of reads to process in each iteration by the internally used FastqStreamer function.
    --overwrite <logical> If TRUE existing file will be overwritten.
    --compress <logical> If TRUE existing file will be overwritten.
    ' -> doc

library(docopt)
suppressPackageStartupMessages(library(systemPipeR))
suppressPackageStartupMessages(library(ShortRead))
opts <- docopt(doc)

print(opts$FileName1)
print(opts$outfile1)
print(opts$FileName2)
print(opts$outfile2)
print(opts$Fct)
batchsize <- as.integer(opts$batchsize)
overwrite <- as.logical(opts$overwrite)
compress <- as.logical(opts$compress)

preprocessReads(FileName1 = opts$FileName1,
                FileName2 = opts$FileName2,
                outfile1 = opts$outfile1,
                outfile2 = opts$outfile2,
                Fct = opts$Fct, batchsize = batchsize,
                overwrite = overwrite,
                compress = compress)

print("preprocessReads executed successfully")

