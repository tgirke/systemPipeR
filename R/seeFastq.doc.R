'seeFastq

Usage: 
      seeFastq.doc.R --targets <file> --param <file> [--batchsize <number> --klength <number>] [--input_file <file>] [--dir_path <string>]
    
Options:
    
    --targets <file>      Named character vector containing paths to FASTQ file in the data fields and sample labels in the name slots.
    --param <file>      Named character vector containing paths to FASTQ file in the data fields and sample labels in the name slots.
    --batchsize <int>   Number of reads to random sample from each FASTQ file that will be considered in the QC analysis. Smaller numbers reduce the memory footprint and compute time.
    --klength <int>     Specifies the k-mer length in the plot for the relative k-mer diversity.
    --input_file <file>      Named character vector containing paths to FASTQ file in the data fields and sample labels in the name slots.
    --dir_path <file>      Named character vector containing paths to FASTQ file in the data fields and sample labels in the name slots.
' -> doc

library(docopt)
library(systemPipeR)


opts <- docopt(doc)
print(opts$input_file)

if(is.null(opts$input_file)==TRUE){
args <- systemArgs(sysma = opts$param, mytargets = opts$targets)
args
} else if(is.null(opts$input_file)==FALSE){
  args <- loadWorkflow(targets = opts$targets, wf_file = opts$param,
                       input_file = opts$input_file, dir_path = opts$dir_path)
  args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH_",
                                       SampleName = "_SampleName_"))
  args
}

.batchsize <- (opts$batchsize)
if (is.null(.batchsize)) {
  .batchsize <- 10000
} else {
  .batchsize <- as.numeric(opts$batchsize)
    }

.klength <- (opts$klength)
if (is.null(.klength)) {
  .klength <- 8
}


print(.batchsize)
fqlist <- seeFastq(fastq = infile1(args), batchsize = .batchsize, klength = .klength)
pdf("./results/fastqReport.pdf", height = 18, width = 4 * length(fqlist))
seeFastqPlot(fqlist)
dev.off()
