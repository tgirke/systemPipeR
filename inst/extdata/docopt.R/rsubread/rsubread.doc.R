'Rsubread_map

Usage: 
      Rsubread.doc.R --output_file <file> --reference <file> --readfile1 <file> (--input_format <gzFASTQ> | <FASTQ> | <FASTA>) (--output_format <SAM> | <BAM>)  --nthreads <int> --indels <int> --TH1 <int> --maxMismatches <int>
    
Options:
    --output_file <file> A character vector specifying names of output files.
    --reference <file> Character string giving the basename of index file. Index files should be located in the current directory.
    --readfile1 <file>  A character vector including names of files that include sequence reads to be aligned.
    --input_format (<gzFASTQ> | <FASTQ> | <FASTA>)
    --output_format (<SAM> | <BAM>)
    --nthreads <int> 
    --indels <int> 
    --TH1 <int>
    --maxMismatches <int>
' -> doc
library(docopt)
library(systemPipeR)
library(Rsubread)

opts <- docopt(doc)
print(opts$readfile1)

align(index=opts$reference, readfile1=opts$readfile1, input_format=opts$input_format, 
      output_file=opts$output_file, output_format=opts$output_format, nthreads=opts$nthreads, 
      indels=opts$indels, TH1=opts$TH1, maxMismatches=opts$maxMismatches)
print("Mapping created successfully")