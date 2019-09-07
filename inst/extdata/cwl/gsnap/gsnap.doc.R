'gsnap_map

Usage: 
      gsnap.doc.R  --readfile1 <file> --readfile2 <file> --output_file <file> (--molecule <DNA> | <RNA>) --max_mismatches <int> --nthreads <int> 
    
Options:
    --output_file <file> A character vector specifying names of output files.
    --readfile1 <file>  A character vector including names of files that include sequence reads to be aligned.
    --readfile2 <file>  A character vector including names of files that include sequence reads to be aligned.
    --molecule (<DNA> | <RNA>)
    --max_mismatches <int> 
    --nthreads <int> 
    ' -> doc

library(docopt)
library(systemPipeR)
library(gmapR)

opts <- docopt(doc)
print(opts$output_file)


load('./data/gmapGenome.RData')

if(exists("gmapGenome")==FALSE){
  stop("No index was detect. An index needs to be built before read mapping can be performed.")
} else {
  p <- GsnapParam(genome=gmapGenome, unique_only=TRUE, molecule=opts$molecule, max_mismatches=opts$max_mismatches)
  o <- gsnap(input_a=opts$readfile1, input_b=opts$readfile2, params=p, output=opts$output_file)
}
print("Mapping created successfully")