'Rsubread_idx

Usage: 
      Rsubread_idx.doc.R --reference <file>
    
Options:
    --reference <file> Named character vector containing paths to the reference sequence.
' -> doc

library(docopt)
library(systemPipeR)
library(Rsubread)

opts <- docopt(doc)
print(opts$reference)

buildindex(basename=opts$reference, reference=opts$reference) 
print("Index created successfully")