'gsnap_idx

Usage: 
      gsnap_idx.doc.R --reference <file>
    
Options:
    --reference <file> Named character vector containing paths to the reference sequence.
' -> doc

library(docopt)
library(systemPipeR)
library(gmapR)

opts <- docopt(doc)
print(opts$reference)

gmapGenome <- GmapGenome(normalizePath(opts$reference), directory="data", name="gmap_tair10chr/", create=TRUE, k=12L)
save(gmapGenome, file = paste0("data/gmapGenome.RData"))
print("Index created successfully")

