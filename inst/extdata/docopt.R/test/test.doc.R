'test

Usage: 
      test.doc.R 
    
' -> doc

library(docopt)
library(systemPipeR)

targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")

write.table(targets, file = "results/targets_test.txt", sep = "\t", row.names = FALSE)

print("File created successfully")