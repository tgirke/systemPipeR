###########
## Class ##
###########

####################
## Main functions ##
####################

###############
## Utilities ##
###############

#########################################################################################
## Function to check if the command line / Software is installed and set in your PATH ##
#########################################################################################
tryCL <- function(command){
  tryCatch({ system(command, ignore.stdout = TRUE, ignore.stderr = TRUE); print("All set up, proceed!") },
            warning = function(w) cat(paste0("ERROR: ", "\n", command, ": COMMAND NOT FOUND. ", '\n',
                                             "Please make sure to configure your PATH environment variable according to the software in use."), "\n"))
}

## Usage:
# tryCL(command="R") 
# tryCL(command="blastp") 
# tryCL(command="hisat2") 

########################################################
## Function to check if the Path (file or dir) exists ##
########################################################
tryPath <- function(path){
  tryCatch(normalizePath(path),
           warning = function(w) message(paste0(path,": ", "No such file or directory")), 
           error = function(e) message(paste0(path,": ", "Please provide a valid Path")))
}

## Usage:
# tryPath(path="./")

###########################################################################################
## Function to evaluate (eval=TRUE) or not evaluate (eval=FALSE) R code in the Rmarkdown ##
###########################################################################################
## infile: Name and path of the infile file, format Rmarkdown.
evalCode <- function(infile, eval=TRUE, output){
  if (!file.exists(infile)) stop("Provide valid 'infile' file. Check the file PATH.")
  if (!"Rmd" %in% .getExt(infile)) stop("Provide Rmarkdown file. Check the file PATH.")
  file <- readLines(infile)
  Rchunk <-  grep("^```[{}]", file)
  if (length(Rchunk) == 0) stop("The file does not have any R Chuck")
  for(i in Rchunk){
    if(grepl("eval", file[i])){
      file[i] <- sub("eval.*E", paste0("eval=", eval), file[i])
    } else {
      file[i] <- sub("[}]", paste0(", eval=", eval, "}"), file[i])
    }
  }
  writeLines(file, output)
}

## Usage:
# file <- system.file("extdata/workflows/rnaseq/", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
# evalCode(infile=file, eval=FALSE, output="test.Rmd")

###############################
## Return the file extension ##
###############################
## Function to return the extension of the file. The argument 'x' is a character vector or an object containing the file PATH.
.getExt <- function(x) {
  ext <- strsplit(basename(x), split="\\.")[[1]]
  ext <- ext[[-1]]
  return(ext)
}

## Usage:
# x <- system.file("extdata", "targets.txt", package = "systemPipeR")
# .getExt(x) ## internal fct: it will delete any suffixes up to the last slash ('/') character and return the 'extension of the file'

