###########
## Class ##
###########


####################
## Main functions ##
####################


###############
## Utilities ##
###############

##########################################
## SYSarg2 OR SYSargsList empty object ##
#########################################
## [Internal] Creation of the SYSargs2 or SYSargsList empty object
SYScreate <- function(class){
  if(class=="SYSargs2"){
    SYS.empty <- list(
      targets=data.frame(),  
      targetsheader=list(),
      modules=list(),
      wf=list(),
      clt=list(),
      yamlinput=list(), 
      cmdlist=list(NULL),
      input=list(),
      output=list(),
      cwlfiles=list(), 
      inputvars=list())
    return((as(SYS.empty, "SYSargs2")))
  } else if(class == "SYSargsList"){
    SYS.empty <- list(
      sysconfig=list(),
      # initWF=list(), 
      codeSteps=list(),
      stepsWF=numeric(),
      # runWF=list(), #should only be a method--think about
      # logload=list(),
      # statusWF=list(),
      # statusWF=list(),
      # plotWF=list(),
      # renderReport=list(), #method?
      SYSargs2_steps=list(),
      statusWF=list(),
      projectWF=list())
    return((as(SYS.empty, "SYSargsList")))
  } else if(!class %in% c("SYSargs2", "SYSargsList")) stop("Argument 'args' needs to be assigned an character 'SYSargs2' OR 'SYSargsList'")
}

## Usage:
# list <- SYScreate("SYSargsList")
# list <- SYScreate("SYSargs2")

#########################################################################################
## Function to check if the command line / Software is installed and set in your PATH ##
#########################################################################################
tryCL <- function(command){
  tryCatch( { system(command, ignore.stdout = TRUE, ignore.stderr = TRUE); print("All set up, proceed!") },
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

#################################
## Unexported helper functions ##
#################################

##################################
## Return the path of the file  ##
##################################
## [x] A character vector or an object containing file PATH.
.getPath <- function(x) {
  x <- normalizePath(x) ## add latter
  if (!file.exists(x)) warning ("No such file or directory. Check the file PATH.")
  path_un <- unlist(strsplit(x, "/"))
  path <- path_un[path_un != basename(x)]
  path <- paste0(path, collapse = "/")
  return(path)
}

## Usage:
# x <- system.file("extdata", "targets.txt", package = "systemPipeR")
# .getPath(x) ## internal fct: it will delete any suffixes up to the last slash ('/') character and return the 'path of the file'
# basename(x) ## BiocGenerics pkg: it will delete any prefix up to the last slash ('/') character and return the 'name of the file'

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

##############################################
## Return the file name, without extension ##
#############################################
## [x] A character vector or an object containing file File name without extension, simmilar with 'basename'
.getFileName <- function(x) {
  if (!file.exists(x)) warning ("No such file or directory. Check the file PATH.")
  filename <- strsplit(basename(x), split="\\.")[[1]]
  filename <- filename[[-2]]
  return(filename)
}

## Usage:
# x <- system.file("extdata", "targets.txt", package = "systemPipeR")
# .getFileName(x) ## internal fct: it will delete any suffixes up to the last slash ('/') character and return the 'filename' of the file'

##########################################################
## Internal function to detect nested of the param list ##
##########################################################
.nest <- function(x) ifelse(is.list(x), 1L + max(sapply(x, .nest)), 0L)

## Usage:
# param <- list(results_new=list(class="Directory", path=8))
# nesting <- .nest(param[1])

############################################################################
## Internal function to replace or add the list values at the input file ##
############################################################################
.replace <- function(input, param) {
  for (i in seq_along(param)) {
    for (j in seq_along(param[[i]])) {
      nesting <- .nest(param[i])
      if (nesting == 1) {
        if (is.numeric(param[[i]][j][[1]])) {
          input[names(param)[i]] <- as.integer(param[[i]][j])
        } else {
          input[names(param)[i]] <- param[[i]][j]
        }
      } else if (nesting > 1) {
        if (is.numeric(param[[i]][j][[1]])) {
          input[[names(param[i])]][[names(param[[i]][j])]] <- as.integer(param[[i]][[j]])
        } else {
          input[[names(param[i])]][[names(param[[i]][j])]] <- (param[[i]][[j]])
          input[[names(param[i])]] <- as.list(input[[names(param[i])]])
        }
      }
    }
  }
  return(input)
}

## Usage:
# input_file="param/cwl/hisat2/hisat2-se/hisat2-mapping-se.yml"
# input <- yaml::read_yaml(file.path(input_file))
# param <- list(thread=14, test="test")
# .replace(input=input, param=param)

##############################
## .sysconfigCheck function ##
##############################
## Function to check with the paths provided on the sysconfig file exists.
.sysconfigCheck <- function(sysconfig){
  if(!file.exists(sysconfig)==TRUE) stop("Provide valid 'sysconfig' file. Check the file PATH.")
  sysconfig <- yaml::read_yaml(sysconfig)
  project <- list(project=sysconfig$project$path, data=sysconfig$data$path, param=sysconfig$param$path,
                  results=sysconfig$results$path)
  for(i in seq_along(project)){
    if(is.null(project[[i]])){
      warning(paste(names(project[i]), "directory is missing..."))
    } else {
      tryPath(path = project[[i]])
    }
  }
  if(!is.null(sysconfig$targets$path)) tryPath(path = sysconfig$targets$path)
  if(is.null(sysconfig$script)) warning("Workflow is missing...")
  if(!is.null(sysconfig$script$path)) tryPath(path = sysconfig$script$path)
}

## Uage: 
# .sysconfigCheck(sysconfig="SYSconfig.yml")

##############################
## .parse_step function ##
##############################
## internal parse function
.parse_step <- function(t_lvl, input_steps){
  t_lvl_name <- names(t_lvl)
  input_steps <- unlist(input_steps %>% stringr::str_remove_all(" ") %>% stringr::str_split(",") %>% list())
  # single steps
  nocolon_steps <- input_steps[stringr::str_which(input_steps, "^[^:]+$")]
  lapply(nocolon_steps, function(x) if (!any(t_lvl_name %in% x)) stop(paste('Step', x, 'is not found')))
  # dash linked steps
  dash_list <- NULL
  for (i in stringr::str_which(input_steps, ":")){
    dash_step <- unlist(stringr::str_split(input_steps[i], ":"))
    dash_parse <- unlist(lapply(dash_step, function(x) {
      which(t_lvl_name %in% x) %>% ifelse(length(.) > 0, ., stop(paste('Step', x, 'is not found')))
    })) %>% {
      t_lvl_name[.[1]: .[2]] 
    }
    dash_list <- append(dash_list, dash_parse)
  }
  # merge
  all_step_name <- unique(append(nocolon_steps, dash_list))
  # if upper level step is selected, all sub-level steps will be added
  unlist(lapply(all_step_name, function(x) stringr::str_which(t_lvl_name, paste0('^', x, '\\..*')))) %>% 
    append(which(t_lvl_name %in% all_step_name)) %>% 
    unique() %>% sort() %>% return()
}
