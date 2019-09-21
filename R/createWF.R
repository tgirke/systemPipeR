###################
## SYSarg2 empty ##
###################
## SYSarg2 empty
SYScreate <- function(){
  SYS.empty <- list(targets=data.frame(),  
            targetsheader=list(),
            modules=list(),
            wf=list(),
            clt=list(),
            yamlinput=list(), 
            cmdlist=list(),
            input=list(),
            output=list(),
            cwlfiles=list(), 
            inputvars=list()) 
 return((as(SYS.empty, "SYSargs2")))
}

## Usage:
# SYScreate()

###################
##   write.cwl   ##
###################
write.clt <- function(commandLine, cwlVersion, class, file.cwl) {
  cwlVersion <- cwlVersion 
  class <- class
  baseCommand <- commandLine$baseCommand[[1]]
  ##requeriments
  if(is.null(commandLine$requeriments)){
    dump <- "do nothing" 
  } else {
    requeriments <- list() ##TODO
  }
  ##ARGUMENTS
  if(is.null(commandLine$arguments)){
    dump <- "do nothing"
  } else {
    arguments <- sapply(seq_along(commandLine$arguments), function(x) list())
    for(i in seq_along(commandLine$arguments)){
      arguments[[i]]["prefix"] <- commandLine$arguments[[i]]$preF
      arguments[[i]]["valueFrom"] <- commandLine$arguments[[i]]$valueFrom
    }
  }
  ##INPUTS
  if(any(names(commandLine$inputs)=="")) stop("Each element of the list 'commandLine' needs to be assigned a name")
  if(is.null(names(commandLine$inputs))) stop("Each element of the list 'commandLine' needs to be assigned a name")
  if(!c("type") %in% names(commandLine$inputs[[i]])) stop("Each element of the sublist 'inputs' in 'commandLine' needs to be defined the type of the argument, for example: type='Directory' or type='File' or type='int' or type='string'")
  inputs <- sapply(names(commandLine$inputs), function(x) list())
  for(i in seq_along(commandLine$inputs)){
    if("type" %in% names(commandLine$inputs[[i]])){
      inputs[[i]]["type"] <-commandLine$inputs[[i]]$type 
    } 
    if("preF" %in% names(commandLine$inputs[[i]])){
      if(commandLine$inputs[[i]]$preF==""){
        inputs[[i]]["inputBinding"] <- list(list(prefix=NULL))
      } else{
        inputs[[i]]["inputBinding"] <- list(list(prefix=commandLine$inputs[[i]]$preF))
      }
    } 
    if(any(c("label", "secondaryFiles", "doc", "default", "format", "streamable") %in% names(commandLine$inputs[[i]]))){
     for(j in which(c("label", "secondaryFiles", "doc", "default", "format", "streamable") %in% names(commandLine$inputs[[i]]))){
       inputs[[i]][names(commandLine$inputs[[i]][j])] <- commandLine$inputs[[i]][names(commandLine$inputs[[i]])][[j]]
      }
    }
  }
  ##OUTPUTS
  outputs <- sapply(names(commandLine$outputs), function(x) list())
  if(!c("type") %in% names(commandLine$inputs[[i]])) stop("Each element of the sublist 'outputs' in 'commandLine' needs to be defined the type of the argument, for example: type='Directory' or type='File'.")
  if(all(!c("File", "Directory") %in% commandLine$outputs[[1]]$type)) stop("Each element of the sublist 'outputs' in 'commandLine' needs to be defined the type = 'Directory', 'File'.")
  for(i in seq_along(commandLine$outputs)){
    outputs[[i]]["type"] <- commandLine$outputs[[i]]$type
    outputs[[i]]["outputBinding"] <- list(list(glob=commandLine$outputs[[i]][[2]]))
  }
  clt <- list(cwlVersion=cwlVersion, class=class)
  if(exists("requeriments")){
    clt <- c(clt, list(requeriments=requeriments))
  }
  if(exists("arguments")){
    clt <- c(clt, list(arguments=arguments))
  }
  clt <- c(clt, list(baseCommand=baseCommand, inputs=inputs, outputs=outputs))
  
  yaml::write_yaml(x=clt, file = file.cwl) 
  clt <- list(clt)
  names(clt) <- baseCommand
  return(clt)
}

## Usage:
# clt_cwl <- write.clt(commandLine, cwlVersion, class, file.cwl) 

###################
##   write.yml   ##
###################
## Write the yaml file  
write.yml <- function(commandLine, file.yml, results_path, module_load){
  inputs <- commandLine$inputs
  if(any(names(inputs)=="")) stop("Each element of the list 'commandLine' needs to be assigned a name")
  if(is.null(names(inputs))) stop("Each element of the list 'commandLine' needs to be assigned a name")
 ##yamlinput_yml 
  yamlinput_yml <- sapply(names(inputs), function(x) list())
  for(i in seq_along(inputs)){
    if(!c("type") %in% names(inputs[[i]])) stop("Each element of the sublist 'inputs' in 'commandLine' needs to be defined the type of the argument, for example: type='Directory' or type='File' or type='int' or type='string'")
    if("type" %in% names(inputs[[i]])){
      if(any(c("File", "Directory") %in% inputs[[i]])){
        yamlinput_yml[[i]]["class"] <- inputs[[i]]$type
        yamlinput_yml[[i]]["path"] <- inputs[[i]]$yml
      } else if (any(c("int", "string") %in% inputs[[i]])){
        yamlinput_yml[[i]]  <- inputs[[i]]$yml
      }
    } else {
      print("do something")
    }
  }
  ## results_path
  yamlinput_yml[["results_path"]]["class"] <- list("Directory")
  yamlinput_yml[["results_path"]]["path"] <- list(results_path)
  ## moduleload
  for(i in seq_along(module_load)){
    yamlinput_yml[["ModulesToLoad"]][paste0("module", i)] <- list(module_load[[i]])
  }
  ## write out the '.yml' file
  yaml::write_yaml(x=yamlinput_yml, file = file.yml) 
  return(yamlinput_yml)
}

## Usage: 
# yamlinput_yml <- write.yml(commandLine, file.yml, results_path, module_load) 

###################################################
##   Create CommandLineTools from Command-line   ##
###################################################
createWF <- function(commandLine, results_path="./results", module_load="baseCommand", file = "default", overwrite = FALSE, cwlVersion = "v1.0", class = "CommandLineTool"){
## TODO if is not default, module load and file
  if(!class(commandLine)=="list") stop("'commandLine' needs to be object of class 'list'.")  
  if(any(!c("baseCommand", "inputs", "outputs") %in% names(commandLine))) stop("Argument 'commandLine' needs to be assigned at least to: 'baseCommand', 'input' or 'output'.")
  if(all(!c("Workflow", "CommandLineTool") %in% class)) stop("Class slot in '<wf_file>.cwl' needs to be 'Workflow' or 'CommandLineTool'.")
  if(dir.exists(results_path)==FALSE) dir.create(path=results_path)
  ## module_load 
  ## 1.  module_load="baseCommand" will use the name of the software 
  ## 2.  module_load=c("ncbi-blast/2.2.30+", "hisat2/2.1.0") will use the specific version and names
  if(module_load == "baseCommand"){
    module_load <- commandLine$baseCommand[[1]]
  } else {
      module_load <- module_load
  }
  ## File Path
  ## 1. file="default"
  ## 2. file = c("test.cwl", "test.yml")
  if("default" %in% file){
    file.cwl <- paste("param/cwl/", commandLine$baseCommand, "/", commandLine$baseCommand, ".cwl", sep="")
    file.yml <- paste("param/cwl/", commandLine$baseCommand, "/", commandLine$baseCommand, ".yml", sep="")
  } else {
    for(i in seq_along(file)){
      extension <- sub('.*\\.', '', file[[i]])
      if(!c("cwl") %in% extension & !c("yml") %in% extension) stop ("Argument 'file' needs to be assigned as a character vector with the names of the two param file. For example, 'test.cwl' and 'test.yml'.")
      if(c("yml") %in% extension){
        file.yml <- file[[i]]
      } else if(c("cwl") %in% extension){
        file.cwl <- file[[i]]
      }
    }
  }
  if(file.exists(file.cwl) & overwrite == FALSE) 
    stop(paste("I am not allowed to overwrite files; please delete existing file:", 
               file, "or set 'overwrite=TRUE', or provide a different name in the 'file' argument"))
  if(file.exists(file.yml) & overwrite == FALSE) 
                 stop(paste("I am not allowed to overwrite files; please delete existing file:", 
                            file, "or set 'overwrite=TRUE', or provide a different name in the 'file' argument"))
  ## class("CommandLineTool", "Workflow")
  ## TODO: Expand to write.WF()
  WF.temp <- SYScreate()
  WF.temp <- as(WF.temp, "list")
  WF.temp$wf <- list()
  WF.temp$clt <- write.clt(commandLine, cwlVersion, class, file.cwl) 
  WF.temp$yamlinput <- write.yml(commandLine, file.yml, results_path, module_load) 
  WF.temp$modules <- WF.temp$yamlinput$ModulesToLoad
  WF.temp$cmdlist <- sapply(names(clt_cwl), function(x) list(NULL))
  WF.temp$input <- sapply(names(clt_cwl), function(x) list(NULL))
  WF.temp$output <- sapply(names(clt_cwl), function(x) list(NULL))
  WF.temp$cwlfiles <- list(cwl=normalizePath(file.path(file.cwl)), yml=normalizePath(file.path(file.yml)), steps=names(clt_cwl))
  return(as(WF.temp, "SYSargs2"))
  ##THE END  
}

## Usage: 
# WF <- createWF(commandLine, results_path="./results", module_load="baseCommand", file = c("test.cwl", "test.yml"),
#               overwrite = FALSE, cwlVersion = "v1.0", class = "CommandLineTool")
# renderWF(WF)
# 
# 
# 
# ## Command-line 
# ## Example 1
# # blastp -query ./data/myseq.fasta -db ./data/halobacterium.faa -outfmt 0 -evalue 1e-6 -out blastp.out
# baseCommand <- list(baseCommand="blastp")
# inputs <- list(
#   "query"=list(type="File", preF="-query", yml="_FASTA_"), 
#   "db"= list(type="string", preF="-db", yml="_database.faa_"), 
#   "outfmt"= list(type="int", preF="-outfmt", yml=0L), 
#   "evalue"=list(type="int", preF="-evalue", yml=1e-6), 
#   "out"=list(type="File", preF="-out", yml="blastp.out"))
# outputs <- list("out"=list(type="File", "results/blastp.out"))
# commandLine <- list(baseCommand=baseCommand, inputs=inputs, outputs=outputs)
# 
# ## Example 2
# # "hisat2 -S ./results/_SampleName_.sam  -x ./data/tair10.fasta  -k 1  --min-intronlen 30  --max-intronlen 3000 --threads 4 -U _FASTQ_PATH1_ "
# 
# baseCommand <- list(baseCommand="hisat2")
# inputs <- list(
#   "S"=list(type="File", preF="-S", yml="./results/_SampleName_.sam"),
#   "x"=list(type="File", preF="-x", yml="./data/tair10.fasta"), 
#   "k"= list(type="int", preF="-k", yml= 1L), 
#   "threads"= list(type="int", preF="-threads", yml=4L),
#   "min-intronlen"= list(type="int", preF="--min-intronlen", yml= 30L),  
#   "max-intronlen"= list(type="int", preF="--max-intronlen", yml=3000L), 
#   "U"=list(type="File", preF="-U", yml="./data/_FASTQ_PATH1_") )
# 
# outputs <- list("hisat2_sam"=list(type="File", "./results/_SampleName_.sam"))
# commandLine <- list(baseCommand=baseCommand, inputs=inputs, outputs=outputs)
# 
# ## Example 3
# # "trim_galore --phred33  --gzip  -o ./results  --basename M1A  --length 25 --quality 20 -a AGATCGGAAGAGC ./data/SRR446027_1.fastq.gz"
# 
# baseCommand <- list(baseCommand="trim_galore")
# inputs <- list(
#   "SampleName"= list(type="string", yml="_SampleName_"), 
#   "length"= list(type="int", preF="--length", yml=25L),
#   "quality"= list(type="int", preF="--quality", yml= 20L),  
#   "adapter_trim_galore"= list(type="string", preF="-a", yml="AGATCGGAAGAGC"), 
#   "fq1"=list(type="File", preF="", yml="_FASTQ_PATH1_"), #achar uma forma de "prefix:" preF="" if ="" do name only
#   "results_path"=list(type="Directory", yml="./results/")
# )
# arguments <- list(
#   "phred33"=list(preF="--phred33"),
#   "gzip"=list(preF="--gzip"),
#   "basename"=list(preF="--basename", valueFrom="$(inputs.SampleName)"),
#   "out"=list(preF="-o", valueFrom="$(inputs.results_path.path)")
# )
# outputs <- list(
#   "trim_galore"=list(type="File", "$(inputs.results_path.path)/$(inputs.SampleName)_trimmed.fq.gz"), 
#   "trim_galore_report"=list(type="File", "$(inputs.results_path.path)/$(inputs.fq1.basename)_trimming_report.txt")
# )
# commandLine <- list(baseCommand=baseCommand, arguments=arguments, inputs=inputs, outputs=outputs)