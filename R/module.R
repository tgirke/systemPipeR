#################################################
## Class and Method Definitions for EnvModules ##
#################################################
## Define EnvModules class
setClass("EnvModules", slots=c(
  available_modules="character",
  loaded_modules="character", 
  default_modules="list", 
  modulecmd="character"))

## Methods to return module components 
setGeneric(name="available_modules", def=function(x) standardGeneric("available_modules"))
setMethod(f="available_modules", signature="EnvModules", definition=function(x) {return(x@available_modules)})
setGeneric(name="loaded_modules", def=function(x) standardGeneric("loaded_modules"))
setMethod(f="loaded_modules", signature="EnvModules", definition=function(x) {return(module.List.Avail(action_type))})
setGeneric(name="default_modules", def=function(x) standardGeneric("default_modules"))
setMethod(f="default_modules", signature="EnvModules", definition=function(x) {return(x@default_modules)})
setGeneric(name="modulecmd", def=function(x) standardGeneric("modulecmd"))
setMethod(f="modulecmd", signature="EnvModules", definition=function(x) {return(x@modulecmd)})

#o method retorna a acao, e a os slots da class retornam as variaveis

## Constructor methods
## List to EnvModules with: as(mylist, "EnvModules")
setAs(from="list", to="EnvModules",  
      def=function(from) {
        new("EnvModules", 
            available_modules=from$available_modules,
            loaded_modules=from$loaded_modules,
            default_modules=from$default_modules,
            modulecmd=from$modulecmd)
      })

## Coerce back to list: as(EnvModules, "list")
setGeneric(name="EnvModules", def=function(x) standardGeneric("EnvModules"))
setMethod(f="EnvModules", signature="EnvModules", definition=function(x) {
  modules <- list(available_modules=x@available_modules, loaded_modules=x@loaded_modules, default_modules=x@default_modules,
                  modulecmd=x@modulecmd)
  return(modules)
})

## EnvModules to list with: as("EnvModules", list)
setAs(from="EnvModules", to="list", 
      def=function(from) {
        EnvModules(from) 
      })

## Define print behavior for EnvModules
setMethod(f="show", signature="EnvModules", 
          definition=function(object) {    
            cat(paste0("Instance of '", class(object), "':"), 
                paste0("   Slot names/accessors: "), 
                paste0("      avail: ", length(object@available_modules), 
                       "      list: ", length(object@loaded_modules)), 
                "\n", sep="\n")
          })

## Extend names() method
setMethod(f="names", signature="EnvModules",
          definition=function(x) {
            return(slotNames(x))
          })

# Behavior of "[" operator for EnvModules
setMethod(f="[", signature="EnvModules", definition=function(x, i, ..., drop) {
  if(is.logical(i)) {
    i <- which(i)
  }
  x@available_modules <- x@available_modules[i]
  x@loaded_modules <- x@loaded_modules[i]
  x@default_modules <- x@default_modules[i]
  x@modulecmd <- x@modulecmd[i]
  return(x)
})

## Behavior of "[[" operator for EnvModules
setMethod(f="[[", signature=c("EnvModules", "ANY", "missing"),
          definition=function(x, i, ..., drop) {
            return(as(x, "list")[[i]]) 
          })

## Behavior of "$" operator for EnvModules
setMethod("$", signature="EnvModules",
          definition=function(x, name) { 
            slot(x, name)
          })

## Replacement method for EnvModules using "[" operator 
setReplaceMethod(f="[[", signature="EnvModules", definition=function(x, i, j, value) {
  if(i==1) x@default_modules <- value
  if(i==2) x@modulecmd <- value
  if(i=="default_modules") x@default_modules <- value
  if(i=="modulecmd") x@modulecmd <- value 
  return(x)
})


#################################
## Module function Constructor ##
#################################
module <- function(action_type, module_name=NULL){
  #  if(!action_type %in% c("load", "unload", "list", "avail", "init", "clear")) stop("That action is not supported.")
  mymodules <- list(
    available_modules = as.character(),
    loaded_modules = as.character(),
    default_modules = default_modules <- myEnvModules(),
    modulecmd = modulecmd_path <- is.modules.avail())
  switch(action_type,
         "load"   = module.Load.Unload(action_type, module_name),
         "unload" = module.Load.Unload(action_type, module_name),
         "list"   = mymodules$loaded_modules <- module.List.Avail(action_type, mymodules$modulecmd),
         "avail"  = mymodules$available_modules <- module.List.Avail(action_type, mymodules$modulecmd),
         "clear"  = module.Clear.Init(action_type, mymodules$modulecmd, mymodules$default_modules),
         "init"   = module.Clear.Init(action_type, mymodules$modulecmd, mymodules$default_modules),
         stop("That action is not supported."))
  return(as(mymodules, "EnvModules"))
}


## Usage: 
# moduleClear()

module.Clear.Init <- function(action_type, modulecmd_path, default_modules){
  if(action_type=="init"){
    ## Make sure to process default modules after the environment is set with the above loop
    for (x in seq_along(default_modules[[1]])){
      print(paste("Loading module", default_modules[[1]][x]))
      try(module.Load.Unload("load", default_modules[[1]][x], modulecmd_path))}
  } else if(action_type=="clear"){
    loaded_modules <-  strsplit(Sys.getenv("LOADEDMODULES"),":")
    if (length(loaded_modules[[1]]) > 0) {
      for (x in seq(1,length(loaded_modules[[1]]))){
        module_name <- loaded_modules[[1]][x]
        print(paste("Unloading module", module_name))
        try(module.Load.Unload("unload", module_name, modulecmd_path))
      }
    }
  }
}

################
## moduleInit ##
################
# Function to load default modules files
moduleInit <- function(action_type){
  return(module("init"))
}
## Usage:
moduleInit()

#################
## moduleClear ##
#################
# Unload all currently loaded modules
moduleClear <- function(){
  return(module("clear"))
}
moduleClear()

################
## moduleAvail ##
################
# Function to check module avail
moduleAvail <- function(){
  return(module("avail"))
}
## Usage:
moduleAvail()

##################
## moduleUnload ##
##################
# Unload module file
moduleUnload <- function(module_name){
  return(module("unload", module_name))
}
## Usage: 
moduleUnload()

#####################
## Legacy Wrappers ##
#####################

## List software available in module system
modulelist <- function() {
  return(module("list"))
}
## Usage: 
modulelist()

## Load software from module system
moduleload <- function(module_name) {
  return(module("load", module_name))
}
## Usage: 
moduleload()

##################################
## module accessories functions ##
#################################

######################
## is.modules.avail ##
######################
# Check if the Environment Modules is available in the system and returns module 
# function PATH.
is.modules.avail <- function(){
  ## Find path for module command
  modulecmd_path <- Sys.getenv('LMOD_CMD')
  if(modulecmd_path =="") {
    try(suppressWarnings(modulecmd_path <- system("which modulecmd", intern=TRUE, ignore.stderr=TRUE)),
        silent=TRUE)
  }
  ## "Environment Modules" is not available
  if(length(modulecmd_path) == 0){
    message("'Environment Modules is not available. Please make sure to configure 
            your PATH environment variable according to the software in use.'", "\n", 
            "Find more detail help('module')")
    ## "Environment Modules" is available and proceed the module action_type
  } else if (length(modulecmd_path) > 0) {
    return(modulecmd_path)
  }
}

## Usage: 
is.modules.avail()

##################
## myEnvModules ##
##################
# It returns a list with all default modules files.
myEnvModules <- function(){
  # Module function assumes are using a bash SHELL
  if(!system("echo $0", intern = TRUE) %in% c("bash", "sh")) warning("'module' function assumes you are using a bash SHELL and may not work as expected in other SHELL types.")
  # Get base environment from login profile
  base_env <- strsplit(system('bash -l -c "env"', intern = TRUE, ignore.stderr=TRUE),'\n')
  base_env <- strsplit(as.character(base_env),'=')
  # Iterate through base environment
  for (x in seq_along(base_env)) {
    # Set environment based on login profile
    if (base_env[[x]][1]=="LOADEDMODULES" || base_env[[x]][1]=="MODULESHOME" || base_env[[x]][1]=="MODULEPATH" || base_env[[x]][1]=="MODULES_DIR" || base_env[[x]][1]=="HPCC_MODULES"){
      if (base_env[[x]][1]=="LOADEDMODULES"){
        default_modules <- strsplit(base_env[[x]][2],":")
        names(default_modules) <- "default_modules"
      } else {
        l <- list(base_env[[x]][2])
        names(l) <- base_env[[x]][1]
        do.call(Sys.setenv, l)
      }
    }
  }
  return(default_modules)
}
## Usage: 
myEnvModules()

########################
## module.Load.Unload ##
########################
# Internal function to module load <module_name> and module unload <module_name>
module.Load.Unload <- function(action_type, module_name="", modulecmd_path){
  module_name <- paste(module_name, collapse=' ')
  #modulecmd_path <- is.modules.avail()
  # Use the low level C binary for generating module environment variables
  try(module_vars <- system(paste(modulecmd_path, 'bash', action_type, module_name), intern = TRUE))
  if (length(module_vars) > 0){
    for (y in seq(1,length(module_vars))) {
      # Separate environment variables
      module_var <- strsplit(module_vars,";")
      # Iterate through all environment variables
      for (x in seq(1,length(module_var[[y]]))) {
        # Isolate key, value pair
        evar <- module_var[[y]][x]
        # Filter export commands
        if (length(grep('^ *export',evar)) == 0 && length(evar) > 0) {
          # Seprate key and value
          evar <- strsplit(as.character(evar),'=')
          # Stip spaces at the end of the value
          evar_val <- gsub('[[:space:]]','',evar[[1]][2])
          # Remove extra backslashes
          l <- list(gsub('\\$','',evar_val))
          # Load dependant modules
          if (length(grep('^ *module',evar[[1]][1])) > 0){
            inner_module <- strsplit(evar[[1]][1]," ")
          }
          # Source environment
          else if (length(grep('^ *source',evar[[1]][1])) > 0){
            warning(paste0("Module uses a bash script to initialize, some software may not function as expected:\n\t",evar[[1]][1]))
          }
          # Unset variables that need to be unset
          else if(length(grep("^ *unset ",evar[[1]][1])) > 0){
            evar <- gsub("^unset (.*)$","\\1",evar[[1]][1])
            Sys.unsetenv(evar)
          } else {
            # Assign names to each value in list
            names(l) <- evar[[1]][1]
            # Set environment variable in current environment
            do.call(Sys.setenv, l)
          }
        }
      }
    }
  }
}

## Usage: 
module.Load.Unload("load", "samtools")
module.Load.Unload("unload", "samtools")

########################
## module.List.Avail ##
########################
# Internal function to module list and module avail
module.List.Avail <- function(action_type, modulecmd_path){
  #modulecmd_path <- is.modules.avail()
  rstudioEnv <- Sys.getenv("RSTUDIO")
  if(rstudioEnv==1){
    try(module_vars <- system2(modulecmd_path, paste("bash", action_type, "-t"), stdout=TRUE, stderr=TRUE))
  } else {
    try(module_vars <- system2("module", action_type, stdout=TRUE, stderr=TRUE))
    # if(length(module_vars) == 0){
    #   try(module_vars <- system2("module", paste(action_type,'-t'), stdout=TRUE, stderr=TRUE))
    # } 
  }
  # Return only the module names
  module <- as.character()
  module_vars <- gsub(".*:|\\(.*?\\)", "", module_vars) ## remove "Currently Loaded Modulefiles:" and "(default)"
  module_vars <- strsplit(module_vars, " ")
  for(i in seq_along(module_vars)){
    for(j in seq_along(module_vars[[i]])){
      module_vars[[i]][j] <-  gsub("[(].*|.*[)]|^---.*", "", module_vars[[i]][j])
      if(nchar(module_vars[[i]][j]) > 0){
        module <- c(module, module_vars[[i]][j])
      }
    }
  }
  return(module)
}

## Usage: 
# module.List.Avail("list")
# module.List.Avail("avail")



#################################
## Access module system from R ##
#################################

# # Define R6 class
# RenvModule <- setRefClass("RenvModule",
#                           fields=list(modulecmd_path="character"),
#                           methods=list(
#                             ## Main function to allow avail, list and list
#                             init=function(){
#                               # Module function assumes MODULEPATH and MODULEDIR are set in login profile
#                               # Get base environment from login profile
#                               base_env <- strsplit(system('bash -l -c "env"',intern = TRUE),'\n')
#                               base_env <- strsplit(as.character(base_env),'=')
# 
#                               # Iterate through base environment
#                               for (x in seq(1,length(base_env))) {
# 
#                                 # Set environment based on login profile
#                                 if (base_env[[x]][1]=="LOADEDMODULES" || base_env[[x]][1]=="MODULESHOME" || base_env[[x]][1]=="MODULEPATH" || base_env[[x]][1]=="MODULES_DIR" || base_env[[x]][1]=="HPCC_MODULES"){
#                                   if (base_env[[x]][1]=="LOADEDMODULES"){
#                                     default_modules <- strsplit(base_env[[x]][2],":")
#                                   }
#                                   else{
#                                     l <- list(base_env[[x]][2])
#                                     names(l) <- base_env[[x]][1]
#                                     do.call(Sys.setenv, l)
#                                   }
#                                 }
#                               }
# 
#                               # Make sure to process default modules after the environment is set with the above loop
#                               for (x in seq_along(default_modules[[1]])){
#                                 print(paste("Loading module", default_modules[[1]][x]))
#                                 try(load_unload("load", default_modules[[1]][x]))
#                               }
#                             },
# 
#                             # Return available modules or currently loaded modules
#                             avail_list=function(action_type){
#                               try(module_vars <- system2("module", action_type, stdout=TRUE, stderr=TRUE))
#                               if(length(module_vars) == 0){
#                                 try(module_vars <- system2("module", paste(action_type,'-t'), stdout=TRUE, stderr=TRUE))
#                               }
#                               # Return only the module names
#                               # return(module_vars[-grep(":$",module_vars)])
#                               module <- as.character()
#                               module_vars <- gsub(".*:|\\(.*?\\)", "", module_vars) ## remove "Currently Loaded Modulefiles:" and "(default)"
#                               module_vars <- strsplit(module_vars, " ")
#                               for(i in seq_along(module_vars)){
#                                 for(j in seq_along(module_vars[[i]])){
#                                   module_vars[[i]][j] <-  gsub("[(].*|.*[)]|^---.*", "", module_vars[[i]][j])
#                                   if(nchar(module_vars[[i]][j]) > 0){
#                                     module <- c(module, module_vars[[i]][j])
#                                   }
#                                 }
#                               }
#                               return(module)
#                             },
# 
#                             # Unload all currently loaded modules
#                             clear=function(action_type){
#                               loaded_modules <-  strsplit(Sys.getenv("LOADEDMODULES"),":")
#                               if (length(loaded_modules[[1]]) > 0) {
#                                 for (x in seq(1,length(loaded_modules[[1]]))){
#                                   module_name <- loaded_modules[[1]][x]
#                                   print(paste("Unloading module",module_name))
#                                   try(load_unload("unload",module_name))
#                                 }
#                               }
#                             },
# 
#                             # Load and unload actions are basically the same, set environment variables given by module command
#                             load_unload=function(action_type, module_name=""){
#                               module_name <- paste(module_name, collapse=' ')
# 
#                               # Use the low level C binary for generating module environment variables
#                               try(module_vars <- system(paste(modulecmd_path,'bash',action_type, module_name),intern = TRUE))
# 
#                               if (length(module_vars) > 0){
#                                 for (y in seq(1,length(module_vars))) {
#                                   # Separate environment variables
#                                   module_var <- strsplit(module_vars,";")
# 
#                                   # Iterate through all environment variables
#                                   for (x in seq(1,length(module_var[[y]]))) {
#                                     # Isolate key, value pair
#                                     evar <- module_var[[y]][x]
# 
#                                     # Filter export commands
#                                     if (length(grep('^ *export',evar)) == 0 && length(evar) > 0) {
#                                       # Seprate key and value
#                                       evar <- strsplit(as.character(evar),'=')
#                                       # Stip spaces at the end of the value
#                                       evar_val <- gsub('[[:space:]]','',evar[[1]][2])
#                                       # Remove extra backslashes
#                                       l <- list(gsub('\\$','',evar_val))
# 
#                                       # Load dependant modules
#                                       if (length(grep('^ *module',evar[[1]][1])) > 0){
#                                         inner_module <- strsplit(evar[[1]][1]," ")
#                                       }
#                                       # Source environment
#                                       else if (length(grep('^ *source',evar[[1]][1])) > 0){
#                                         warning(paste0("Module uses a bash script to initialize, some software may not function as expected:\n\t",evar[[1]][1]))
#                                       }
#                                       # Unset variables that need to be unset
#                                       else if(length(grep("^ *unset ",evar[[1]][1])) > 0){
#                                         evar <- gsub("^unset (.*)$","\\1",evar[[1]][1])
#                                         Sys.unsetenv(evar)
#                                       }
#                                       else {
#                                         # Assign names to each value in list
#                                         names(l) <- evar[[1]][1]
#                                         # Set environment variable in current environment
#                                         do.call(Sys.setenv, l)
#                                       }
#                                     }
#                                   }
#                                 }
#                               }
#                             }
#                           )
# )
# 
# #Define what happens based on action
# module <- function(action_type,module_name=""){
#   # Create instance of RenvModule
#   myEnvModules <- RenvModule()
# 
#   # Find path for module command
#   myEnvModules$modulecmd_path <- Sys.getenv('LMOD_CMD')
#   if(myEnvModules$modulecmd_path =="") {
#     #    if (nchar(sub("\\s+", "", myEnvModules$modulecmd_path)) == 0) {
#     try(
#       suppressWarnings(myEnvModules$modulecmd_path <- system("which modulecmd",intern=TRUE,ignore.stderr=TRUE)),
#       silent=TRUE
#     )
#   }
# 
#   # Only initialize module system if it has not yet been initialized and the module command exists
#   if ( Sys.getenv('MODULEPATH') == "" && length(myEnvModules$modulecmd_path) > 0) {
#     myEnvModules$init()
#   } else if (Sys.getenv('MODULEPATH') == "" && length(myEnvModules$modulecmd_path) == 0) {
#     stop("Cound not find the path of Environment Modules command \"modulecmd\" nor LMOD_CMD")
#   }
# 
#   switch(action_type,
#          "load"   = myEnvModules$load_unload(action_type,module_name),
#          "unload" = myEnvModules$load_unload(action_type,module_name),
#          "list"   = return(myEnvModules$avail_list(action_type)),
#          "avail"  = return(myEnvModules$avail_list(action_type)),
#          "clear"  = myEnvModules$clear(action_type),
#          "init"   = myEnvModules$init(),
#          stop("That action is not supported.")
#   )
# }
## Usage:
# module("load","tophat")
# module("load","tophat/2.1.1")
# module("list")
# module("avail")
# module("init")
# module("unload", "tophat")
# module("unload", "tophat/2.1.1")

# #####################
# ## Legacy Wrappers ##
# #####################
# ## List software available in module system
# modulelist <- function() {
#   print(module("avail"))
#   warning("The function modulelist will be deprecated in future releases, please refer to the documentation for proper useage.")
# }
# 
# ## Load software from module system
# moduleload <- function(module,envir="PATH") {
#   module("load",module)
#   warning("The function moduleload will be deprecated in future releases, please refer to the documentation for proper useage.")
# }
