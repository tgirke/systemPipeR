###############################################
## Class and Method Definitions for SYSargs2 ##
###############################################
## Define SYSargs2 class
setClass("SYSargs2", slots = c(
  targets = "list",
  targetsheader = "list",
  modules = "list",
  wf = "list",
  clt = "list",
  yamlinput = "list",
  cmdlist = "list",
  input = "list",
  output = "list",
  cwlfiles = "list",
  inputvars = "list", 
  cmdToCwl = "list"
))
## Methods to return SYSargs2 components
setGeneric(name = "targets", def = function(x) standardGeneric("targets"))
setMethod(f = "targets", signature = "SYSargs2", definition = function(x) {
  return(x@targets)
})

setMethod(f = "targetsheader", signature = "SYSargs2", definition = function(x) {
  return(x@targetsheader)
})
setMethod(f = "modules", signature = "SYSargs2", definition = function(x) {
  return(setNames(as.character(x@modules), names(x@modules)))
})
setGeneric(name = "wf", def = function(x) standardGeneric("wf"))
setMethod(f = "wf", signature = "SYSargs2", definition = function(x) {
  return(x@wf)
})
setGeneric(name = "clt", def = function(x) standardGeneric("clt"))
setMethod(f = "clt", signature = "SYSargs2", definition = function(x) {
  return(x@clt)
})
setGeneric(name = "yamlinput", def = function(x) standardGeneric("yamlinput"))
setMethod(f = "yamlinput", signature = "SYSargs2", definition = function(x) {
  return(x@yamlinput)
})
setGeneric(name = "cmdlist", def = function(x, ...) standardGeneric("cmdlist"))
setMethod(f = "cmdlist", signature = "SYSargs2", definition = function(x) {
  return(x@cmdlist)
})
setGeneric(name = "input", def = function(x) standardGeneric("input"))
setMethod(f = "input", signature = "SYSargs2", definition = function(x) {
  return(x@input)
})
setGeneric(name = "output", def = function(x) standardGeneric("output"))
setMethod(f = "output", signature = "SYSargs2", definition = function(x) {
  return(x@output)
})
setGeneric(name = "cwlfiles", def = function(x) standardGeneric("cwlfiles"))
setMethod(f = "cwlfiles", signature = "SYSargs2", definition = function(x) {
  return(x@cwlfiles)
})
setGeneric(name = "inputvars", def = function(x) standardGeneric("inputvars"))
setMethod(f = "inputvars", signature = "SYSargs2", definition = function(x) {
  return(x@inputvars)
})

setGeneric(name = "cmdToCwl", def = function(x) standardGeneric("cmdToCwl"))
setMethod(f = "cmdToCwl", signature = "SYSargs2", definition = function(x) {
  return(x@cmdToCwl)
})


## Constructor methods
## List to SYSargs2 with: as(mylist, "SYSargs2")
setAs(
  from = "list", to = "SYSargs2",
  def = function(from) {
    new("SYSargs2",
        targets = from$targets,
        targetsheader = from$targetsheader,
        modules = from$modules,
        wf = from$wf,
        clt = from$clt,
        yamlinput = from$yamlinput,
        cmdlist = from$cmdlist,
        input = from$input,
        output = from$output,
        cwlfiles = from$cwlfiles,
        inputvars = from$inputvars, 
        cmdToCwl = from$cmdToCwl
    )
})

## Coerce back to list: as(SYSargs2, "list")
setGeneric(name = "sysargs2", def = function(x) standardGeneric("sysargs2"))
setMethod(f = "sysargs2", signature = "SYSargs2", definition = function(x) {
  sysargs2 <- list(targets = x@targets, targetsheader = x@targetsheader, modules = x@modules, wf = x@wf,
                   clt = x@clt, yamlinput = x@yamlinput, cmdlist = x@cmdlist, input = x@input, output = x@output, 
                   cwlfiles = x@cwlfiles, inputvars = x@inputvars, cmdToCwl = x@cmdToCwl)
  return(sysargs2)
})

## SYSargs2 to list with: as(SYSargs2, "list")
setAs(from = "SYSargs2", to = "list",  def = function(from) {
    sysargs2(from)
})

## Define print behavior for SYSargs2
setMethod(
  f = "show", signature = "SYSargs2",
  definition = function(object) {
    cat(crayon::green$bold(paste0("Instance of '", class(object), "':")),
        paste0("   Slot names/accessors: "),
        paste0(
          "      targets: ", length(object@targets),
          " (", head(names(object@targets), 1), "...",
          tail(names(object@targets), 1), ")",
          ", targetsheader: ", length(unlist(object@targetsheader)), " (lines)"
        ),
        paste0("      modules: ", length(object@modules)),
        paste0(
          "      wf: ", length(object@wf),
          ", clt: ", length(object@clt),
          ", yamlinput: ", length(object@yamlinput), " (components)"
        ),
        paste0(
          "      input: ", length(object@input),
          ", output: ", length(object@output)
        ),
        paste0("      cmdlist: ", length(object@cmdlist)),
        "   WF Steps:",
        paste0(
          "      ", seq_along(object@clt), ". ", object@cwlfiles$steps,
          " (rendered: ", length(object@cmdlist[[1]]) != 0, ")"
        ),
        "\n",
        sep = "\n"
    )
})

## Extend names() method
setMethod(f = "names", signature = "SYSargs2",
  definition = function(x) {
    return(slotNames(x))
})

## Extend infile1() method
setMethod(f = "infile1", signature = "SYSargs2", definition = function(x) {
  subset_input <- input(x)
  subset_sample <- sapply(names(subset_input), function(x) list(NULL))
  for (i in seq_along(names(subset_input))) {
    if ("FileName" %in% names(subset_input[[i]])) {
      subset_sample[[i]] <- normalizePath(subset_input[[i]][["FileName"]])
    } else {
      subset_sample[[i]] <- normalizePath(subset_input[[i]][["FileName1"]])
    }
    subset_sample <- as.character(subset_sample)
    names(subset_sample) <- names(subset_input)
  }
  return(subset_sample)
})

## Extend infile2() method
setMethod(f = "infile2", signature = "SYSargs2", definition = function(x) {
  subset_input <- input(x)
  subset_sample <- sapply(names(subset_input), function(x) list(NULL))
  for (i in seq_along(names(subset_input))) {
    if ("FileName2" %in% names(subset_input[[i]])) {
      subset_sample[[i]] <- normalizePath(subset_input[[i]][["FileName2"]])
    }
    subset_sample <- as.character(subset_sample)
    names(subset_sample) <- names(subset_input)
  }
  return(subset_sample)
})

## Extend length() method
setMethod(
  f = "length", signature = "SYSargs2", definition = function(x) {
    return(length(x@targets))
})

## Convert targets data.frame to list
targets.as.list <- function(x) {
  targetslist <- yaml::yaml.load(yaml::as.yaml(x, column.major = FALSE))
  names(targetslist) <- x$SampleName
  return(targetslist)
}

## Usage:
# targets <- read.delim("targets.txt", comment.char = "#")
# targetslist <- targets.as.list(x=targets)

## Convert targets list to data.frame
targets.as.df <- function(x) {
  targetstmp <- sapply(x, as.character, simplify = FALSE)
  targetsDF <- as.data.frame(do.call("rbind", targetstmp))
  rownames(targetsDF) <- NULL
  colnames(targetsDF) <- names(x[[1]])
  return(S4Vectors::DataFrame(targetsDF))
}

## Usage:
# targets.as.df(x=targetslist)

## targets slot from a SYSargs2 obj to df with: as(SYSargs2, "data.frame")
setAs(from = "SYSargs2", to = "data.frame", def = function(from) {
    targets.as.df(targets(from))
})

# Behavior of "[" operator for SYSargs2
setMethod(f = "[", signature = "SYSargs2", definition = function(x, i, ..., drop) {
  if (is.logical(i)) {
    i <- which(i)
  }
  x@targets <- x@targets[i]
  x@input <- x@input[i]
  x@output <- x@output[i]
  x@cmdlist <- x@cmdlist[i]
  return(x)
})

## Behavior of "[[" operator for SYSargs2
setMethod(f = "[[", signature = c("SYSargs2", "ANY", "missing"),
  definition = function(x, i, ..., drop) {
    return(as(x, "list")[[i]])
})

## Behavior of "$" operator for SYSargs2
setMethod("$", signature = "SYSargs2",
          definition = function(x, name) {
            slot(x, name)
})

setGeneric(name = "baseCommand", def = function(x) standardGeneric("baseCommand"))
setMethod("baseCommand", signature = "SYSargs2", definition = function(x) {
  return(x@clt[[1]]$baseCommand)
})

## Replacement method for SYSargs2 using "[" operator
setReplaceMethod(f = "[[", signature = "SYSargs2", definition = function(x, i, j, value) {
  if (i == 1) x@targets <- value
  if (i == 2) x@targetsheader <- value
  if (i == 3) x@modules <- value
  if (i == 4) x@wf <- value
  if (i == 5) x@clt <- value
  if (i == 6) x@yamlinput <- value
  if (i == 7) x@cmdlist <- value
  if (i == 8) x@input <- value
  if (i == 9) x@output <- value
  if (i == 10) x@cwlfiles <- value
  if (i == "targets") x@targets <- value
  if (i == "targetsheader") x@targetsheader <- value
  if (i == "modules") x@modules <- value
  if (i == "wf") x@wf <- value
  if (i == "clt") x@clt <- value
  if (i == "yamlinput") x@yamlinput <- value
  if (i == "cmdlist") x@cmdlist <- value
  if (i == "input") x@input <- value
  if (i == "output") x@output <- value
  if (i == "cwlfiles") x@cwlfiles <- value
  if (i == "cmdToCwl") x@cmdToCwl <- value
  return(x)
})

## Replacement method
setGeneric(name="yamlinput<-", def=function(x, ..., value) standardGeneric("yamlinput<-"))
setReplaceMethod("yamlinput", c("SYSargs2"), function(x, ..., value) {
  x <- as(x, "list")
  x$yamlinput[[...]] <- value
  x <- as(x, "SYSargs2")
  x <- updateWF(x)
  x
})

setGeneric(name="cmdToCwl<-", def=function(x, ..., value) standardGeneric("cmdToCwl<-"))
setReplaceMethod("cmdToCwl", c("SYSargs2"), function(x,..., value) {
  x@cmdToCwl <- value
  x
})

####################################################
## Class and Method Definitions for SYSargsList ##
###################################################
## Define SYSargsList class
setClass("SYSargsList", slots = c(
  sprconfig = "list",
  stepsWF = "list",
  statusWF = "list",
  dependency = "list",
  projectWF = "list",
  targetsWF = "list",
  SEobj = "list",
  outfiles = "list",
  targets_connection="list"
))

## Methods to return SYSargsList components
setGeneric(name = "sprconfig", def = function(x) standardGeneric("sprconfig"))
setMethod(f = "sprconfig", signature = "SYSargsList", definition = function(x) {
  return(x@sprconfig)
})
setGeneric(name = "stepsWF", def = function(x) standardGeneric("stepsWF"))
setMethod(f = "stepsWF", signature = "SYSargsList", definition = function(x) {
  # if(length(x)==1){
  #   return(x@stepsWF[[1]])
  # } else {
    return(x@stepsWF)
  # }

})
setGeneric(name = "statusWF", def = function(x) standardGeneric("statusWF"))
setMethod(f = "statusWF", signature = "SYSargsList", definition = function(x) {
  return(x@statusWF)
})
setGeneric(name = "dependency", def = function(x) standardGeneric("dependency"))
setMethod(f = "dependency", signature = "SYSargsList", definition = function(x) {
  return(x@dependency)
})
setGeneric(name = "projectWF", def = function(x) standardGeneric("projectWF"))
setMethod(f = "projectWF", signature = "SYSargsList", definition = function(x) {
  return(x@projectWF)
})
setGeneric(name = "targetsWF", def = function(x) standardGeneric("targetsWF"))
setMethod(f = "targetsWF", signature = "SYSargsList", definition = function(x) {
  return(x@targetsWF)
})

setGeneric(name = "SEobj", def = function(x) standardGeneric("SEobj"))
setMethod(f = "SEobj", signature = "SYSargsList", definition = function(x) {
  return(x@SEobj)
})

setGeneric(name = "outfiles", def = function(x) standardGeneric("outfiles"))
setMethod(f = "outfiles", signature = "SYSargsList", definition = function(x) {
  return(x@outfiles)
})


## Coerce back to list: as(SYSargsList, "list")
setGeneric(name = "sysargslist", def = function(x) standardGeneric("sysargslist"))
setMethod(f = "sysargslist", signature = "SYSargsList", definition = function(x) {
  sysargslist <- list(sprconfig = x@sprconfig, stepsWF = x@stepsWF,
                      statusWF = x@statusWF, dependency = x@dependency,
                      projectWF = x@projectWF, targetsWF = x@targetsWF, SEobj = x@SEobj,
                      outfiles = x@outfiles, targets_connection = x@targets_connection
  )
  return(sysargslist)
})

## Constructor methods
## List to SYSargsList
setAs(from = "list", to = "SYSargsList",
      def = function(from) {
        new("SYSargsList",
            sprconfig = from$sprconfig, stepsWF = from$stepsWF,
            statusWF = from$statusWF, dependency = from$dependency, projectWF = from$projectWF,
            targetsWF = from$targetsWF, SEobj = from$SEobj, outfiles = from$outfiles,
            targets_connection = from$targets_connection
        )
      })

## SYSargsList to list with: as("SYSargsList", list)
setAs(from = "SYSargsList", to = "list",
      def = function(from) {
        sysargslist(from)
      })


## Define print behavior for SYSargsList
setMethod(f = "show", signature = "SYSargsList",
          definition = function(object) {
            if(length(object)>0){
              status_1 <- as.character()
              for(i in seq_along(object@stepsWF)){
                status_1 <- c(status_1, paste0(
                  "      ", i, ". ", names(object@stepsWF)[i],
                  " (Status: ", object@statusWF[[i]]$status, ")", sep="\n"))
              }
              status <- append("   WF Steps:\n", status_1)
            } else {
              status <- "No workflow steps added"
            }
            cat(crayon::blue$bold(paste0("Instance of '", class(object), "':")), "\n",
                status, "\n")
          })

## Extend names() method
setMethod(f = "names", signature = "SYSargsList",
          definition = function(x) {
            return(slotNames(x))
          })

## Extend length() method
setMethod(f = "length", signature = "SYSargsList",
          definition = function(x) {
            return(length(x@stepsWF))
          })

# Behavior of "[" operator for SYSargsList
setMethod(f = "[", signature = "SYSargsList", definition = function(x, i, ..., drop) {
  if(missing(i)){
    i <- 1:length(x)
  }
    if (is.logical(i)) {
    i <- which(i)
  }
  for(s in seq_along(i)){
    tryCatch({
      if(i[s]<0) { ii <- i[s]*-1 } else {ii=i[s]}
      x@stepsWF[[ii]]
    }, error=function(e) {
      e$message <- paste0("ERROR: ", '\n',
                          "Subset out of bounds. Please subset accordingly with the 'length(x)'.")
      stop(e)
    })
  }
  x@stepsWF <- x@stepsWF[i]
  x@statusWF <- x@statusWF[i]
  x@dependency <- x@dependency[i]
  x@targetsWF <- x@targetsWF[i]
  x@SEobj <- x@SEobj[i]
  x@outfiles <- x@outfiles[i]
  return(x)
})

## Behavior of "subsetSamples" operator for SYSargsList
setGeneric(name = "subsetSamples", def = function(x, samples) standardGeneric("subsetSamples"))
setMethod(f = "subsetSamples", signature = "SYSargsList", definition = function(x, samples) {
  if(missing(samples)){
    samples <- 1:max(sapply(stepsWF(x), function(x) length(x)))
  }
  for(s in seq_along(x)){
    x@stepsWF[[s]] <- x@stepsWF[[s]][samples]
    x@targetsWF[[s]] <- x@targetsWF[[s]][samples,]
    #x@SEobj[[s]] <- x@SEobj[[s]][i]
    x@outfiles[[s]] <- x@outfiles[[s]][samples,]
  }
  
  x
})

# ## Behavior of "[[" operator for SYSargsList
# setMethod(f = "[[", signature = "SYSargsList", definition = function(x, i, ..., drop) {
#   #setMethod(f = "[[", signature = c("SYSargsList","ANY", "missing"), definition = function(x, i, ..., drop) {
#   for(s in seq_along(x)){
#     x@stepsWF[[s]] <- x@stepsWF[[s]][i]
#     x@targetsWF[[s]] <- x@targetsWF[[s]][i,]
#     #x@SEobj[[s]] <- x@SEobj[[s]][i]
#     x@outfiles[[s]] <- x@outfiles[[s]][i,]
#   }
# 
#   x
# })

## Behavior of "[[" operator for SYSargsList
setMethod(f = "[[", signature = c("SYSargsList", "ANY", "missing"), definition = function(x, i, ..., drop) {
            return(as(x, "list")[[i]])
})



## Behavior of "$" operator for SYSargsList
setMethod("$", signature = "SYSargsList",
          definition = function(x, name) {
            slot(x, name)
          })


setMethod(f = "cmdlist", signature = "SYSargsList", definition = function(x, sample=NULL) {
  cmd <- sapply(names(x$stepsWF), function(x) list(NULL))
  for(i in seq_along(x)){
    if(nchar(cmdlist(x$stepsWF[[i]])[[1]][[1]])>0){
      cmd_list <- cmdlist(x$stepsWF[[i]])
      if(!is.null(sample)){
        cmd_list <- cmd_list[sample]
      }
      cmd[[i]] <- cmd_list
    }
  }
  return(cmd)
})

setMethod(f = "yamlinput", signature = "SYSargsList", definition = function(x) {
  if(length(x) > 1) stop("`x` needs to have length(x)==1")
  if(inherits(stepsWF(x), "LineWise")) stop("Provide a stepWF with a 'SYSargs2' class")
  stepsWF(x)[[1]]$yamlinput
})




# setGeneric(name="yamlinput<-", def=function(x, ..., value) standardGeneric("yamlinput<-"))
setReplaceMethod("yamlinput", c("SYSargsList"), function(x, paramName, value) {
  x <- x
  args <- x@stepsWF[[1]]
  yamlinput(args, paramName) <- value
  x <- sysargslist(x)
  x$stepsWF[[1]] <- args
  x <- as(x, "SYSargsList")
  x
})

## Replacement method for SYSargsList using "[[" operator
setReplaceMethod(f = "[[", signature = "SYSargsList", definition = function(x, i, j, value) {
  if (i == 1) x@sprconfig <- value
  if (i == 2) x@stepsWF <- value
  if (i == 3) x@statusWF <- value
  if (i == 4) x@dependency <- value
  if (i == 5) x@projectWF <- value
  if (i == 6) x@targetsWF <- value
  if (i == 7) x@SEobj <- value
  if (i == "sprconfig") x@sprconfig <- value
  if (i == "stepsWF") x@stepsWF <- value
  if (i == "statusWF") x@statusWF <- value  
  if (i == "dependency") x@dependency <- value 
  if (i == "projectWF") x@projectWF <- value
  if (i == "targetsWF") x@targetsWF <- value
  if (i == "SEobj") x@SEobj <- value
  return(x)
})

## Replacement method

setGeneric(name="sprconfig<-", def=function(x, ..., value) standardGeneric("sprconfig<-"))
setReplaceMethod("sprconfig", c("SYSargsList"), function(x,..., value) {
  x@sprconfig <- value
  x
})

setGeneric(name="appendStep<-", def=function(x, after=length(x), ..., value) standardGeneric("appendStep<-"))
setReplaceMethod("appendStep", c("SYSargsList"), function(x, after=length(x), ..., value) {
  lengx <- length(x)
  after <- after
  if(inherits(value, "SYSargsList")){
    value <- .validationStepConn(x, value)
    x <- sysargslist(x)
    if(names(value$stepsWF)=="Step_x"){
      step_name <- paste0("Step_", after +1L)
      renameStep(value, 1) <- step_name
    }
    if (!after) {
      x$stepsWF <- c(value$stepsWF, x$stepsWF)
      x$targetsWF <- c(targetsWF(value), x$targetsWF)
      x$statusWF <- c(statusWF(value), x$statusWF)
      x$dependency <- c(dependency(value), x$dependency)
      x$outfiles <- c(outfiles(value), x$outfiles)
    } else if (after >= lengx) {
      x$stepsWF <- c(x$stepsWF, value$stepsWF)
      x$targetsWF <- c(x$targetsWF, targetsWF(value))
      x$statusWF <- c(x$statusWF, statusWF(value))
      x$dependency <- c(x$dependency, dependency(value))
      x$outfiles <- c(x$outfiles, outfiles(value))
    } else {
      x$stepsWF <- c(x$stepsWF[1L:after], value$stepsWF, x$stepsWF[(after + 1L):lengx])
      x$targetsWF <- c(x$targetsWF[1L:after], targetsWF(value), x$targetsWF[(after + 1L):lengx])
      x$statusWF <- c(x$statusWF[1L:after], statusWF(value), x$statusWF[(after + 1L):lengx])
      x$dependency <- c(x$dependency[1L:after], dependency(value), x$dependency[(after + 1L):lengx])
      x$outfiles <- c(x$outfiles[1L:after], outfiles(value), x$outfiles[(after + 1L):lengx])
    }
    x <- as(x, "SYSargsList")
  } else if(inherits(value, "LineWise")){
    if(value$stepName=="Step_x"){
      step_name <- paste0("Step_", after +1L)
    } else {
      step_name <- value$stepName
    }
    x <- sysargslist(x)
    if (!after) {
      x$stepsWF <- c(value, x$stepsWF)
      x$targetsWF <- c(DataFrame(), x$targetsWF)
      x$statusWF <- c(list(list(status="Pending")), x$statusWF)
      x$dependency <- c(NA, x$dependency)
      x$outfiles <- c(DataFrame(), x$outfiles)
    } else if (after >= lengx) {
      x$stepsWF <- c(x$stepsWF, value)
      x$targetsWF <- c(x$targetsWF, DataFrame())
      x$statusWF <- c(x$statusWF, list(list(status="Pending")))
      x$dependency <- c(x$dependency, NA)
      x$outfiles <- c(x$outfiles, DataFrame())
    } else {
      x$stepsWF <- c(x$stepsWF[1L:after], value, x$stepsWF[(after + 1L):lengx])
      x$targetsWF <- c(x$targetsWF[1L:after], DataFrame() , x$targetsWF[(after + 1L):lengx])
      x$statusWF <- c(x$statusWF[1L:after], list(list(status="Pending")), x$statusWF[(after + 1L):lengx])
      x$dependency <- c(x$dependency[1L:after], NA, x$dependency[(after + 1L):lengx])
      x$outfiles <- c(x$outfiles[1L:after], DataFrame(), x$outfiles[(after + 1L):lengx])
    }
    names(x$stepsWF)[after+1L] <- step_name
    names(x$statusWF)[after+1L] <- step_name
    names(x$dependency)[after+1L] <- step_name
    names(x$targetsWF)[after+1L] <- step_name
    names(x$outfiles)[after+1L] <- step_name
    x <- as(x, "SYSargsList")
  } else stop("Argument 'value' needs to be assigned an object of class 'SYSargsList' OR 'LineWise'.")
  if(any(duplicated(names(stepsWF(x))))) warning("Duplication is found in names(stepsWF(x)). Consider renaming the steps.")
  x
})

.validationStepConn <- function(x, value){
  if (length(value$targets_connection) > 0){
    step <- value$targets_connection[[1]]
    targets_name <- paste(colnames(targetsWF(x)[step][[1]]),collapse="|")
    new_targets_col <- value$targets_connection[[2]][-c(which(grepl(targets_name, value$targets_connection[[2]])))]
    if(!step %in% names(stepsWF(x))) stop(paste0("'targets' argument needs to be assigned as valid targets file OR the names of a previous step, for example: ", "\n",
                                                 paste0(names(stepsWF(x)), collapse = " OR ")))
    if(all(!new_targets_col %in% colnames(x$outfiles[[1]]))) stop(paste0("'targets_column' argument needs to be assigned as valid column names of a previous step, for example: ", "\n",
                                                                         paste0(colnames(x$outfiles[[1]]), collapse = " OR \n")))
    if(is.null(value$targets_connection$rm_targets_col)){
      old_targets <- x$targetsWF[[step]]
    } else {
      old_targets <- x$targetsWF[[step]][-c(which(grepl(paste(value$targets_connection[[3]], collapse="|"), colnames(x$targetsWF[[step]]))))]
    }
    new_targets <- cbind(x$outfiles[[1]][new_targets_col], old_targets)
    WF <- value$stepsWF[[1]]
    #inputvars_v <- unlist(WF$inputvars)
    ## TODO: check inputvars...
    WF2 <- updateWF(WF, new_targets= targets.as.list(data.frame(new_targets)), inputvars=WF$inputvars, write.yaml = FALSE)
    value <- sysargslist(value)
    value$stepsWF[[1]] <- WF2
    value$targetsWF[[1]] <- as(WF2, "data.frame")
    value$outfiles <- output(WF2)
    value <- as(value, "SYSargsList")
  }
  if(all(!is.na(dependency(value)))){
    dep <- dependency(value)[[1]]
    if(inherits(dep, "character")){
      if(all(!dep %in% names(stepsWF(x)))) stop(
        "'dependency' argument needs to be assigned as valid previous Step Name, for example: ", "\n",
        paste0(names(stepsWF(x)), collapse = " OR "))
    } else {
      if(inherits(dep, "numeric")){
        if(all(!dep %in% 1:length(stepsWF(x)))) stop(
          "'dependency' argument needs to be assigned as valid previous Step Index, for example: ", "\n",
          paste0(1:length(stepsWF(x)), collapse = " OR "))
      }
    }
  }
  return(value)
}

## USage:
# appendStep(sal) <- SYSargsList(WF)
# appendStep(sal, after=0) <- SYSargsList(WF)
# sal
# appendStep(sal, after=0, step_index="test_11") <- SYSargsList(WF)

setGeneric(name="replaceStep<-", def=function(x, step, step_name="default", value) standardGeneric("replaceStep<-"))
setReplaceMethod("replaceStep", c("SYSargsList"), function(x, step, step_name="default", value) {
  if(!is(value, "SYSargsList")) stop("Argument 'value' needs to be assigned an object of class 'SYSargsList'.")
  if(length(value) > 1) stop("Argument 'value' cannot have 'length(value) > 1")
  if(is(step, "numeric")){
    if(step > length(x)) stop(paste0("Argument 'step' cannot be greater than ", length(x)))
  } else if(is(step, "character")){
    if(!step %in% names(stepsWF(x))) stop(paste0("Argument 'step' needs to be assigned one of the following: ",
                                                 paste(names(stepsWF(x)), collapse = " OR ")))
  }
  x <- sysargslist(x)
  x$stepsWF[step] <- value$stepsWF
  if(step_name=="default"){
    name <- names(value$stepsWF)
    if(name %in% names(x$stepsWF)){
      names(x$stepsWF)[step] <- paste0("Step_", step)
      names(x$statusWF)[step] <-paste0("Step_", step)
      names(x$dependency)[step] <- paste0("Step_", step)
      names(x$targetsWF)[step] <- paste0("Step_", step)
      names(x$outfiles)[step] <- paste0("Step_", step)
      cat(paste0("Index name of x", "[", step, "]", " was rename to ", paste0("Step_", step), " to avoid duplications."))
    } else {
      names(x$stepsWF)[step] <- names(value$stepsWF)
      names(x$statusWF)[step] <- names(value$stepsWF)
      names(x$dependency)[step] <- names(value$stepsWF)
      names(x$targetsWF)[step] <- names(value$stepsWF)
      names(x$outfiles)[step] <- names(value$stepsWF)
    }
  } else {
    names(x$stepsWF)[step] <- step_name
    names(x$statusWF)[step] <- step_name
    names(x$dependency)[step] <- step_name
    names(x$targetsWF)[step] <- step_name
    names(x$outfiles)[step] <- step_name
  }
  x <- as(x, "SYSargsList")
  if(any(duplicated(names(stepsWF(x))))) warning("Duplication is found in names(stepsWF(x)). Consider renaming the steps.")
  x
})

# replaceStep(sal, 1) <- sal[1]
# sal
# replaceStep(sal, 3, step_name="test") <- sal[1]

setGeneric(name="stepsWF<-", def=function(x, step, ..., value) standardGeneric("stepsWF<-"))
setReplaceMethod("stepsWF", c("SYSargsList"), function(x, step, ..., value) {
    x@stepsWF[[step]] <- value
  x
})

setGeneric(name="renameStep<-", def=function(x, step, ..., value) standardGeneric("renameStep<-"))
setReplaceMethod("renameStep", c("SYSargsList"), function(x, step, ..., value) {
  if(length(step)!=length(value)) stop("value argument needs to be the same length of the step for rename")
  if(is(value, "character")){
    names(x@stepsWF)[step] <- value
    names(x@statusWF)[step] <- value
    names(x@dependency)[step] <- value
    names(x@targetsWF)[step] <- value
    names(x@outfiles)[step] <- value
    #names(x@SEobj)[step] <- value
    #names(x@targets_connection)[step] <- value
    
  }  else {
    stop("Replace value needs to be assigned an 'character' name for the workflow step.")
  }
  x
})

setGeneric(name="statusWF<-", def=function(x, step, ..., value) standardGeneric("statusWF<-"))
setReplaceMethod("statusWF", c("SYSargsList"), function(x, step, ..., value) {
  if(is(value, "character")){
    x@statusWF[step] <- value
  }
  x
})
