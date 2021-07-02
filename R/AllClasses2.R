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
  files = "list",
  inputvars = "list", 
  cmdToCwl = "list", 
  status = "list", 
  internal_outfiles = "list" 
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
setGeneric(name = "yamlinput", def = function(x, ...) standardGeneric("yamlinput"))
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
setGeneric(name = "files", def = function(x) standardGeneric("files"))
setMethod(f = "files", signature = "SYSargs2", definition = function(x) {
  return(x@files)
})
setGeneric(name = "inputvars", def = function(x) standardGeneric("inputvars"))
setMethod(f = "inputvars", signature = "SYSargs2", definition = function(x) {
  return(x@inputvars)
})

setGeneric(name = "cmdToCwl", def = function(x) standardGeneric("cmdToCwl"))
setMethod(f = "cmdToCwl", signature = "SYSargs2", definition = function(x) {
  return(x@cmdToCwl)
})

setGeneric(name = "status", def = function(x) standardGeneric("status"))
setMethod(f = "status", signature = "SYSargs2", definition = function(x) {
  return(x@status)
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
        files = from$files,
        inputvars = from$inputvars, 
        cmdToCwl = from$cmdToCwl, 
        status = from$status, 
        internal_outfiles = from$internal_outfiles
        
    )
})

## Coerce back to list: as(SYSargs2, "list")
setGeneric(name = "sysargs2", def = function(x) standardGeneric("sysargs2"))
setMethod(f = "sysargs2", signature = "SYSargs2", definition = function(x) {
  sysargs2 <- list(targets = x@targets, targetsheader = x@targetsheader, modules = x@modules, wf = x@wf,
                   clt = x@clt, yamlinput = x@yamlinput, cmdlist = x@cmdlist, input = x@input, output = x@output, 
                   files = x@files, inputvars = x@inputvars, cmdToCwl = x@cmdToCwl,
                   status = x@status, internal_outfiles = x@internal_outfiles)
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
          "      wf: ", length(object@wf$steps), 
          ", clt: ", length(object@clt),
          ", yamlinput: ", length(object@yamlinput), " (inputs)"
        ),
        paste0(
          "      input: ", length(object@input),
          ", output: ", length(object@output)
        ),
        paste0("      cmdlist: ", length(object@cmdlist)),
        "   Sub Steps:",
        paste0(
          "      ", seq_along(object@clt), ". ", object@files$steps,
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
    return(length(x@cmdlist))
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
  return(targetsDF)
}

## Usage:
# targets.as.df(x=targetslist)

## targets slot from a SYSargs2 obj to df with: as(SYSargs2, "data.frame")
setAs(from = "SYSargs2", to = "DataFrame", def = function(from) {
  S4Vectors::DataFrame(targets.as.df(targets(from)))
})

# Behavior of "[" operator for SYSargs2
setMethod(f = "[", signature = "SYSargs2", definition = function(x, i, ..., drop) {
  if (is.logical(i)) {
    i <- which(i)
  }
  x@targets <- x@targets[i]
  x@input <- x@input[i]
  x@output <- x@output[i]
  x@internal_outfiles <- x@internal_outfiles[i]
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
  return(x@clt[[1]]$baseCommand[[1]])
})

setMethod("SampleName", signature = "SYSargs2", definition = function(x) {
  targets_x <- targets(x)
  if(length(targets_x) > 0){
    sample_name_x <- as(x, "DataFrame")
    return(sample_name_x$SampleName)
  } else if(length(targets_x)==0){
    message("This step doesn't contain multiple samples.")
  }
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
  if (i == 10) x@files <- value
  if (i == 11) x@status <- value
  if (i == 12) x@internal_outfiles <- value
  if (i == "targets") x@targets <- value
  if (i == "targetsheader") x@targetsheader <- value
  if (i == "modules") x@modules <- value
  if (i == "wf") x@wf <- value
  if (i == "clt") x@clt <- value
  if (i == "yamlinput") x@yamlinput <- value
  if (i == "cmdlist") x@cmdlist <- value
  if (i == "input") x@input <- value
  if (i == "output") x@output <- value
  if (i == "files") x@files <- value
  if (i == "cmdToCwl") x@cmdToCwl <- value
  if (i == "status") x@status <- value
  if (i == "internal_outfiles") x@internal_outfiles <- value 
    return(x)
})

## Replacement method
setGeneric(name="yamlinput<-", def=function(x, paramName, ..., value) standardGeneric("yamlinput<-"))
setReplaceMethod("yamlinput", c("SYSargs2"), function(x, paramName, value) {
  x <- as(x, "list")
  ## Check paramName
  if(!paramName %in% names(x$yamlinput)) stop ("'paramName' argument need to be one of following")
  ## Check class of value
  if(!identical(class(x$yamlinput[[paramName]]), class(value))) stop("message")
  x$yamlinput[[paramName]] <- value
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
  stepsWF = "list",
  statusWF = "list",
  targetsWF = "list",
  outfiles = "list",
  SEobj = "list",
  dependency = "list",
  targets_connection="list",
  projectInfo = "list",
  runInfo = "list"
))

## Methods to return SYSargsList components
setGeneric(name = "stepsWF", def = function(x) standardGeneric("stepsWF"))
setMethod(f = "stepsWF", signature = "SYSargsList", definition = function(x) {
    return(x@stepsWF)
})
setGeneric(name = "statusWF", def = function(x) standardGeneric("statusWF"))
setMethod(f = "statusWF", signature = "SYSargsList", definition = function(x) {
  return(lapply(lapply(x@statusWF, "[[", 2), function(y) S4Vectors::DataFrame(y, check.names=FALSE)))
  #return(x@statusWF)
})
setGeneric(name = "targetsWF", def = function(x) standardGeneric("targetsWF"))
setMethod(f = "targetsWF", signature = "SYSargsList", definition = function(x) {
  return(x@targetsWF)
})
setGeneric(name = "outfiles", def = function(x) standardGeneric("outfiles"))
setMethod(f = "outfiles", signature = "SYSargsList", definition = function(x) {
  return(x@outfiles)
})
setGeneric(name = "SEobj", def = function(x) standardGeneric("SEobj"))
setMethod(f = "SEobj", signature = "SYSargsList", definition = function(x) {
  return(x@SEobj)
})
setGeneric(name = "dependency", def = function(x) standardGeneric("dependency"))
setMethod(f = "dependency", signature = "SYSargsList", definition = function(x) {
  return(x@dependency)
})
setGeneric(name = "projectInfo", def = function(x) standardGeneric("projectInfo"))
setMethod(f = "projectInfo", signature = "SYSargsList", definition = function(x) {
  return(x@projectInfo)
})
setGeneric(name = "runInfo", def = function(x) standardGeneric("runInfo"))
setMethod(f = "runInfo", signature = "SYSargsList", definition = function(x) {
  return(x@runInfo)
})

## Coerce back to list: as(SYSargsList, "list")
setGeneric(name = "sysargslist", def = function(x) standardGeneric("sysargslist"))
setMethod(f = "sysargslist", signature = "SYSargsList", definition = function(x) {
  sysargslist <- list(stepsWF = x@stepsWF,
                      statusWF = x@statusWF, 
                      targetsWF = x@targetsWF,
                      outfiles = x@outfiles,
                      SEobj = x@SEobj,
                      dependency = x@dependency,
                      targets_connection = x@targets_connection,
                      projectInfo = x@projectInfo,
                      runInfo = x@runInfo
                      )
  return(sysargslist)
})

## Constructor methods
## List to SYSargsList
setAs(from = "list", to = "SYSargsList",
      def = function(from) {
        new("SYSargsList",
            stepsWF = from$stepsWF,
            statusWF = from$statusWF, 
            targetsWF = from$targetsWF,
            outfiles = from$outfiles,
            SEobj = from$SEobj,
            dependency = from$dependency, 
            targets_connection = from$targets_connection,
            projectInfo = from$projectInfo,
            runInfo = from$runInfo
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
              status_1 <- .showHelper(object)
              status <- append("   WF Steps:\n", status_1)
            } else {
              status <- "No workflow steps added"
            }
            cat(crayon::blue$bold(paste0("Instance of '", class(object), "':")), "\n",
                status, "\n")
          })

.showHelper <- function(object){
  status_1 <- as.character()
  for (i in seq_along(object@stepsWF)) {
    status_color <- switch(
      tolower(object@statusWF[[i]]$status.summary),
      "pending" = crayon::blue$bold,
      "warning" = crayon::make_style("orange")$bold,
      "error" = crayon::red$bold,
      "success" = crayon::green$bold
    )
    if(inherits(object@stepsWF[[i]], "SYSargs2")){
      status_1 <- c(
        status_1,
        c(
          paste0("      ", i, ". ", crayon::blue(names(object@stepsWF)[i], "--> Status: "), status_color(object@statusWF[[i]]$status.summary), sep = ""), "\n",
          paste0(c("          Total Files: ", "Existing: ", "Missing: "), colSums(object@statusWF[[i]][[2]][2:4]), collapse = " | "), "\n",
          paste0(
            # "         Sub Steps:", "\n",
            paste0(
              "        ", i, ".", seq_along(object@stepsWF[[i]]@clt), ". ", crayon::green(object@stepsWF[[i]]$files$steps)
            ), "\n",
            paste0("             cmdlist: ", length(object@stepsWF[[i]]), " | "),
            sapply(as.list(object@stepsWF[[i]]$files$steps), function(x) {
              paste0(
                paste0(names(table(unlist(object@statusWF[[i]][[2]][object@stepsWF[[i]]$files$steps][x]))), ": ",
                       table(unlist(object@statusWF[[i]][[2]][object@stepsWF[[i]]$files$steps][x])),
                       collapse = " | "
                )
              )
            }),
            "\n"
          )
        )
      )
    } else if(inherits(object@stepsWF[[i]], "LineWise")){
      status_1 <- c(status_1,
                    paste0("      ", i, ". ", crayon::blue(names(object@stepsWF)[i], "--> Status: "), status_color(object@statusWF[[i]]$status.summary), "\n"))
    }
  }
  return(status_1)
}

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
      e$message <- paste0('\n',
                          "Step number is out of range. Please subset accordingly with the 'length(x)'", 
                          '\n',
                          paste0(1:length(x), collapse = ", "))
      stop(e)
    })
  }
  x@stepsWF <- x@stepsWF[i]
  names_tc <- names(x@stepsWF)
  x@statusWF <- x@statusWF[i]
  x@targetsWF <- x@targetsWF[i]
  x@outfiles <- x@outfiles[i]
  x@SEobj <- x@SEobj[i]
  x@dependency <- x@dependency[i]
  x@targets_connection <- x@targets_connection[names(x@targets_connection) %in% names_tc]
  x@projectInfo <- x@projectInfo
  x@runInfo$directory <- x@runInfo$directory[i]
  x <- .check_write_SYSargsList(x)
  return(x)
})

## Behavior of "subsetTargets" method for SYSargsList
setGeneric(name = "subsetTargets", def = function(x, subset_steps, input_targets, keep_steps=TRUE) standardGeneric("subsetTargets"))
setMethod(f = "subsetTargets", signature = "SYSargsList", definition = function(x, subset_steps, input_targets,  keep_steps=TRUE) {
  x_sub <- x[subset_steps]
  ## check subset_steps length
  if(length(unique(sapply(stepsWF(x_sub), function(x) length(x)))) > 1) stop("All 'subset_steps' should contain the same length.")
  if(missing(input_targets)){
    input_targets <- 1:max(sapply(stepsWF(x_sub), function(x) length(x)))
  } 
  # Check targets index, names
  if(inherits(input_targets, "numeric")){
    if(!all(input_targets %in% sapply(stepsWF(x_sub), function(x) 1:length(x)))) stop("Please select the number of 'input_targets' accordingly, options are: ",
                                                                          paste0(1:length(x_sub@stepsWF[[1]]), collapse=", "))
  } else if(inherits(input_targets, "character")){
    if(!all(input_targets %in% sapply(stepsWF(x_sub), function(x) SampleName(x)))) stop("Please select the number of 'input_targets' accordingly, options are: ",
                                                                                        paste0(SampleName(x_sub@stepsWF[[1]]), collapse=", "))
    input_targets <- which(SampleName(x_sub$stepsWF[[1]]) %in% input_targets)
  }
  
  if(keep_steps==FALSE){
    x <- x_sub
    subset_steps <- 1:length(x_sub)
  }
  for(s in subset_steps){
    x@stepsWF[[s]] <- x@stepsWF[[s]][input_targets]
    x@statusWF[[s]]$status.completed <- x@statusWF[[s]]$status.completed[input_targets,]
    x@statusWF[[s]]$status.time <- x@statusWF[[s]]$status.time[input_targets,]
    x@targetsWF[[s]] <- x@targetsWF[[s]][input_targets,]
    out <- DataFrame(x@outfiles[[s]][input_targets,])
    colnames(out) <- colnames(x@outfiles[[s]])
    x@outfiles[[s]] <- out
    #x_sub@SEobj[[s]] <- x_sub@SEobj[[s]][i]
    x@dependency <- x@dependency
    x@targets_connection <- x@targets_connection
    x@projectInfo <- x@projectInfo
    x@runInfo$directory <- x@runInfo$directory
  }
  x <- .check_write_SYSargsList(x)
  x
})

setMethod("SampleName", signature = "SYSargsList", definition = function(x, step) {
  ## Check steps
  if(inherits(step, "numeric")){
    if(!step %in% 1:length(x)) stop("We can not find this step in the Workflow")
  } else if(inherits(step, "character")){
    if(!step %in% stepName(x)) stop("We can not find this step in the Workflow")
  }
  if(!is.null(targetsWF(x)[[step]]$SampleName)){
    return(targetsWF(x)[[step]]$SampleName)
  } else if(is.null(targetsWF(x)[[step]]$SampleName)){
    message("This step doesn't contain multiple samples.")
  }
})

setGeneric(name = "stepName", def = function(x) standardGeneric("stepName"))
setMethod("stepName", signature = "SYSargsList", definition = function(x) {
    return(names(stepsWF(x)))
})


setMethod("targetsheader", signature = "SYSargsList", definition = function(x, step) {
  return(stepsWF(x)[[step]]$targetsheader)
})

setGeneric(name = "getColumn", def = function(x, step, df, column=1, names=SampleName(x, step)) standardGeneric("getColumn"))
setMethod("getColumn", signature = "SYSargsList", definition = function(x, step, df=c("targetsWF", "outfiles"), column=1, names=SampleName(x, step)) {
  ## assertions 
  stopifnot(inherits(x, "SYSargsList"))
  stopifnot(length(step) == 1)
  stopifnot(length(column) == 1)
  ## Check steps
  if(inherits(step, "numeric")){
    if(!step %in% 1:length(x)) stop("We can not find this step in the Workflow")
  } else if(inherits(step, "character")){
    if(!step %in% stepName(x)) stop("We can not find this step in the Workflow")
  }
  ## Check column
  if(inherits(column, "numeric")){
    if(!column %in% 1:ncol(x[[df]][[step]])) stop("We can not find this column in the Workflow")
  } else if(inherits(column, "character")){
    if(!column %in% colnames(x[[df]][[step]])) stop("We can not find this column in the Workflow")
  }
  ## Check names
  if(!length(names) ==  length(x[[df]][[step]][[column]])) stop("'names' argument needs to have the same length of desired output")
  ## 
  if(!is.null(x[[df]][[step]][[column]])){
    subset <- x[[df]][[step]][[column]]
    names(subset) <- names
  } else {
    message("This step doesn't contain expected outfiles.")
  }
  return(subset)
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

## viewEnvir() methods for SYSargslist
setGeneric(name = "viewEnvir", def = function(x) standardGeneric("viewEnvir"))
setMethod(f = "viewEnvir", signature = "SYSargsList", definition = function(x) {
  print(x@runInfo$env)
  print(ls(x@runInfo$env, all.names = TRUE))
})

## copyEnvir() methods for SYSargslist
setGeneric(name = "copyEnvir", def = function(x, list=character(), new.env=globalenv()) standardGeneric("copyEnvir"))
setMethod(f = "copyEnvir", signature = "SYSargsList", definition = function(x, list=character(), new.env) {
  envir <- x@runInfo$envir
  print(envir)
  if(length(list)==0){
    list <- ls(envir, all.names=TRUE)
  } else {
    list <- list
  }
  for(l in list) {
    assign(l, get(l, envir), new.env)
  }
  cat(paste0("Copying to 'new.env': ", "\n", paste0(list, collapse = ", ")))
})

## cmdlist method for SYSargslist
setMethod(f = "cmdlist", signature = "SYSargsList", definition = function(x, targets=NULL) {
  cmd <- sapply(names(x$stepsWF), function(x) list(NULL))
  for(i in seq_along(x)){
    if(nchar(cmdlist(x$stepsWF[[i]])[[1]][[1]])>0){
      cmd_list <- cmdlist(x$stepsWF[[i]])
      if(!is.null(targets)){
        cmd_list <- cmd_list[targets]
      }
      cmd[[i]] <- cmd_list
    }
  }
  return(cmd)
})

setMethod(f = "yamlinput", signature = "SYSargsList", definition = function(x, step) {
  if(inherits(stepsWF(x)[[step]], "LineWise")) stop("Provide a stepWF with a 'SYSargs2' class")
  stepsWF(x)[[step]]$yamlinput
})

# setGeneric(name="yamlinput<-", def=function(x, paramName, value, ...) standardGeneric("yamlinput<-"))
setReplaceMethod("yamlinput", c("SYSargsList"), function(x, step, paramName, value) {
  x_sub <- x[step]
  args <- x_sub@stepsWF[[1]]
  yamlinput(args, paramName) <- value
  x <- sysargslist(x)
  x$stepsWF[[1]] <- args
  x <- as(x, "SYSargsList")
  x <- .check_write_SYSargsList(x)
  x
})

## Replacement method for SYSargsList using "[[" operator
setReplaceMethod(f = "[[", signature = "SYSargsList", definition = function(x, i, j, value) {
  if (i == 1) x@stepsWF <- value
  if (i == 2) x@statusWF <- value
  if (i == 3) x@targetsWF <- value
  if (i == 4) x@outfiles <- value
  if (i == 5) x@SEobj <- value
  if (i == 6) x@dependency <- value
  if (i == 7) x@targets_connection <- value
  if (i == 8) x@projectInfo <- value
  if (i == 9) x@runInfo <- value
  if (i == "stepsWF") x@stepsWF <- value
  if (i == "statusWF") x@statusWF <- value  
  if (i == "targetsWF") x@targetsWF <- value
  if (i == "outfiles") x@outfiles <- value
  if (i == "SEobj") x@SEobj <- value
  if (i == "dependency") x@dependency <- value 
  if (i == "targets_connection") x@targets_connection <- value
  if (i == "projectInfo") x@projectInfo <- value
  if (i == "runInfo") x@runInfo <- value
  return(x)
})

## Replacement method

setGeneric(name="appendStep<-", def=function(x, after=length(x), ..., value) standardGeneric("appendStep<-"))
setReplaceMethod("appendStep", c("SYSargsList"), function(x, after=length(x), ..., value) {
  lengx <- length(x)
  after <- after
  if(stepName(value) %in% stepName(x)) stop("Steps Names need to be unique.")
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
      x$statusWF <- c(value$statusWF, x$statusWF)
      x$dependency <- c(dependency(value), x$dependency)
      x$outfiles <- c(outfiles(value), x$outfiles)
      x$targets_connection <- c(value$targets_connection, x$targets_connection)
      x$runInfo$directory <- c(value$runInfo, x$runInfo$directory)
    } else if (after >= lengx) {
      x$stepsWF <- c(x$stepsWF, value$stepsWF)
      x$targetsWF <- c(x$targetsWF, targetsWF(value))
      x$statusWF <- c(x$statusWF, value$statusWF)
      x$dependency <- c(x$dependency, dependency(value))
      x$outfiles <- c(x$outfiles, outfiles(value))
      x$targets_connection <- c(x$targets_connection, value$targets_connection)
      x$runInfo$directory <- c(x$runInfo$directory, value$runInfo)
    } else {
      after_tc <- names(x$stepsWF)[1L:after]
      before_tc <- names(x$stepsWF)[(after + 1L):lengx]
      x$targets_connection <- c(x$targets_connection[names(x$targets_connection) %in% after_tc], value$targets_connection, x$targets_connection[names(x$targets_connection) %in% before_tc])
      x$stepsWF <- c(x$stepsWF[1L:after], value$stepsWF, x$stepsWF[(after + 1L):lengx])
      x$targetsWF <- c(x$targetsWF[1L:after], targetsWF(value), x$targetsWF[(after + 1L):lengx])
      x$statusWF <- c(x$statusWF[1L:after], value$statusWF, x$statusWF[(after + 1L):lengx])
      x$dependency <- c(x$dependency[1L:after], dependency(value), x$dependency[(after + 1L):lengx])
      x$outfiles <- c(x$outfiles[1L:after], outfiles(value), x$outfiles[(after + 1L):lengx])
      x$runInfo$directory <- c(x$runInfo$directory[1L:after], value$runInfo, x$runInfo$directory[(after + 1L):lengx])
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
      x$targetsWF <- c(list(DataFrame()), x$targetsWF)
      x$statusWF <- c(list(value$status), x$statusWF)
      x$dependency <- c(value$dependency, x$dependency)
      x$outfiles <- c(list(DataFrame()), x$outfiles)
      x$targets_connection <- c(list(NULL), x$targets_connection)
      x$runInfo$directory <- c(list(NULL), x$runInfo$directory)
    } else if (after >= lengx) {
      x$stepsWF <- c(x$stepsWF, value)
      x$targetsWF <- c(x$targetsWF, list(DataFrame()))
      x$statusWF <- c(x$statusWF, list(value$status))
      x$dependency <- c(x$dependency, value$dependency)
      x$outfiles <- c(x$outfiles, list(DataFrame()))
      x$targets_connection <- c(x$targets_connection, list(NULL))
      x$runInfo$directory <- c(x$runInfo$directory, list(NULL))
    } else {
      after_tc <- names(x$stepsWF)[1L:after]
      before_tc <- names(x$stepsWF)[(after + 1L):lengx]
      x$targets_connection <- c(x$targets_connection[names(x$targets_connection) %in% after_tc], list(NULL), x$targets_connection[names(x$targets_connection) %in% before_tc])
      x$stepsWF <- c(x$stepsWF[1L:after], value, x$stepsWF[(after + 1L):lengx])
      x$targetsWF <- c(x$targetsWF[1L:after], list(DataFrame()) , x$targetsWF[(after + 1L):lengx])
      x$statusWF <- c(x$statusWF[1L:after], list(value$status), x$statusWF[(after + 1L):lengx])
      x$dependency <- c(x$dependency[1L:after], value$dependency, x$dependency[(after + 1L):lengx])
      x$outfiles <- c(x$outfiles[1L:after], list(DataFrame()), x$outfiles[(after + 1L):lengx])
      x$runInfo$directory <- c(x$runInfo$directory[1L:after], list(NULL), x$runInfo$directory[(after + 1L):lengx])
    }
    names(x$stepsWF)[after+1L] <- step_name
    names(x$statusWF)[after+1L] <- step_name
    names(x$dependency)[after+1L] <- step_name
    names(x$targetsWF)[after+1L] <- step_name
    names(x$outfiles)[after+1L] <- step_name
    names(x$targets_connection)[after+1L] <- step_name
    names(x$runInfo$directory)[after+1L] <- step_name
    x <- as(x, "SYSargsList")
  } else stop("Argument 'value' needs to be assigned an object of class 'SYSargsList' OR 'LineWise'.")
 # if(any(duplicated(names(stepsWF(x))))) warning("Duplication is found in names(stepsWF(x)). Consider renaming the steps.")
  x <- .check_write_SYSargsList(x)
  x
})

.check_write_SYSargsList <- function(x){
  if(!inherits(x, "SYSargsList")) stop("Argument 'x' needs to be assigned an object of class 'SYSargsList'.")
  sys.file <- projectInfo(x)$sysargslist
  if(is.null(sys.file)){
    if (interactive()) {
      init <- readline(cat(
        cat(crayon::bgMagenta("** Object 'x' was NOT initialized with SPRproject function! **\n")),
        "Do you like to initialize the project now using the default values?", "Type a number: \n",
        "\n 1. Please create the project \n 2. I don't need the project information \n 3. Quit"
      ))
    } else {
      ## For an non-interactive session
      init <- "1"
    }
    if (init == "1") {
      init <- SPRproject()
      x[["projectInfo"]] <- init$projectInfo
      sys.file <- projectInfo(x)$sysargslist
      write_SYSargsList(x, sys.file, silent=TRUE)
      return(x)
    } else if (init == "2") {
     print("For more details, check help(SRRproject)")
      return(x)
    } else if (init == "3") {
      stop("Quiting...")
    }
  } else if(!is.null(sys.file)){
    sys.file <- projectInfo(x)$sysargslist
    write_SYSargsList(x, sys.file, silent=TRUE)
    return(x)
  }
}

.validationStepConn <- function(x, value){
  ## Check outfiles names
  if(any(
    duplicated(unlist(append(sapply(outfiles(x), function(y) names(y)), sapply(outfiles(value), function(y) names(y)))))))
    stop("'outfiles' columns names needs to be unique")
  ## Check value length
  if(length(value) > 1) stop("One step can be appended in each operation.")
  targesCon <- value$targets_connection[[1]]
  if (!is.null(targesCon[[1]])){
    step <- targesCon[[1]][[1]]
    if(any(!step %in% names(stepsWF(x)))) stop(paste0("'targets' argument needs to be assigned as valid targets file OR the names of a previous step, for example: ", "\n",
                                                      paste0(names(stepsWF(x)), collapse = " OR ")))
    if(length(step)==1){
      targets_name <- paste(colnames(targetsWF(x)[step][[1]]), collapse="|")
      new_targets_col <- targesCon[[2]][[1]][-c(which(grepl(targets_name, targesCon[[2]][[1]])))]
      ## addd skip
      if(all(!new_targets_col %in% colnames(x$outfiles[[step]]))) stop(paste0("'targets_column' argument needs to be assigned as valid column names of a previous step, for example: ", "\n",
                                                                           paste0(colnames(x$outfiles[[step]]), collapse = " OR \n")))
      if(is.null(targesCon[[3]][[1]])){
        old_targets <- x$targetsWF[[step]]
      } else {
        old_targets <- x$targetsWF[[step]][-c(which(grepl(paste(targesCon[[3]][[1]], collapse="|"), colnames(x$targetsWF[[step]]))))]
      }
      new_targets <- cbind(x$outfiles[[step]][new_targets_col], old_targets)
      new_targetsheader <- targetsheader(x, step)
      ## DOUBLE CONNECTION
    } else if(length(step) > 1){
      targets_list <- sapply(step, function(y) targetsWF(x)[[y]])
      targets_list_name <- unique(unlist(lapply(targets_list, function(y) names(y))))
      old_targets <- Reduce(function(x, y) merge(x, y, by=targets_list_name, all=TRUE), targets_list)
      targets_name <- paste(targets_list_name, collapse="|")
      new_targets_col <- targesCon[[2]][[1]][-c(which(grepl(targets_name, targesCon[[2]][[1]])))]
      colnames_outfiles <- sapply(outfiles(x), function(y) names(y))
      if(!all(new_targets_col %in% colnames_outfiles)) stop(paste0("'targets_column' argument needs to be assigned as valid column names of a previous step, for example: ", "\n",
                                                                   paste0(colnames_outfiles, collapse = " OR \n")))
      if(is.null(targesCon[[3]][[1]])){
        old_targets <- old_targets
      } else {
        old_targets <- old_targets[-c(which(grepl(paste(targesCon[[3]][[1]], collapse="|"), colnames(old_targets))))]
      }
      new_col_list <- lapply(step, function(y) outfiles(x)[[y]])
      new_targets <- cbind(new_col_list, old_targets)
      new_targetsheader <- sapply(step, function(y) targetsheader(x, y))[1]
      names(new_targetsheader) <- "targetsheader"
    }
    WF <- value$stepsWF[[1]]
    #inputvars_v <- unlist(WF$inputvars)
    ## TODO: check inputvars...
    WF2 <- updateWF(WF, new_targets= targets.as.list(data.frame(new_targets)), new_targetsheader=new_targetsheader, inputvars=WF$inputvars, write.yaml = FALSE)
    value <- sysargslist(value)
    value$stepsWF[[1]] <- WF2
    value$targetsWF[[1]] <- as(WF2, "DataFrame")
    value$outfiles[[1]] <- output.as.df(WF2)
    value$statusWF[[1]] <- WF2$status
    value <- as(value, "SYSargsList")
  } #else if (is.null(targesCon[[1]])){
  #   value <- sysargslist(value)
  #   value$targets_connection <- list(NULL)
  #   names(value$targets_connection) <- names(value$stepsWF)
  #   value <- as(value, "SYSargsList")
  # }
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
  ## runInfo
  if("env" %in% names(value$runInfo)){
    value[["runInfo"]] <- value$runInfo$directory
  }
  
  return(value)
}

## USage:
# appendStep(sal) <- SYSargsList(WF)
# appendStep(sal, after=0) <- SYSargsList(WF)
# sal
# appendStep(sal, after=0, step_index="test_11") <- SYSargsList(WF)

output.as.df <- function(x) {
  out_x <- output(x)
  out_x <- S4Vectors::DataFrame(matrix(unlist(out_x), nrow=length(out_x), byrow=TRUE))
  colnames(out_x) <- x$files$output_names
  return(out_x)
}


setGeneric(name="replaceStep<-", def=function(x, step, step_name="default", value) standardGeneric("replaceStep<-"))
setReplaceMethod("replaceStep", c("SYSargsList"), function(x, step, step_name="default", value) {
  if(!inherits(value, "SYSargsList")) stop("Argument 'value' needs to be assigned an object of class 'SYSargsList'.")
  if(length(value) > 1) stop("Argument 'value' cannot have 'length(value) > 1")
  if(inherits(step, "numeric")){
    if(step > length(x)) stop(paste0("Argument 'step' cannot be greater than ", length(x)))
  } else if(inherits(step, "character")){
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
  x <- .check_write_SYSargsList(x)
  x
})

# replaceStep(sal, 1) <- sal[1]
# sal

setGeneric(name="stepsWF<-", def=function(x, step, ..., value) standardGeneric("stepsWF<-"))
setReplaceMethod("stepsWF", c("SYSargsList"), function(x, step, ..., value) {
    x@stepsWF[[step]] <- value
    x <- .check_write_SYSargsList(x)
    x
})

setGeneric(name="renameStep<-", def=function(x, step, ..., value) standardGeneric("renameStep<-"))
setReplaceMethod("renameStep", c("SYSargsList"), function(x, step, ..., value) {
  if(length(step)!=length(value)) stop("value argument needs to be the same length of the step for rename")
  if(inherits(value, "character")){
    names(x@stepsWF)[step] <- value
    names(x@statusWF)[step] <- value
    names(x@dependency)[step] <- value
    names(x@targetsWF)[step] <- value
    names(x@outfiles)[step] <- value
    #names(x@SEobj)[step] <- value
    names(x@targets_connection)[step] <- value
    if(!is.null(x$runInfo$directory)){
      names(x@runInfo$directory)[step] <- value
    }
  }  else {
    stop("Replace value needs to be assigned an 'character' name for the workflow step.")
  }
  if(!is.null(x$projectInfo$sysargslist)){
    x <- .check_write_SYSargsList(x)
  }
  x
})

setGeneric(name="statusWF<-", def=function(x, step, ..., value) standardGeneric("statusWF<-"))
setReplaceMethod("statusWF", c("SYSargsList"), function(x, step, ..., value) {
    x@statusWF[[step]] <- value
    x <- .check_write_SYSargsList(x)
    x
})
