####################################################################
## Functions to construct SYSargsList objects and other utilities ##
####################################################################

##########################
## SPRproject function ##
##########################
## Detection and creation of the directory of the project
SPRproject <- function(projPath = "./", overwrite = FALSE, silent = FALSE, 
                       tempdir = FALSE) {
  if (!class(projPath) == "character") stop("Provide valid 'projPath' PATH.")
  if (tempdir == TRUE) {
    projPath <- tempdir()
    projPath <- file.path(projPath, ".SPRproject")
  } else if (tempdir == FALSE) {
    projPath <- file.path(projPath, ".SPRproject")
  }
  if (file.exists(projPath) == FALSE) {
    dir.create(projPath, recursive = TRUE)
    if (silent != TRUE) cat("Creating directory '", normalizePath(projPath), "'", sep = "")
  } else if (file.exists(projPath) == TRUE) {
    if (overwrite == FALSE) {
      if (file.exists(projPath)) {
        stop("/.SPRproject already exists. '/.SPRproject' can be overwritten OR restart the project where it was stopped.")
      }
    } else if (overwrite == TRUE) {
      unlink(projPath, recursive = TRUE)
      dir.create(projPath, recursive = TRUE)
      if (silent != TRUE) cat("Directory '", normalizePath(projPath), "' was overwritten", sep = "", "\n")
    }
  }
  return(normalizePath(projPath))
}

## Usage
# projPath = "./"; overwrite = FALSE; silent = FALSE; tempdir = FALSE
# SPRproj <- SPRproject()
# SRPproj <- SPRproject(projPath="./", overwrite=FALSE, tempdir=TRUE)

##########################
## initProject function ##
##########################
## This function detects an existing project or creates a project structure on the path provide
initProject <- function(projPath = "./", data = "data", param = "param", results = "results",
                        script = NULL, targets = NULL, filename = "SPRconfig.yml",
                        overwrite = FALSE, silent = FALSE, tempdir=FALSE) {
  project <- list(
    project = projPath, 
    data = file.path(normalizePath(projPath), data),
    param = file.path(normalizePath(projPath), param), 
    results = file.path(normalizePath(projPath), results)
  )
  path <- sapply(project, function(x) suppressMessages(tryPath(x)))
  create <- NULL
  for (i in seq_along(path)) {
    if (is.null(path[[i]])) create <- c(create, project[i])
  }
  ## Question
  if (!is.null(create)) {
    ## For an interactive() session
    if (interactive()) {
      dir_create <- readline(cat(
        "There is no directory called", "\n", paste(names(create), collapse = " OR ", sep="\n"), "\n",
        "Would you like to create this directory now? Type a number: \n 1. Yes \n 2. No \n"
      ))
    } else {
      ## For an non-interactive session
      dir_create <- "1"
    }
    for (i in seq_along(create)) {
      if (dir_create == "1") {
        dir.create(create[[i]], recursive = TRUE)
        print(paste("Creating directory:", create[[i]]))
      } else if (dir_create == 2) {
        stop("Aborting project creation.")
      }
    }
  }
  path <- sapply(project, function(x) suppressMessages(tryPath(x)))
  project_path <- list(class = "Directory", path = normalizePath(path[1]))
  data_path <- list(class = "Directory", path = normalizePath(path[2]))
  param_path <- list(class = "Directory", path = normalizePath(path[3]))
  results_path <- list(class = "Directory", path = normalizePath(path[4]))
  targets <- list(class = "File", path = targets)
  script <- list(class = "File", path = script)
  initProj <- list(
    project = project_path, data = data_path,
    param = param_path, results = results_path, targets = targets,
    script = script
  )
  SPRproj <- SPRproject(projPath = projPath, overwrite = overwrite, silent = silent, tempdir = tempdir)
  if (all(file.exists(file.path(projPath, filename)) & overwrite == FALSE)) 
    stop(paste("I am not allowed to overwrite files; please delete existing file:", filename, "or set 'overwrite=TRUE'"))
  yaml::write_yaml(x = initProj, file = file.path(projPath, filename))
  if (silent != TRUE) cat("\n", "Project started with success: ./SPRproject and", filename, "were created.")
  initProj <- c(initProj, list(SPRproj = SPRproj, SPRconfig = normalizePath(file.path(projPath, filename))))
  return(initProj)
}

## Usage
# projPath = "./"; data = "data"; param = "param"; results = "results"; script = NULL; targets = NULL;
# subProj = FALSE; filename = "SPRconfig.yml"; overwrite = FALSE; silent = FALSE
# SPRconfig <- initProject(projPath="./", targets="targets.txt", script="systemPipeRNAseq.Rmd", overwrite=TRUE, silent=FALSE)
# SPRconfig <- initProject(projPath="./", targets="targets.txt", script="systemPipeRNAseq.Rmd", overwrite=TRUE, silent=TRUE)
# SPRconfig <- initProject(projPath="./", targets="targets.txt", script="systemPipeRNAseq.Rmd", overwrite=FALSE, silent=FALSE)
# SPRconfig <- initProject(projPath="./", data="./data2", param="./param2", results="./results")

#####################
## initWF function ##
#####################
initWF <- function(sprconfig = NULL, script=NULL, targets = NULL, silent = FALSE,
                   overwrite = FALSE) {
  ## Validations
  if (all(!is.null(sprconfig) && !file.exists(sprconfig) == TRUE)) stop("Provide valid 'sprconfig' file. Check the file PATH.")
  if (all(!is.null(script) && !file.exists(script) == TRUE)) stop("Provide valid 'script' file. Check the file PATH.")
  if (all(!is.null(targets) && !file.exists(targets) == TRUE)) stop("Provide valid 'script' file. Check the file PATH.")
  ## Building an 'SYSargsList' empty container
  init <- as(SYScreate("SYSargsList"), "list")
  ## detects an existing project or creates a project structure on getwd() OR it uses all the info from the sprconfig argument
  if (is.null(sprconfig)) {
    sprconfig <- initProject(projPath = "./", targets = targets, script = script, overwrite = overwrite, silent = silent)
  } else if (!is.null(sprconfig)) {
    tryCatch(.sprconfigCheck(sprconfig),
             warning = function(w) {
               w$message <- paste("Please check the", sprconfig, "file. Some file path is missing. For more details: 'help(initWF)'")
               stop(w)
             }
    )
    sprconfig <- yaml::read_yaml(sprconfig, eval.expr = TRUE)
  }
  ## TODO: subproject structure
  init$projectWF <- list(
    project = sprconfig$project$path, data = sprconfig$data$path, param = sprconfig$param$path,
    results = sprconfig$results$path, logs=NA
  )
  init$sprconfig <- sprconfig
  init <- as(init, "SYSargsList")
  ## TODO: if init$sprconfig$script$path==TRUE, creates LINEwise object...
  return(as(init, "SYSargsList"))
}

## Usage:
# sprconfig = NULL; script=NULL; targets = NULL; silent = FALSE; overwrite = FALSE
# sal <- initWF(script="systemPipeRNAseq.Rmd", targets = "targets.txt", overwrite = TRUE)
# sal <- initWF(sprconfig = "SPRconfig.yml")
# sal
# sal <- initWF(sprconfig = NULL, overwrite = TRUE)
# sal <- initWF()

########################################################
## SYSargsList ##
########################################################
SYSargsList <- function(args=NULL, 
                        targets=NULL, wf_file=NULL, input_file=NULL, dir_path=".", inputvars=NULL,
                        step_index="default",
                        silent = FALSE, restartProject = TRUE, 
                        rm_targets_col = NULL, dependency=NULL) {
  if(all(is.null(args) && is.null(wf_file) && is.null(input_file))){
    sal <- initWF(sprconfig = NULL, script=NULL, targets = NULL, silent = silent,
                  overwrite = restartProject)
    return(sal)
  } else if (!is.null(args)){
    if(any(length(cmdlist(args)[[1]])==0)) stop("Argument 'args' needs to be assigned an object of class 'SYSargs2'.") 
    sal <- SYScreate("SYSargsList"); sal <- sysargslist(sal)
    if(step_index=="default"){
      step_index <- "Step_1"
    } else {
      step_index <- step_index
    }
    sal$stepsWF <- list(args)
    names(sal$stepsWF) <- step_index
  } else if(all(!is.null(wf_file) && !is.null(input_file))){
    ## targets
    if(is.null(targets)){
      targets <- targets
    } else if(is(targets, "character")){
      if(file.exists(targets)){
        targets <- targets
      } else {
        targets_step <- targets
        targets <- NULL
      }
    }
    WF <- loadWorkflow(targets=targets, wf_file=wf_file, 
                       input_file=input_file, dir_path=dir_path)
    WF <- renderWF(WF, inputvars=inputvars)
    sal <- SYScreate("SYSargsList"); sal <- sysargslist(sal)
    if(step_index=="default"){
      step_index <- "Step_1"
    } else {
      step_index <- step_index
    }
    sal$stepsWF <- list(WF)
    names(sal$stepsWF) <- step_index
    if(exists("targets_step")){
      sal$targets_connection <- list(targets_step=targets_step)
      new_targets_col <- names(inputvars)
      if(is.null(new_targets_col)) 
        stop("inputvars argument need to be assigned to the output column names from the previous step specified on the targets argument")
      sal$targets_connection["new_targets_col"] <- list(new_targets_col)
      if(!is.null(rm_targets_col))
        sal$targets_connection["rm_targets_col"] <- list(rm_targets_col)
    }
  }
  ## dependency
  if(is.null(dependency)){
    sal$dependency <- list(NA)
  } else {
    sal$dependency <- list(dependency)
  }
  names(sal$dependency) <- names(sal$stepsWF)
  ## statusWF
  sal$statusWF <- list(list(status="Pending"))
  names(sal$statusWF) <- names(sal$stepsWF)
  ## outfiles
  if(length(sal$stepsWF) >0) {
    sal$outfiles <- outList2DF(sal)
    ## targets
    if(length(targets(sal$stepsWF[[1]])) > 0 ) {
      sal$targetsWF <- list(as(sal$stepsWF[[1]], "data.frame"))
      names(sal$targetsWF) <- names(sal$stepsWF)
    } else {
      sal$targetsWF <- list(NULL)
      names(sal$targetsWF) <- names(sal$stepsWF)
    }
  }
  sal <- as(sal, "SYSargsList")
  return(sal)
}

########################
## internal function ##
########################

outList2DF <- function(sal){
  sal <- as(sal, "SYSargsList")
  out <- sapply(names(stepsWF(sal)), function(x) list(NULL))
  for(i in seq_along(stepsWF(sal))){
    l_out <- output(stepsWF(sal)[[i]])
    out[[i]] <- S4Vectors::DataFrame(matrix(unlist(l_out), nrow=length(l_out), byrow=TRUE))
    if(length(stepsWF(sal)[[i]]$clt) > 1){
      colnames <- sapply(1:length(stepsWF(sal)[[i]]$clt), function(x) 
        names(stepsWF(sal)[[i]]$clt[[x]]$outputs), simplify = TRUE)
    } else{
      colnames <- names(stepsWF(sal)[[i]]$clt[[1]]$outputs)
    }
    colnames(out[[i]]) <- colnames
  }
  return(out)
}

.outputTargets <- function(args, fromStep, index=1, toStep, replace=c("FileName")){
  if(!is(args, "SYSargsList")) stop("Argument 'args' needs to be assigned an object of class 'SYSargsList'")
  outputfiles <- outfiles(args[fromStep])[[1]]
  if(length(targetsWF(args)[[fromStep]])>0){
    df <- targetsWF(args)[[fromStep]]
    df[replace] <- outputfiles
  }
  return(df)
}

########################
## ConfigWF function ##
########################
configWF <- function(x, input_steps = "ALL", exclude_steps = NULL, silent = FALSE, ...) {
  ## Validations
  if (class(x) != "SYSargsList") stop("Argument 'x' needs to be assigned an object of class 'SYSargsList'")
  capture.output(steps_all <- subsetRmd(Rmd = x$sysconfig$script$path), file = ".SYSproject/.NULL") ## TODO: refazer
  if ("ALL" %in% input_steps) {
    input_steps <- paste0(steps_all$t_number[1], ":", steps_all$t_number[length(steps_all$t_number)])
    save_rmd <- FALSE
    Rmd_outfile <- NULL
    Rmd_path <- x$sysconfig$script$path
  } else {
    input_steps <- input_steps
    save_rmd <- TRUE
    Rmd_outfile <- paste0(
      .getPath(x$sysconfig$script$path), "/", paste0(format(Sys.time(), "%b%d%Y")),
      "_", basename(x$sysconfig$script$path)
    )
    Rmd_path <- Rmd_outfile
  }
  capture.output(steps_all <- subsetRmd(
    Rmd = x$sysconfig$script$path, input_steps = input_steps,
    exclude_steps = exclude_steps, save_Rmd = save_rmd, Rmd_outfile = Rmd_outfile
  ), file = ".NULL") ## TODO: refazer
  steps_number <- steps_all$t_number[steps_all$selected]
  ## Save input_steps in the SYSconfig.yml
  if (!any(names(x$sysconfig) %in% c("input_steps", "exclude_steps"))) {
    input_file <- x$sysconfig$SYSconfig
    param <- list(input_steps = input_steps, exclude_steps = exclude_steps)
    input <- config.param(input_file = input_file, param, file = "append", silent = TRUE)
    x <- as(x, "list")
    x$sysconfig <- input
    x <- as(x, "SYSargsList")
  }
  ## STEPS
  if (is.null(x$sysconfig$script)) {
    cat("empty list")
  } else {
    steps <- paste0(
      strrep("    ", (as.numeric(steps_all$t_lvl[steps_all$selected]) - 1)), steps_all$t_number[steps_all$selected],
      ".", steps_all$t_text[steps_all$selected]
    )
  }
  names(steps_number) <- steps
  ## CODESteps
  code <- steps_all$code[steps_all$selected]
  names(code) <- steps
  x <- as(x, "list")
  script_path <- x$sysconfig$script$path
  x$sysconfig$script[["path"]] <- Rmd_path
  x$sysconfig$script[["orginal"]] <- script_path
  x$stepsWF <- steps_number
  x$codeSteps <- code
  x$dataWF <- steps_all
  # x$sysconfig$script$path <- Rmd_outfile
  return(as(x, "SYSargsList"))
}

## Usage:
# sysargslist <- configWF(x=sysargslist)

#####################
## runWF function ##
#####################

#####################
## runWF function ##
#####################
runWF <- function(args, warning.stop=FALSE, error.stop=TRUE, silent=FALSE) {
  # Validations
  if (!inherits(args, "SYSargsList")) stop("Argument 'args' needs to be assigned an object of class 'SYSargsList'")
  sysproj <- paste(projectWF(args)$project, ".SPRproject", sep = "/")
  if (!file.exists(sysproj) == TRUE) stop("Project was not initialized with the 'initWF' function.")
  ## check dependency
  for(i in seq_along(dependency(args))){
    if(all(!is.na(dependency(args)[[i]]))){
      dep_names <- dependency(args)[[i]]
      if(any(!dep_names %in% names(stepsWF(args)))) stop(
        "'args' has dependency on the following steps:", "\n",
        paste0(dep_names, collapse = " AND "))
    }
  }
  ## Check status, rerun if warning or error...
  file_log <- file.path(sysproj, paste0("_logWF_", format(Sys.time(), "%b%d%Y_%H%M")))
  for (i in seq_along(stepsWF(args))){
    step <- names(stepsWF(args)[i])
    print(paste0("Running Step: ", step))
    df.out <- outfiles(args)[[i]]
    cat(c(paste0("# ", step), "\n",
          paste0("Time: ", paste0(format(Sys.time(), "%b%d%Y_%H%Ms%S")))), file = file_log, sep = "\n", append = TRUE)
    run_status <- as.character()
    status_files_check <- sapply(names(cmdlist(args)[[i]]), function(x) list(NULL))
    for (j in seq_along(cmdlist(args)[[i]])) {
      cat(c(paste0("## ", names(cmdlist(args)[[i]][j])),
            paste0("Time: ", paste0(format(Sys.time(), "%b%d%Y_%H%Ms%S"))), "\n",
            "### Code: ",
            "```{r, eval=FALSE} ",
            cmdlist(args)[[i]][[j]][[1]],
            "```",
            "### Stdout: ",
            "```{r, eval=FALSE}" ), file = file_log, sep = "\n", append = TRUE)

      ## Create an object for executable
      command <- gsub(" .*", "", as.character(cmdlist(args)[[i]][[j]][[1]]))
      commandargs <- gsub("^.*? ", "",as.character(cmdlist(args)[[i]][[j]][[1]]))
      ## Check if the command is in the PATH
      if(!command == c("bash")){
        tryCatch(system(command, ignore.stdout = TRUE, ignore.stderr = TRUE), warning=function(w) cat(paste0("ERROR: ", "\n", command, ": command not found. ", '\n', "Please make sure to configure your PATH environment variable according to the software in use."), "\n"))
      }
      #             ## Run executable
      # if(command %in% "bwa") {
      #     stdout <- system2(command, args=commandargs, stdout=TRUE, stderr=FALSE)
      # } else if(command %in% c("bash")) {
      #     stdout <- system(paste(command, commandargs))
      # } else if(isTRUE(grep('\\$', command)==1)) {
      #     stdout <- system(paste(command, commandargs))
      # } else {
      #     stdout <- system2(command, args=commandargs, stdout=TRUE, stderr=TRUE)
      # }
      #
      capture.output(capture.output(.tryCmd(x=system2(command, args=commandargs, stdout=TRUE, stderr=TRUE), file_log), file=file_log, append = TRUE),
                     file=file.path(sysproj, "err"), type = "message", append = FALSE)
      err <- readLines(file.path(sysproj, "err"))
      if(err=="WARNING"){
        if(warning.stop==TRUE) {
          message("Caught an warning!")
          break() }
      } else if(err=="ERROR"){
        if(error.stop==TRUE) {
          message("Caught an error!")
          break()
        }
      }
      ## Printing sample name
      print(paste0("Sample ", names(cmdlist(args)[[i]][j]), ": ", err))
      cat("```", file = file_log, sep = "\n", append = TRUE)
      cat("\n", file = file_log, append = TRUE)
      run_status <- c(run_status, err)
      unlink(file.path(sysproj, "err"))
      status_files <- outfiles(args)[[i]][j,]
      if(length(status_files)==1){
        status_files_check[[j]] <- file.exists(as.character(as.data.frame(status_files[1])))
      } else{
        status_files_check[[j]] <- file.exists(as.character(as.data.frame(status_files[1,])))
      }
    }
    args <- as(args, "list")
    df.status <- data.frame(matrix(unlist(status_files_check), nrow=length(status_files_check), byrow=T))
    df.status1 <- rowSums(df.status)
    df.status2 <- cbind(Running_Status=run_status, Existing_Files = df.status1, Missing_Files=length(df.status)-df.status1)
    rownames(df.status2) <- names(args$stepsWF[[i]]$cmdlist)
    args$statusWF[[i]] <- list(status=c(unique(run_status)))
    args$statusWF[[i]]$Complete <- df.status2
    args <- as(args, "SYSargsList")
  }

  if(silent != TRUE) cat("\t", "Written content of 'stdout and stderr' to file:", file_log, "\n")
  args <- as(args, "list")
  args$projectWF[["logs"]] <- file_log
  return(as(args, "SYSargsList"))
}


#############################
## .tryCatch function ##
#############################
.tryCmd <- function(x, file=NULL) {
  if(is.null(file)) file=tempfile()
  tryCatch(
    expr = {
      # cat("## Output", append = TRUE, file = file, "\n")
      capture.output(out <- eval(x, envir = globalenv()), file = file, append = TRUE)
      message("DONE")
      return(out)
    },
    error = function(e) {
      print(e)
      message("ERROR")
      return("Caught an error!")
    },
    warning = function(w) {
      print(w)
      message("WARNING")
      return("Caught an warning!")
    }
  )
}


# runWF <- function(sysargslist, steps = "ALL", warning.stop=FALSE, error.stop=TRUE, silent=FALSE) {
#   ## Validations
#   if (class(sysargslist) != "SYSargsList") stop("Argument 'sysargslist' needs to be assigned an object of class 'SYSargsList'")
#   sysproj <- paste(sysargslist$projectWF$project, ".SYSproject", sep = "/")
#   if (!file.exists(sysproj) == TRUE) stop("Project was not initialized with the 'initWF' function.")
#   ## Storing steps list
#   codeList <- sysargslist$codeSteps
#   if ("ALL" %in% steps) {
#     stepslist <- seq_along(sysargslist$stepsWF)
#     level <- sysargslist$stepsWF[`stepslist`]
#   } else {
#     stepslist <- steps
#     t_lvl <- sysargslist$dataWF$t_lvl
#     names(t_lvl) <- sysargslist$dataWF$t_number
#     stepslist <- .parse_step(t_lvl, input_steps = steps)
#     level <- sysargslist$stepsWF[`stepslist`]
#   }
#   file_log <- file.path(sysproj, paste0("_logWF_", format(Sys.time(), "%b%d%Y_%H%M")))
#   status <- as.character()
#   for (i in stepslist) {
#     cat(c(paste0("# ", sub(".*[.]\\s", "", names(codeList[i]))),
#           paste0("Time: ", paste0(format(Sys.time(), "%b%d%Y_%H%Ms%S"))), "\n",
#           "## Code: ",
#           "```{r, eval=FALSE} ",
#           codeList[[i]],
#           "```",
#           "## Output: ",
#           "```{r, eval=FALSE}" ), file = file_log, sep = "\n", append = TRUE)
#     capture.output(capture.output(.tryCatch(codeList[[i]], file_log), file=file_log, append = TRUE),
#                    file=file.path(sysproj, "err"), type = "message", append = FALSE)
#     err <- readLines(file.path(sysproj, "err"))
#     if(err=="WARNING"){
#       if(warning.stop==TRUE) {
#         message("Caught an warning!")
#         break() }
#     } else if(err=="ERROR"){
#       if(error.stop==TRUE) {
#         message("Caught an error!")
#         break()
#       }
#     }
#     err <- paste0("Step ", names(level)[i], ": ",  err)
#     print(err)
#     cat("```", file = file_log, sep = "\n", append = TRUE)
#     cat("\n", file = file_log, append = TRUE)
#     status <- c(status, err)
#     unlink(file.path(sysproj, "err"))
#   }
#   if(silent != TRUE) cat("\t", "Written content of 'stdout and stderr' to file:", file_log, "\n")
#   sysargslist <- as(sysargslist, "list")
#   sysargslist$statusWF <- list(statusWF = status)
#   sysargslist$projectWF[["SYSproject"]] <- sysproj
#   sysargslist$projectWF[["logs"]] <- file_log
#   return(as(sysargslist, "SYSargsList"))
# }

## Usage:
# script <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
# sysargslist <- initWF(script="systemPipeRNAseq.Rmd", overwrite = TRUE)
# sysargslist <- configWF(x=sysargslist, input_steps = "1:2")
# sysargslist <- runWF(sysargslist = sysargslist, steps = "ALL")
# sysargslist <- runWF(sysargslist = sysargslist, steps = "1:2")
# sysargslist <- initWF(script="systemPipeRNAseq.Rmd", overwrite = T) %>%
#   configWF( input_steps = "1:5") %>%
#   runWF(steps = 1:2)
# sysargslist <- renderReport(sysargslist)
# sysargslist <- renderLogs(sysargslist)

###########################
## renderReport function ##
###########################
## type: c("pdf_document", "html_document")
renderReport <- function(sysargslist, type = c("html_document"), silent=FALSE) {
  file <- sysargslist$sysconfig$script$path
  if (!file.exists(file) == TRUE) stop("Provide valid 'sysargslist' object. Check the initialization of the project.")
  evalCode(infile = file, eval = FALSE, output = file)
 # rmarkdown::render(input = file, c(paste0("BiocStyle::", type)), quiet = TRUE, envir = new.env())
  rmarkdown::render(input = file, c(paste0(type)), quiet = TRUE, envir = new.env())
  file_path <- .getPath(file)
  file_out <- .getFileName(file)
  ext <- strsplit(basename(type), split = "\\_")[[1]][1]
  sysargslist <- as(sysargslist, "list")
  sysargslist$projectWF[["Report"]] <- normalizePath(paste0(file_path, "/", file_out, ".", ext))
  if(silent != TRUE) cat("\t", "Written content of 'Report' to file:", paste0(file_out, ".", ext), "\n")
  return(as(sysargslist, "SYSargsList"))
}

## Usage
# renderReport(sysargslist)

###########################
## renderLogs function ##
###########################
## type: c("pdf_document", "html_document")
renderLogs <- function(sysargslist, type = "html_document", fileName="default", silent=FALSE) {
  file <-  projectWF(sysargslist)$logs
  if (!file.exists(file) == TRUE) stop("Provide valid 'sysargslist' object. Check the initialization of the project.")
  if(fileName=="default"){
    fileName <- file.path(projectWF(sysargslist)$project, paste0("logs_", format(Sys.time(), "%b%d%Y_%H%M"),".Rmd"))
  } else {
    fileName <- fileName
  }
  log <- readLines(file)
  writeLines(c(
    "---",
    "title: 'Report'",
    paste0("date: 'Last update: ", format(Sys.time(), '%d %B, %Y'), "'"),
    "output:",
    paste0("  BiocStyle::", type, ":"),
    "    toc_float: true",
    "    code_folding: show",
    "package: systemPipeR",
    "fontsize: 14pt",
    "---", 
    log), 
    con = fileName)
  #rmarkdown::render(input = fileName, c(paste0("BiocStyle::", type)), quiet = TRUE, envir = new.env())
  rmarkdown::render(input = fileName, c(paste0(type)), quiet = TRUE, envir = new.env())
  file_path <- .getPath(fileName)
  file_out <- .getFileName(fileName)
  if(type=="html_document"){
    ext <- "html"
  } else {
    ext <- "pdf"
  }
  sysargslist <- as(sysargslist, "list")
  sysargslist$projectWF[["Report_Logs"]] <- file.path(file_path, paste(file_out, ext, sep="."))
  if(silent != TRUE) cat("\t", "Written content of 'Report' to file:", paste(file_out, ext, sep="."), "\n")
  return(as(sysargslist, "SYSargsList"))
}

## Usage
# sysargslist <- renderLogs(sysargslist)

########################
## subsetRmd function ##
########################
## This function allows the user to subset the Rmarkdown and select the corresponding text and R chunk code.
## The steps definition is based on the block-level elements, defined by the `#`.
## For example: # First-level header represents step number 1;
## For example: ## Second-level header represents step number 1.1;
## For example: ### Third-level header represents step number 1.1.1
## Arguments:
# Rmd: a character vector name and path of the Rmd file
# Rmd_outfile: a character vector name and path of the output Rmd file
# input_steps: a character vector of all steps desires to preserve on the output file.
# Default is NULL. If no input_steps defined, it will return only list steps on Rmd and datagram, however, all titles are unselected (FALSE) and no output file will be saved.
# It can be used the symbol ":" to select many steps in sequence, for example, input_steps=1:5.2, from step 1 to step 5.2.
# The symbol "." represents the substeps and symbol "," is used to separate selections.
# Jump from a major step to sub-step is supported, but if a major step is selected/excluded, all sub-steps of this major step will be selected/excluded. Repeatedly selected steps will only result in a unique step.  It is recommended to put major steps in `input_steps`, like '1:4, 6:8, 10'; and unwanted sub-steps in `exclude_steps`, like '1.1, 3.1.1-3.1.3, 6.5'.
# Reverse selecting is supported e.g. '10:1'.
# exclude_steps: a character vector of all steps desires to remove on the output file.
# save_Rmd: logical, default TRUE. If set FALSE, list new selected tiles and exit.
## Value
## It returns a data frame of title levels, title numbers, title text, whether it is selected, and R code under this title.

subsetRmd <- function(Rmd, input_steps = NULL, exclude_steps = NULL, Rmd_outfile = NULL, save_Rmd = TRUE) {
  # function start, check inputs
  assertthat::assert_that(file.exists(Rmd))
  if (assertthat::not_empty(input_steps)) assertthat::assert_that(assertthat::is.string(input_steps))
  if (assertthat::not_empty(Rmd_outfile)) assertthat::assert_that(file.exists(dirname(Rmd_outfile)))
  if (assertthat::not_empty(exclude_steps)) assertthat::assert_that(assertthat::is.string(exclude_steps))
  # default out behavior, in ISO 8601 time format
  if (is.null(Rmd_outfile)) Rmd_outfile <- paste0("new", format(Sys.time(), "%Y%m%d_%H%M%S"), basename(Rmd))
  # read file
  file <- readLines(Rmd)
  # check for proper start and end
  t_start <- file %>% stringr::str_which("^#")
  if (length(t_start) == 0) stop("This Rmd does not have any '#' titles")
  # get code chunks
  chunk_start <- file %>% stringr::str_which("^```\\{.*\\}.*")
  chunk_end <- file %>% stringr::str_which("^```[[:blank:]]{0,}$")
  if (length(chunk_start) != length(chunk_end)) stop("unmatched number of code chunk starts and ends")
  for (i in seq_along(chunk_start)[-length(chunk_end)]) {
    if (chunk_start[i + 1] <= chunk_end[i]) stop(paste("A code chunk does not end: chunk line", chunk_start[i + 1]))
  }
  # remove '#' titles in code chunk
  t_start <- t_start[!unlist(lapply(t_start, function(x) any(x >= chunk_start & x <= chunk_end)))]
  # get end
  t_end <- append((t_start - 1)[c(-1)], length(file))
  # get # levels and text
  t_text <- file[t_start] %>% stringr::str_remove("^[#]+")
  t_lvl <- file[t_start] %>%
    stringr::str_extract("^[#]+") %>%
    nchar()
  # parse levels
  for (lvl in unique(t_lvl)[-length(t_lvl)]) {
    if (lvl == min(unique(t_lvl))) {
      step_main <- which(t_lvl == lvl)
      names(t_lvl)[step_main] <- names(step_main) <- seq_along(step_main)
      step_main <- append(step_main, 9999)
    }
    sub_lvl <- lvl
    while (sub_lvl <= max(t_lvl)) {
      step_sub <- which(t_lvl == sub_lvl + 1)
      if (length(step_sub) < 1) {
        sub_lvl <- sub_lvl + 1
      } else {
        break()
      }
    }
    jump_step_glue <- if (sub_lvl - lvl == 0) {
      "."
    } else {
      rep(".1.", sub_lvl - lvl) %>%
        paste0(collapse = "") %>%
        stringr::str_replace_all("\\.\\.", "\\.")
    }
    for (i in seq_along(step_main[-1])) {
      subs <- step_sub[step_sub > step_main[i] & step_sub < step_main[i + 1]]
      names(t_lvl)[subs] <- names(step_sub)[step_sub %in% subs] <- paste0(names(step_main[i]), jump_step_glue, seq_along(subs))
    }
    step_main <- append(step_sub, 9999)
  }
  # get code in lists
  code_list <- lapply(seq_along(t_start), function(t_index) {
    code_start <- chunk_start[chunk_start %in% (t_start[t_index]:t_end[t_index])]
    code_end <- chunk_end[chunk_end %in% (t_start[t_index]:t_end[t_index])]
    code_lines <- lapply(seq_along(code_start), function(code_index) {
      (code_start[code_index] + 1):(code_end[code_index] - 1)
    }) %>% unlist()
    file[code_lines]
  })
  # create a df to store everything
  rmd_df <- data.frame(
    t_lvl = t_lvl, t_number = names(t_lvl),
    t_text = t_text, selected = FALSE,
    row.names = NULL, stringsAsFactors = FALSE
  )
  rmd_df$code <- code_list
  # add sample run/success, step link cols
  rmd_df$no_run <- NA
  rmd_df$no_run <- ifelse(sapply(rmd_df$code, length) == 0, 0, NA)
  rmd_df$no_success <- NA
  rmd_df$no_success <- ifelse(sapply(rmd_df$code, length) == 0, 0, NA)
  rmd_df$link_to <- NA
  rmd_df$link_to[1:(nrow(rmd_df) - 1)] <- rmd_df$t_number[2:nrow(rmd_df)]
  # list all steps if no input_steps
  if (!assertthat::not_empty(input_steps)) {
    cat("No input_steps is given, list all sections and exit\n")
    cat("This file contains following sections\n")
    stringr::str_replace(t_text, "^", paste0(strrep("    ", (t_lvl - 1)), names(t_lvl), " ")) %>%
      paste0(., collapse = "\n") %>%
      stringr::str_replace("$", "\n") %>%
      cat()
    return(rmd_df)
  }
  # parse steps
  index_select <- .parse_step(t_lvl, input_steps)
  index_exclude <- .parse_step(t_lvl, exclude_steps)
  index_final <- index_select[!index_select %in% index_exclude]
  rmd_df$selected[index_final] <- TRUE
  # print again what will be write in the new file
  cat("The following sections are selected\n")
  stringr::str_replace(
    t_text[index_final], "^",
    paste0(
      strrep("    ", (t_lvl[index_final] - 1)),
      names(t_lvl[index_final]), " "
    )
  ) %>%
    paste0(collapse = "\n") %>%
    stringr::str_replace("$", "\n") %>%
    cat()
  # to print new titles and return
  if (save_Rmd == FALSE) {
    return(rmd_df)
  }
  # sebset lines
  t_start[index_final]
  t_end[index_final]
  final_lines <- mapply(seq, t_start[index_final], t_end[index_final]) %>%
    unlist() %>%
    append(1:(t_start[1] - 1), .) %>%
    unique()
  writeLines(file[final_lines], Rmd_outfile)
  cat(paste("File write to", normalizePath(Rmd_outfile), "\n"))
  return(rmd_df)
}

# # Usage:
# Rmd <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
# newRmd <- subsetRmd(Rmd=Rmd, input_steps="1:2.1, 3.2:4, 4:6", exclude_steps="3.1", Rmd_outfile="test_out.Rmd", save_Rmd=TRUE)

######################
## plotWF function ##
######################
plotWF <- function(sysargslist, plot_style = "detect", out_type = "html", out_path = "default", height = NULL, width = NULL) {
  ## Validation for 'sysargslist'
  if (class(sysargslist) != "SYSargsList") stop("Argument 'sysargslist' needs to be assigned an object of class 'SYSargsList'")
  df_wf <- sysargslist$dataWF
  # pre checks
  assert_that(out_type %in% c("html", "png", "svg", "shiny"), msg = "out_type needs to be one of 'html', 'png', 'svg'")
  assert_that(plot_style %in% c("detect", "none", "linear"), msg = "out_type needs to be one of 'detect', 'none', 'linear'")
  assert_that(is.data.frame(df_wf))
  all(c("t_lvl", "t_number", "t_text", "selected", "no_run", "no_success", "link_to") %in% names(df_wf)) %>%
    assert_that(msg = 'One of following columns is missing: "t_lvl" "t_number" "t_text" "selected" "no_run" "no_success" "link_to"')
  if (out_path == "default" & !out_type %in% c("html", "shiny")) {
    assert_that(is.writeable(out_path))
    assert_that(is.count(height) | is.null(height))
    assert_that(is.count(width) | is.null(width))
    out_path <- switch(out_type,
      "svg" = paste0("wfplot", format(Sys.time(), "%Y%m%d%H%M%S"), ".svg"),
      "png" = paste0("wfplot", format(Sys.time(), "%Y%m%d%H%M%S"), ".png")
    )
  }
  df_wf <- df_wf[df_wf$selected == TRUE, ]
  if (nrow(df_wf) == 0) {
    return(cat("no step is selected"))
  }

  wf <- .make_plot(df_wf, plot_style, is_main_branch = FALSE)
  wf <- append(wf, .make_plot(df_wf, plot_style, is_main_branch = TRUE), after = length(wf) - 1)
  # special case for detection style plotting, need to move unneeded nodes out of main branch
  if (plot_style == "detect") wf <- .change_branch(df_wf, wf)
  # collapse entire graph
  # return(wf)
  wf <- paste0(wf, collapse = "")
  # plot
  plot <- switch(out_type,
    "shiny" = dot(wf, return = "verbatim"),
    "html" = dot(wf),
    "svg" = dot(wf, return = "verbatim") %>% rsvg_svg(file = out_path, height = height, width = width),
    "png" = dot(wf, return = "verbatim") %>% charToRaw() %>% rsvg_png(file = out_path, height = height, width = width)
  )
  return(plot)
}

## Usage:
# df_wf <- dataWF(sysargslist)
# df_wf$no_success[3:8] <- 1
# df_wf$no_run[3:5] <- 10
# df_wf$no_run[6:8] <- 1
# df_wf$selected[1:35] <- TRUE
# df_wf$link_to <- NA
# df_wf$link_to[1:(nrow(df_wf) - 1)] <- df_wf$t_number[2:nrow(df_wf)]
# df_wf$link_to[3] <- NA
# df_wf$link_to[1] <- "1.1, 2"
# df_wf$link_to[4] <- "2.1, 3"
# df_wf$link_to[14] <- NA
# df_wf <- df_wf[1:17,]
# df_wf$link_to[8] <- "3, 2.5"
# plotWF(df_wf, plot_style = "linear")

###########################
## config.param function ##
###########################
config.param <- function(input_file = NULL, param, file = "default", silent = FALSE) {
  ## In the case of 'input_file' == character (file)
  if (class(input_file) == "character") {
    if (!file.exists(input_file)) {
      stop("Provide valid 'input_file' file. Check the file PATH.")
    }
    input <- yaml::read_yaml(file.path(input_file), eval.expr = TRUE)
    input <- out_obj <- .replace(input = input, param = param)
    path_file <- normalizePath(input_file)
    out_msg <- c("input_file")
  } else if (class(input_file) == "list") {
    if (is.null(names(param))) {
      stop("for each element of the 'param' list need to assign a name.")
    }
    input <- out_obj <- .replace(input = input_file, param = param)
    path_file <- normalizePath(file) ## TODO find a better solution!
    out_msg <- c("input_file")
  } else if (class(input_file) == "SYSargs2") {
    input <- .replace(input = yamlinput(input_file), param = param)
    dir_path <- .getPath(cwlfiles(input_file)[['yml']])
    if(is.na(cwlfiles(input_file)[['targets']])) {
      targets <- NULL
    } else {
      targets <- targets(input_file)
    }
    args1 <- loadWorkflow(
      targets = targets, wf_file = basename(cwlfiles(input_file)[['cwl']]), 
      input_file = basename(cwlfiles(input_file)[['yml']]), dir_path = dir_path
    )
    args1 <- as(args1, "list")
    args1$yamlinput <- input
    if("ModulesToLoad" %in% names(param))
      for (i in seq_along(param$ModulesToLoad)){
        args1$modules[names(param$ModulesToLoad[i])] <- param$ModulesToLoad[[i]]
      }
    args1 <- out_obj <- as(args1, "SYSargs2")
    out_msg <- c("yamlinput(args1)")
    path_file <- cwlfiles(input_file)[['yml']]
  } else if (class(input_file) == "SYSargsList") { ##TODO
    input <- out_obj <- .replace(input = input_file$sysconfig, param = param)
    path_file <- input_file$projectWF$project
    out_msg <- c("input_file")
  } else if (all(is.null(input_file))) {
    stop("'input_file' need to be defenid as '.yml' OR 'SYSargs2 class' OR 'SYSargsList class' OR 'list class' ")
  }
  ## File and Path to be written
  if (file == "default") {
    run <- dir(.getPath(path_file))
    fileName <- paste0("_run", length(grep(
      paste0(gsub("\\..*", "", basename(path_file)), "_run"),
      run
    )) + 1)
    # path <- paste0(.getPath(path_file), '/', paste(format(Sys.time(), '%b%d%Y_%H%M%S')), '_',
    # basename(path_file), collapse = '/')
    path <- paste0(.getPath(path_file), "/", gsub("\\..*", "", basename(path_file)), fileName, ".",
      gsub(".*\\.", "", basename(path_file)),
      collapse = "/"
    )
  } else if (file == "append") {
    path <- path_file
  } else {
    path <- file
  }
  ## Rename original file
  # file.rename(from = path_file, to = path)
  # if (silent != TRUE) cat("\t", "The original file was renamed to:", "\n", paste0(path), "\n")
  ## Write YML file
  yaml::write_yaml(x = input, file = path)
  if (silent != TRUE) cat("\t", "All the new param + ", out_msg, "were written to:", "\n", paste0(path), "\n")
  if (class(input_file) == "SYSargs2") {
    args1 <- as(args1, "list")
    args1$cwlfiles$yml <- path
    args1 <- as(args1, "SYSargs2")
    out_obj <- args1
  }
  return(out_obj)
}

## Usage:
# Config param files
# targets <- system.file("extdata", "targets.txt", package="systemPipeR")
# input_file <- system.file("extdata", "cwl/hisat2/hisat2-se/hisat2-mapping-se.yml", package="systemPipeR")
# param <- list(thread=10, fq=list(class="File", path="./results2"))
# input <- config.param(input_file=input_file, param, file="default") # Example with 'input_file'

###############
## Utilities ##
###############

########################################################
## SYScreate ##
########################################################
SYScreate <- function(class) {
  if (class == "SYSargs2") {
    SYS.empty <- list(
      targets = data.frame(),
      targetsheader = list(),
      modules = list(),
      wf = list(),
      clt = list(),
      yamlinput = list(),
      cmdlist = list(NULL),
      input = list(),
      output = list(),
      cwlfiles = list(),
      inputvars = list(), 
      cmdToCwl = list()
    )
    return(as(SYS.empty, "SYSargs2"))
  } else if (class == "SYSargsList") {
    SYS.empty <- list(
      sprconfig = list(),
      stepsWF = list(),
      statusWF = list(),
      dependency = list(),
      projectWF = list(),
      targetsWF = list(), 
      SEobj = list(), 
      outfiles = list(), 
      targets_connection=list()
    )
    return(as(SYS.empty, "SYSargsList"))
  } else if (!class %in% c("SYSargs2", "SYSargsList")) 
    stop("Argument 'args' needs to be assigned an character 'SYSargs2' OR 'SYSargsList'")
}
## Usage:
# list <- SYScreate("SYSargsList")
# list <- SYScreate("SYSargs2")

#########################################################################################
## Function to check if the command line / Software is installed and set in your PATH ##
#########################################################################################
tryCL <- function(command) {
  if(command=="fastqc") command <- "fastqc --version"
  tryCatch(
    {
      system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
      print("All set up, proceed!")
    },
    warning = function(w) {
      cat(paste0(
        "ERROR: ", "\n", command, ": COMMAND NOT FOUND. ", "\n",
        "Please make sure to configure your PATH environment variable according to the software in use."
      ), "\n")
    }
  )
}

## Usage:
# tryCL(command="R")
# tryCL(command="blastp")
# tryCL(command="fastqc")

########################################################
## Function to check if the Path (file or dir) exists ##
########################################################
tryPath <- function(path) {
  tryCatch(normalizePath(path),
           warning = function(w) message(paste0(path, ": ", "No such file or directory")),
           error = function(e) message(paste0(path, ": ", "Please provide a valid Path"))
  )
}
## Usage:
# tryPath(path="./")

###########################################################################################
## Function to evaluate (eval=TRUE) or not evaluate (eval=FALSE) R code in the Rmarkdown ##
###########################################################################################
## infile: Name and path of the infile file, format Rmarkdown.
evalCode <- function(infile, eval = TRUE, output) {
  if (!file.exists(infile)) stop("Provide valid 'infile' file. Check the file PATH.")
  if (!"Rmd" %in% .getExt(infile)) stop("Provide Rmarkdown file. Check the file PATH.")
  file <- readLines(infile)
  Rchunk <- grep("^```[{}]", file)
  if (length(Rchunk) == 0) stop("The file does not have any R Chuck")
  for (i in Rchunk) {
    if (grepl("eval", file[i])) {
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
  if(!file.exists(x)) warning("No such file or directory. Check the file PATH.")
  x <- normalizePath(x) 
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
  ext <- strsplit(basename(x), split = "\\.")[[1]]
  ext <- ext[[length(ext)]]
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
  if (!file.exists(x)) warning("No such file or directory. Check the file PATH.")
  filename <- strsplit(basename(x), split = "\\.")[[1]]
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
.sprconfigCheck <- function(sysconfig) {
  if (!file.exists(sysconfig) == TRUE) stop("Provide valid 'sysconfig' file. Check the file PATH.")
  sysconfig <- yaml::read_yaml(sysconfig, eval.expr = TRUE)
  project <- list(
    project = sysconfig$project$path, data = sysconfig$data$path, param = sysconfig$param$path,
    results = sysconfig$results$path
  )
  for (i in seq_along(project)) {
    if (is.null(project[[i]])) {
      warning(paste(names(project[i]), "directory is missing..."))
    } else {
      tryPath(path = project[[i]])
    }
  }
  if (!is.null(sysconfig$targets$path)) tryPath(path = sysconfig$targets$path)
  if (is.null(sysconfig$script)) warning("Workflow is missing...")
  if (!is.null(sysconfig$script$path)) tryPath(path = sysconfig$script$path)
}

## Uage:
# .sysconfigCheck(sysconfig="SYSconfig.yml")

##########################
## .parse_step function ##
##########################
## Internal parse function used in the subsetRmd function
.parse_step <- function(t_lvl, input_steps) {
  t_lvl_name <- names(t_lvl)
  input_steps <- unlist(input_steps %>% stringr::str_remove_all(" ") %>% stringr::str_split(",") %>% list())
  # single steps
  nocolon_steps <- input_steps[stringr::str_which(input_steps, "^[^:]+$")]
  lapply(nocolon_steps, function(x) if (!any(t_lvl_name %in% x)) stop(paste("Step", x, "is not found")))
  # dash linked steps
  dash_list <- NULL
  for (i in stringr::str_which(input_steps, ":")) {
    dash_step <- unlist(stringr::str_split(input_steps[i], ":"))
    dash_parse <- unlist(lapply(dash_step, function(x) {
      which(t_lvl_name %in% x) %>% ifelse(length(.) > 0, ., stop(paste("Step", x, "is not found")))
    })) %>% {
      t_lvl_name[.[1]:.[2]]
    }
    dash_list <- append(dash_list, dash_parse)
  }
  # merge
  all_step_name <- unique(append(nocolon_steps, dash_list))
  # if upper level step is selected, all sub-level steps will be added
  unlist(lapply(all_step_name, function(x) stringr::str_which(t_lvl_name, paste0("^", x, "\\..*")))) %>%
    append(which(t_lvl_name %in% all_step_name)) %>%
    unique() %>%
    sort() %>%
    return()
}

################################
## .find_long_branch function ##
################################
## Internal parse function used in the plotWF function
.find_long_branch <- function(t_number, link_to) {
  track_back <- function(t_number, link_to, track_list) {
    for (each_track_n in seq_along(track_list)) {
      each_track <- track_list[[each_track_n]] %>% unlist()
      previous_t_number <- names(link_to_list[which(sapply(link_to_list, function(x) any(x == t_number[each_track[1]])))])
      for (each_num in seq_along(previous_t_number)) {
        previous_link <- which(t_number == previous_t_number[each_num])
        newtrack <- append(previous_link, each_track)
        if (each_num < 2) {
          track_list[[each_track_n]] <- newtrack
        } else {
          track_list[[each_track_n + each_num - 1]] <- newtrack
        }
      }
    }
    if (length(previous_t_number) == 0) {
      return(track_list)
    }
    track_list <- track_back(t_number, link_to, track_list)
    return(track_list)
  }
  link_to_list <- str_split(link_to, ",") %>% lapply(function(x) str_remove_all(x, " "))
  names(link_to_list) <- t_number
  last_step <- list(length(t_number))
  track_list <- track_back(t_number, link_to, last_step)
  long <- sapply(track_list, function(x) all(c(1, length(t_number)) %in% x)) %>%
    track_list[.] %>%
    sapply(length) %>%
    which.max() %>%
    track_list[.] %>%
    unlist()
  return(long)
}

#########################
## .make_plot function ##
#########################
## Internal parse function used in the plotWF function
.make_plot <- function(df_wf, plot_style, is_main_branch = TRUE) {
  if (is_main_branch) {
    # graph start
    wf <- "subgraph { rank=same;\n"
    df_wf <- switch(plot_style,
      "detect" = df_wf[.find_long_branch(df_wf$t_number, df_wf$link_to), ],
      "none" = df_wf[0, ],
      "linear" = {
        df_wf$link_to[1:(nrow(df_wf) - 1)] <- df_wf$t_number[2:nrow(df_wf)]
        df_wf
      }
    )
  } else {
    wf <- switch(plot_style,
      "detect" = "digraph { rankdir=LR;\n",
      "none"   = "digraph { rankdir=TB;\n",
      "linear" = "digraph { rankdir=LR;\n"
    )
    df_wf <- switch(plot_style,
      "detect" = df_wf[-.find_long_branch(df_wf$t_number, df_wf$link_to), ],
      "none"   = df_wf,
      "linear" = df_wf[0, ]
    )
  }
  if (nrow(df_wf) == 0) {
    return(c(wf, "}"))
  }

  steps <- df_wf$t_number
  step_text <- str_replace_all(df_wf$t_text, '[\'\"]', "\\\\'")
  link_to <- df_wf$link_to
  # reslove 1 to n links
  link_to <- str_split(link_to, ",") %>% lapply(function(x) str_remove_all(x, " "))
  # set up colors
  step_color <- ifelse(df_wf$no_run == 0, "gray",
    ifelse(is.na(df_wf$no_run) | is.na(df_wf$no_success), "black",
      ifelse(df_wf$no_run != df_wf$no_success, "red", "green")
    )
  )
  # dot language
  # add steps
  for (t in seq_along(steps)) {
    if (!is.na(link_to[t]) & t < length(steps)) {
      for (nlink in link_to[t]) {
        wf <- append(
          wf,
          paste0(
            "  n", str_replace_all(steps[t], "\\.", "_"), " -> ",
            "  n", str_replace_all(nlink, "\\.", "_"), ";\n"
          )
        )
      }
    }
  }
  # add color, text
  wf <- c(wf, paste0(
    "  n", str_replace_all(steps, "\\.", "_"),
    ' [label=\"',
    steps, step_text, " ", ifelse(df_wf$no_run == 0, "", paste0(df_wf$no_success, "/", df_wf$no_run)),
    '\"', "fontcolor=", step_color,
    " color=white",
    "];\n"
  ))
  # end graph
  wf <- paste0(c(wf, "}"))
  return(wf)
}

#############################
## .change_branch function ##
#############################
## Internal parse function used in the plotWF function
.change_branch <- function(df_wf, wf) {
  long <- .find_long_branch(df_wf$t_number, df_wf$link_to)
  plot_start <- wf %>% str_which("digraph")
  sub_start <- wf %>% str_which("subgraph ")
  sub_steps_lines <- wf %>%
    str_which(" -> [^\\[]") %>%
    .[. > sub_start]
  sub_number <- wf[sub_steps_lines] %>%
    str_remove_all("[->;\n]") %>%
    str_remove("^.*[ ]+") %>%
    str_remove("n") %>%
    str_replace_all("_", "\\.") %>%
    sapply(function(x) df_wf$t_number[df_wf$t_number == x]) %>%
    unlist()
  move_line_num <- sub_steps_lines[!sub_number %in% df_wf$t_number[long]]
  if (length(move_line_num) == 0) {
    return(wf)
  }
  move_lines <- wf[move_line_num]
  wf <- wf %>%
    .[-move_line_num] %>%
    append(move_lines, after = plot_start)
  return(wf)
}

#############################
## .tryCatch function ##
#############################
.tryCatch <- function(x, file=NULL) {
  if(is.null(file)) file=tempfile()
  tryCatch(
    expr = {
      cat("## Output", append = TRUE, file = file, "\n")
      capture.output(out <- eval(parse(text = x), envir = globalenv()), file = file, append = TRUE)
      message("DONE")
      return(out)
    },
    error = function(e) {
      print(e)
      message("ERROR")
      return("Caught an error!")
    },
    warning = function(w) {
      print(w)
      message("WARNING")
      return("Caught an warning!")
    }
  )
}

## Usage:
# .tryCatch(x=codeList[[1]])


