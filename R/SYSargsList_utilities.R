####################################################################
## Functions to construct SYSargsList objects and other utilities ##
####################################################################

##########################
## SPRproject function ##
##########################
## Detection and creation of the logs directory of the project
## This function detects an existing project or creates a project structure on the path provide
SPRproject <- function(projPath = getwd(), data = "data", param = "param", results = "results",
                       logs.dir= ".SPRproject", sys.file=".SPRproject/SYSargsList.yml",
                       envir=new.env(), restart = FALSE, load.envir = FALSE,
                       overwrite=FALSE, silent = FALSE){
  if (!file.exists(projPath)) stop("Provide valid 'projPath' PATH.")
  ## Main folder structure
  dirProject <- .dirProject(projPath=projPath, data = data, param = param, results = results, silent = silent)
  ## sys.file full path
  sys.file <- file.path(projPath, sys.file)
  ## log Folder
  logs.dir <- file.path(projPath, logs.dir)
  if (file.exists(logs.dir) == FALSE) {
    dir.create(logs.dir, recursive = TRUE)
    if (silent != TRUE) cat("Creating directory '", normalizePath(logs.dir), "'", sep = "", "\n")
  } else if (file.exists(logs.dir) == TRUE) {
    if (restart == FALSE) {
      if (file.exists(logs.dir)) {
        if(overwrite==FALSE){
          stop(paste0(logs.dir, " already exists.", "\n", "The Workflow can be restart where it was stopped, using the argument 'restart=TRUE'."))
        } else if(overwrite==TRUE){
          unlink(logs.dir, recursive = TRUE)
          dir.create(logs.dir, recursive = TRUE)
          if (silent != TRUE) cat("Creating directory '", normalizePath(logs.dir), "'", sep = "", "\n")
        }
      }
    } else if (restart == TRUE) {
      if(!file.exists(sys.file)) stop("Provide valid 'sys.file' PATH")
      sal <- read_SYSargsList(sys.file)
      sal <- as(sal, "list")
      if(load.envir==TRUE){
        if(is.null(sal$projectInfo$envir)){
          message("We could not find any environment saved...")
        } else {
          envir <- readRDS(sal$projectInfo$envir)
          sal$runInfo$envir <- envir
        }
      } else {
        sal$runInfo$envir <- envir
        sal$projectInfo$envir <- NULL
      }
      sal <- as(sal, "SYSargsList")
      write_SYSargsList(sal, sys.file, silent=FALSE)
      return(sal)
    }
  }
  ## Return SYSargsList obj - empty
  yaml::write_yaml("", file= sys.file)
  dirProject <- c(dirProject, logsDir=logs.dir, sysargslist=normalizePath(sys.file))
  init <- as(SYScreate("SYSargsList"), "list")
  init$projectInfo <- dirProject
  init$runInfo <- list(env=envir)
  init <- as(init, "SYSargsList")
  write_SYSargsList(init, sys.file, silent=silent)
  return(init)
}

## Usage: 
# sal <- SPRproject()
# sal <- SPRproject(projPath = tempdir())
# sal <- SPRproject(restart = TRUE)
# sal <- SPRproject(overwrite = TRUE)

################################
## write_SYSargsList function ##
################################
write_SYSargsList <- function(args, sys.file=".SPRproject/SYSargsList.yml", silent=FALSE){
  if(!inherits(args, "SYSargsList")) stop("args needs to be object of class 'SYSargsList'.") 
  args2 <- sysargslist(args)
  args_comp <- sapply(args2, function(x) list(NULL))
  steps <- names(stepsWF(args))
  ## special case for "runInfo" slot
  yaml_slots <- c("runInfo")
  for(i in yaml_slots){
    for(j in steps){
    args_comp[[i]] <- yaml::as.yaml(args2[[i]]$directory[[j]])
    }
  }
 ## Simple yaml slots
  yaml_slots <- c("projectInfo", "SEobj")
  for(i in yaml_slots){
    args_comp[[i]] <- yaml::as.yaml(args2[[i]])
  }
  ## Yaml Slots + steps
  yaml_slots_S <- c("statusWF", "dependency","targets_connection")
  for(i in yaml_slots_S){
    steps_comp <- sapply(steps, function(x) list(NULL))
    for(j in steps){
      steps_comp[j] <- yaml::as.yaml(args2[[i]][j])
    }
    args_comp[[i]] <- steps_comp
  }
  ## DataFrame Slots
  df_slots <- c("targetsWF", "outfiles")
  for(i in df_slots){
    #  args_comp[[i]] <- yaml::as.yaml(as.data.frame(args2[[i]]$Mapping))
    steps_comp <- sapply(steps, function(x) list(NULL))
    for(j in steps){
      steps_comp[j] <- yaml::as.yaml(data.frame(args2[[i]][[j]], check.names=FALSE))
    }
    args_comp[[i]] <- steps_comp
  }
 ## SYSargs2 and LineWise
  steps_comp <- sapply(steps, function(x) list(NULL))
  for(j in steps){
    if(inherits(args2[["stepsWF"]][[j]], "SYSargs2")){
      step_obj <- sysargs2(args2[["stepsWF"]][[j]])
      steps_comp[[j]] <- yaml::as.yaml(step_obj)
    } else if(inherits(args2[["stepsWF"]][[j]], "LineWise")){
      step_obj <- linewise(args2[["stepsWF"]][[j]])
      step_obj$codeLine <- as.character(step_obj$codeLine)
      steps_comp[[j]] <- yaml::as.yaml(step_obj)
    }
  }
  args_comp[["stepsWF"]] <- steps_comp
  ## Save file
  yaml::write_yaml(args_comp, sys.file)
  if (silent != TRUE) cat("Creating file '", normalizePath(sys.file), "'", sep = "", "\n")
  return(normalizePath(sys.file))
}

# ## Usage: 
# write_SYSargsList(sal, sys.file, silent=FALSE)

################################
## read_SYSargsList function ##
################################
read_SYSargsList <- function(sys.file){
  args_comp_yml <- yaml::read_yaml(sys.file)
  args_comp <- sapply(args_comp_yml, function(x) list(NULL))
  steps <- names(args_comp_yml$stepsWF)
  ## Simple yaml slots
  yaml_slots <- c("projectInfo", "runInfo", "SEobj")
  for(i in yaml_slots){
    args_comp[[i]] <- yaml::yaml.load(args_comp_yml[i])
  }
  ## Yaml Slots + steps
  yaml_slots_S <- c("dependency","targets_connection")
  for(i in yaml_slots_S){
    steps_comp <- sapply(steps, function(x) list(NULL))
    for(j in steps){
      steps_comp[j] <- yaml::yaml.load(args_comp_yml[[i]][j])
    }
    args_comp[[i]] <- steps_comp
  }
  ## targetsWF Slots
  df_slots <- c("targetsWF", "outfiles")
  for(i in df_slots){
    steps_comp <- sapply(steps, function(x) list(NULL))
    for(j in steps){
      steps_comp[[j]] <- DataFrame(yaml::yaml.load(args_comp_yml[[i]][[j]]), check.names = FALSE)
    }
    args_comp[[i]] <- steps_comp
  }
  ## SYSargs2 and LineWise
  if(length(args_comp_yml$stepsWF)>=1){
    steps_comp <- sapply(steps, function(x) list(NULL))
    for(j in steps){
      if("codeLine" %in% names(yaml::yaml.load(args_comp_yml[["stepsWF"]][[j]]))){
        args <- yaml::yaml.load(args_comp_yml[["stepsWF"]][[j]])
        if(length(args$codeLine)>=1) { args$codeLine <- parse(text=args$codeLine)}
        if(length(args$codeChunkStart)==0) args$codeChunkStart <- integer()
        if(length(args$rmdPath)==0) args$rmdPath <- character()
        if(length(args$dependency)==0) args$dependency <- character()
        args <- as(args, "LineWise")
        steps_comp[[j]] <- args
      } else {
        args <- yaml::yaml.load(args_comp_yml[["stepsWF"]][[j]])
        args[["status"]][[2]] <- data.frame(args[["status"]][[2]], check.names = FALSE)
        args[["status"]][[3]] <- data.frame(args[["status"]][[3]])
        steps_comp[[j]] <- as(args, "SYSargs2")
      }
      args_comp[["stepsWF"]] <- steps_comp
    }
  } else if (length(args_comp_yml$stepsWF)>=0){
    args_comp[["stepsWF"]] <- list()
  }
  ## status
  yaml_slots_Status <- c("statusWF")
  for(i in yaml_slots_Status){
    steps_comp <- sapply(steps, function(x) list(NULL))
    for(j in steps){
      steps_comp[j] <- yaml::yaml.load(args_comp_yml[[i]][j])
      steps_comp[j][[1]][[2]] <- data.frame(steps_comp[j][[1]][[2]], check.names = FALSE)
      steps_comp[j][[1]][[3]] <- data.frame(steps_comp[j][[1]][[3]])
    }
    args_comp[[i]] <- steps_comp
  }
  
  return(as(args_comp,"SYSargsList"))
}

# ## Usage: 
# sys.file=".SPRproject/SYSargsList.yml"
# sal3 <- read_SYSargsList(sys.file)

################################
## .dirProject function ##
################################
.dirProject <- function(projPath, data, param, results, silent){
  project <- list(
    project = projPath, 
    data = file.path(projPath, data),
    param = file.path(projPath, param), 
    results = file.path(projPath, results)
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
        "There is no directory called", "\n", paste(names(create), collapse = " OR ", sep="\n"), "\n", "\n",
        "Would you like to create this directory now? Type a number: \n 1. Yes \n 2. No \n"
      ))
    } else {
      ## For an non-interactive session
      dir_create <- "1"
    }
    for (i in seq_along(create)) {
      if (dir_create == "1") {
        dir.create(create[[i]], recursive = TRUE)
        if (silent != TRUE) cat(paste("Creating directory:", create[[i]]), "\n")
      } else if (dir_create == 2) {
        stop("Aborting project creation.")
      }
    }
  }
  return(project)
}

########################################################
## SYSargsList ##
########################################################
SYSargsList <- function(args=NULL, step_name="default",
                        targets=NULL, wf_file=NULL, input_file=NULL, dir_path=".", inputvars=NULL,
                        rm_targets_col = NULL, 
                        dir=TRUE,
                        dependency=NULL, silent = FALSE) {
  sal <- as(SYScreate("SYSargsList"), "list")
  if(all(is.null(args) && is.null(wf_file) && is.null(input_file))){
    sal <- sal ## This will not create a SPRproject. 
    #message("please use 'SPRproject()' function") 
  } else if (!is.null(args)){
    if(inherits(args, "SYSargs2")){
      if(length(cmdlist(args))==0) stop ("Argument 'args' needs to be assigned an object of class 'SYSargs2' after 'renderWF()'.")
      # sal <- as(SYScreate("SYSargsList"), "list")
      if(step_name=="default"){
        step_name <- "Step_x"
      } else {
        step_name <- step_name
      }
      sal$stepsWF <- list(args)
      ## Targets
      if(length(targets(args))==0){
        sal$targetsWF <- list(NULL)
      } else {
        sal$targetsWF <- list(as(args, "DataFrame"))
      }
      ## Status
      if(length(status(args))==0){
        sal$statusWF <- list(.statusPending(args))
      } else {
        sal$statusWF <- list(status(args))
      }
      sal$dependency <- list(dependency)
      sal$outfiles <- list(.outList2DF(args))
      sal$targets_connection <- list(NULL)
      sal$runInfo <- list(directory=dir)
      names(sal$stepsWF) <- names(sal$targetsWF) <- names(sal$statusWF) <- names(sal$dependency) <- names(sal$outfiles) <- names(sal$targets_connection) <- names(sal$runInfo) <- step_name
    } else stop("Argument 'args' needs to be assigned an object of class 'SYSargs2'.")
  } else if(all(!is.null(wf_file) && !is.null(input_file))){
    ## targets
    if(is.null(targets)){
      targets <- targets
    } else if(inherits(targets, "character")){
        if(all(all(file.exists(targets)) && length(targets)==1)){
          targets <- targets
        } else {
          targets_step <- targets
          targets <- NULL
        }
    }
    WF <- loadWorkflow(targets=targets, wf_file=wf_file, 
                       input_file=input_file, dir_path=dir_path)
    WF <- renderWF(WF, inputvars=inputvars)
    if(step_name=="default"){
      step_name <- "Step_x"
    } else {
      step_name <- step_name
    }
    sal$stepsWF <- list(WF)
    names(sal$stepsWF) <- step_name
    ## Connection to previous outfiles
    if(exists("targets_step")){
      targets_step_list <- list(targets_step=targets_step)
      new_targets_col <- names(inputvars)
      if(is.null(new_targets_col)) 
        stop("inputvars argument need to be assigned to the output column names from the previous step specified on the targets argument")
      new_col <- list(new_targets_col=new_targets_col)
      if(!is.null(rm_targets_col)){
        rm_col <- list(rm_targets_col=rm_targets_col)
      } else {
        rm_col <- list(rm_targets_col=NULL)
      }
      sal$targets_connection <- list(list(targets_step=targets_step_list, new_targets_col=new_col, rm_targets_col=rm_col))
      names(sal$targets_connection) <- step_name
    }
    if(length(sal$targets_connection)==0){
    sal$targets_connection <- list(NULL)
    names(sal$targets_connection) <- step_name
    } 
    ## dependency
    if(is.null(dependency)){
      sal$dependency <- list(NA)
    } else {
      sal$dependency <- list(dependency)
    }
    ## statusWF
    sal$statusWF <- list(.statusPending(WF))
    # directory structure
    sal$runInfo <- list(directory=dir)
    ## names
    names(sal$statusWF) <- names(sal$dependency) <- names(sal$runInfo) <- step_name
    ## outfiles
    if(length(sal$stepsWF) > 0) {
      sal$outfiles <- .outList2DF(sal)
      ## targets
      if(length(targets(sal$stepsWF[[1]])) > 0 ) {
        sal$targetsWF <- list(as(sal$stepsWF[[1]], "DataFrame"))
      } else {
        sal$targetsWF <- list(S4Vectors::DataFrame())
      }
    }
    names(sal$targetsWF) <- step_name
  }
  sal <- as(sal, "SYSargsList")
  return(sal)
}

########################
## internal function ##
########################
.statusPending <- function(args){
  status.pending <- check.output(args)
  if(inherits(args, "SYSargsList")){
    for (i in seq_along(status.pending)){
      pending <- sapply(stepsWF(args)[[i]]$files$steps, function(x) list(x="Pending"))
      pending <- data.frame(matrix(unlist(pending), ncol = length(pending), byrow=TRUE), stringsAsFactors = TRUE)
      colnames(pending) <- stepsWF(args)[[i]]$files$steps
      status.pending[[i]] <- cbind(status.pending[[i]], pending)
    }
  } else if(inherits(args, "SYSargs2")){
    pending <- sapply(args$files$steps, function(x) list(x="Pending"))
    pending <- data.frame(matrix(unlist(pending), ncol = length(pending), byrow=TRUE), stringsAsFactors = TRUE)
    colnames(pending) <- args$files$steps
    status.pending <- cbind(status.pending, pending)
  }
  status.pending[c(2:4)] <- sapply(status.pending[c(2:4)], as.numeric)
  pendingList <- list(status.summary= .statusSummary(status.pending), 
                      status.completed=status.pending, status.time=data.frame())
  return(pendingList)
}

.outList2DF <- function(args) {
  if (inherits(args, "list")) {
    args <- as(args, "SYSargsList")
    out <- sapply(names(stepsWF(args)), function(x) list(NULL))
    for (i in seq_along(stepsWF(args))) {
      l_out <- output(stepsWF(args)[[i]])
      out[[i]] <- S4Vectors::DataFrame(matrix(unlist(l_out), nrow = length(l_out), byrow = TRUE))
      colnames(out[[i]]) <- stepsWF(args)[[i]]$files$output_names
    }
  } else if (inherits(args, "SYSargs2")) {
    l_out <- output(args)
    out <- S4Vectors::DataFrame(matrix(unlist(l_out), nrow = length(l_out), byrow = TRUE))
    colnames(out) <- args$files$output_names
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
runWF <- function(args, force=FALSE, saveEnv=TRUE, 
                  warning.stop=FALSE, error.stop=TRUE, silent=FALSE, ...) {
  # Validations
  if (!inherits(args, "SYSargsList")) stop("Argument 'args' needs to be assigned an object of class 'SYSargsList'")
  if (is.null(projectInfo(args)$project)) stop("Project was not initialized with the 'SPRproject' function.")
  sysproj <- projectInfo(args)$logsDir
  ## check dependency
  for(i in seq_along(dependency(args))){
    if(all(!is.na(dependency(args)[[i]]))){
      dep_names <- unlist(dependency(args)[[i]])
      if(any(!dep_names %in% names(stepsWF(args))))
        stop("'args' has dependency on the following steps:", "\n",
             paste0("      ", paste0(dep_names, collapse = " AND ")), "\n", 
             "Please make sure that this step is present.")
    }
  }
  ## Logs
  file_log <- file.path(sysproj, paste0("_logWF_", format(Sys.time(), "%b%d%Y_%H%M")))
  args[["projectInfo"]]$logsFile <- file_log
  ## steps loop
  args2 <- args
  for (i in seq_along(stepsWF(args2))){
    cat("# ", names(stepsWF(args2)[i]), "\n", file=file_log, fill=TRUE, append=TRUE)
    ## SYSargs2 STEP
    if(inherits(stepsWF(args2)[[i]], "SYSargs2")){
      step <- names(stepsWF(args2)[i])
      cat(crayon::bgMagenta(paste0("Running Step: ", step)), "\n")
      args.run <- stepsWF(args2)[[i]]
      ## runC arguments
      dir <- args2$runInfo$directory[[i]]
      dir.name <- step
      args.run <- runCommandline(args.run, dir = dir, dir.name = dir.name, force=force)
      cat(readLines(args.run$files$log), file=file_log, sep = "\n", append=TRUE)
      ## update object
      step.status.summary <- status(args.run)$status.summary
      statusWF(args2, i) <- args.run$status
      stepsWF(args2, i) <- args.run
      args2[["outfiles"]][[i]] <- .outList2DF(args.run)
      args2 <- .updateAfterRunC(args2, step)
      assign(x=as.character(as.list(match.call())$args), args2, envir = args2$runInfo$env)
      ## Stop workflow
      if(is.element("Warning", unlist(step.status.summary))){
        if(warning.stop==TRUE) {
          on.exit(return(args2))
          stop("Caught an warning, stop workflow!")
           }
      } else if(is.element("Error", unlist(step.status.summary))){
        if(error.stop==TRUE) {
          on.exit(return(args2))
          stop("Caught an error, stop workflow!")
        }
      }
      cat(crayon::blue(paste0("Step Status: ", step.status.summary), "\n"))
    } else if(inherits(stepsWF(args2)[[i]], "LineWise")){
      step <- names(stepsWF(args2)[i])
      cat(crayon::bgMagenta(paste0("Running Step: ", step)), "\n")
      args.run <- stepsWF(args2)[[i]]
      envir <- args2$runInfo$env
      args.run <- runRcode(args.run, step, file_log, envir, force=force)
      stepsWF(args2, i) <- args.run
      statusWF(args2, i) <- args.run$status
      ## Stop workflow
      if(is.element("Warning", unlist(args.run$status$status.summary))){
        if(warning.stop==TRUE) {
          on.exit(return(args2))
          stop("Caught an warning, stop workflow!")
        }
      } else if(is.element("Error", unlist(args.run$status$status.summary))){
        if(error.stop==TRUE) {
          on.exit(return(args2))
          stop("Caught an error, stop workflow!")
        }
      }
      cat(crayon::blue(paste0("Step Status: ", args.run$status$status.summary), "\n"))
    }
  }
  if(saveEnv==TRUE){
    envPath <- file.path(sysproj, "sysargsEnv.rds")
    saveRDS(args2$runInfo$env, envPath)
    args2[["projectInfo"]][["envir"]] <- envPath
  }
  args2 <- .check_write_SYSargsList(args2)
  return(args2)
}


runRcode <- function(args.run, step, file_log, envir, force=FALSE){
  ## Validation for 'args.run'
  if(!inherits(args.run, "LineWise")) stop("Argument 'args.run' needs to be assigned an object of class 'LineWise'")
  pb <- txtProgressBar(min = 0, max = length(args.run), style = 3)
  ## log_file
  cat(c(
    paste0("Time: ", paste0(format(Sys.time(), "%b%d%Y_%H%Ms%S"))), "\n",
    "## Code: ",
    "```{r, eval=FALSE} ",
    capture.output(codeLine(args.run)),
    "```", "\n",
    "## Stdout: ",
    "```{r, eval=FALSE}" ), file = file_log, sep = "\n", append = TRUE)
  ## Check status of step
  if(all("Success" %in% status(args.run)[[step]]$status.summary && force==FALSE)){
    cat("The step status is 'Success' and it was skipped.", file=file_log, fill=TRUE, append=TRUE)
  } else {
    ## Status and time register
    step_status <- list()
    time_status <- data.frame(Step=step, time_start=NA, time_end=NA)
    time_status$time_start <- Sys.time()
    ## Running the code
    stdout <- .tryRcode(args.run$codeLine, envir = envir)
    ## save stdout to file
    capture.output(stdout$stdout, file=file_log, append=TRUE)
    ## save error and warning messages
    if(!is.null(stdout$error)) {
      cat("## Error", file=file_log, sep = "\n", append=TRUE)
      cat(stdout$error, file=file_log, sep = "\n", append=TRUE)
      step_status[["status.summary"]] <- "Error"
    } else if(!is.null(stdout$warning)) {
      cat("## Warning", file=file_log, sep = "\n", append=TRUE)
      cat(stdout$warning, file=file_log, sep = "\n", append=TRUE)
      step_status[["status.summary"]] <- "Warning"
    } else if(all(is.null(stdout$error) && is.null(stdout$warning))){
      step_status[["status.summary"]] <- "Success"
    }
    ## Saving the new status
    step_status[["status.completed"]] <- data.frame(Step=step, Status=step_status[[1]])
    time_status$time_end <- Sys.time()
    step_status[["status.time"]] <- time_status
    args.run[["status"]] <- step_status
  }
  setTxtProgressBar(pb, length(args.run))
  # close R chunk
  cat("```", file=file_log, sep = "\n", append=TRUE)
  close(pb) 
  return(args.run)
}

.tryRcode <- function(command, envir){
  warning <- error <- NULL
  value <- withCallingHandlers(
    tryCatch(
      eval(command, envir = envir), 
      error = function(e) {
        error <<- conditionMessage(e)
        NULL
      }), 
    warning = function(w) {
      warning <<- append(warning, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(stdout = value, warning = warning, error = error)
}

.updateAfterRunC <- function(args, step){
  # conList <- lapply(args$targets_connection, function(x) if(!is.null(x)){
  #   x$targets_step })
  conList <- args$targets_connection[lengths(args$targets_connection) != 0]
  conList_step <- sapply(conList, "[[", 1)
  for(l in seq_along(conList_step)){
    if(step %in% conList_step[[l]]){
      requiredUP <- names(conList)[[l]]
      for(s in requiredUP){
        WF <- args[s]
        WFstep <- names(stepsWF(WF))
        
        new_targets <- WF$targetsWF[[1]]
        #targets_name <- paste(colnames(targets_WF), collapse="|")
        col_out <- lapply(outfiles(args), function(x) colnames(x))
        col_out <- col_out[col_out %in% WF$targets_connection[[WFstep]]$new_targets_col[[1]]]
        col_out_df <- lapply(names(col_out), function(x) getColumn(args, step=x, df = "outfiles", column = col_out[[x]]))
        names(col_out_df) <- col_out
        new_targets[as.character(col_out)] <- col_out_df
# 
#         WF$targets_connection[[WFstep]]$new_targets_col[[1]]
#         
#         targesCon <- args$targets_connection[names(args$targets_connection)==WFstep][[1]]
#         targets_name <- paste(colnames(targetsWF(args)[step][[1]]), collapse="|")
#         new_targets_col <- targesCon[[2]][[1]][-c(which(grepl(targets_name, targesCon[[2]][[1]])))]
#         if(is.null(targesCon[[3]][[1]])){
#           old_targets <- args$targetsWF[[step]]
#         } else {
#           old_targets <- args$targetsWF[[step]][-c(which(grepl(paste(targesCon[[3]][[1]], collapse="|"), colnames(args$targetsWF[[step]]))))]
#         }
#         new_targets <- cbind(args$outfiles[[step]][new_targets_col], old_targets)
        WF2 <- stepsWF(WF)[[1]]
        WF2 <- updateWF(WF2, new_targets= targets.as.list(data.frame(new_targets)), inputvars=WF2$inputvars, write.yaml = FALSE)
        ## Preserve outfiles
        WF2[["output"]] <- WF$stepsWF[[s]]$output
        args <- sysargslist(args)
        args$stepsWF[[WFstep]] <- WF2
        args$targetsWF[[WFstep]] <- as(WF2, "DataFrame")
        args$outfiles[[WFstep]] <- output.as.df(WF2)
        args$statusWF[[WFstep]] <- WF2$status
        args <- as(args, "SYSargsList")
      }
    }
  }
  return(args)
}

.statusSummary <- function(args){
  if(inherits(args, "SYSargs2")) {
    step.status.summary <- args$status$status
  } else if(inherits(args, "data.frame")){
    step.status.summary <- args[5:ncol(args)]
    } else stop("Argument 'args' needs to be assigned an object of class 'SYSargs2' or 'data.frame'.")
  if ("Error" %in% unlist(unique(step.status.summary))){
    step.status <- "Error"
  } else  if ("Warning" %in% unlist(unique(step.status.summary))){
    step.status <- "Warning"
  } else if ("Success" %in% unlist(unique(step.status.summary))){
    step.status <- "Success"
  } else if ("Pending" %in% unlist(unique(step.status.summary))){
    step.status <- "Pending"
  } else if (is.null(step.status.summary)){
    step.status <- "Pending"
  }
  return(step.status)
}

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
  sysargslist$projectInfo[["Report"]] <- normalizePath(paste0(file_path, "/", file_out, ".", ext))
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
  file <-  projectInfo(sysargslist)$logsFile
  if (!file.exists(file) == TRUE) stop("Provide valid 'sysargslist' object. Check the initialization of the project.")
  if(fileName=="default"){
    fileName <- file.path(projectInfo(sysargslist)$project, paste0("logs_", format(Sys.time(), "%b%d%Y_%H%M"),".Rmd"))
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
  sysargslist$projectInfo[["Report_Logs"]] <- file.path(file_path, paste(file_out, ext, sep="."))
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
    dir_path <- .getPath(files(input_file)[['yml']])
    if(is.na(files(input_file)[['targets']])) {
      targets <- NULL
    } else {
      targets <- targets(input_file)
    }
    args1 <- loadWorkflow(
      targets = targets, wf_file = basename(files(input_file)[['cwl']]), 
      input_file = basename(files(input_file)[['yml']]), dir_path = dir_path
    )
    args1 <- as(args1, "list")
    args1$yamlinput <- input
    if("ModulesToLoad" %in% names(param))
      for (i in seq_along(param$ModulesToLoad)){
        args1$modules[names(param$ModulesToLoad[i])] <- param$ModulesToLoad[[i]]
      }
    args1 <- out_obj <- as(args1, "SYSargs2")
    out_msg <- c("yamlinput(args1)")
    path_file <- files(input_file)[['yml']]
  } else if (class(input_file) == "SYSargsList") { ##TODO
    input <- out_obj <- .replace(input = input_file$sysconfig, param = param)
    path_file <- input_file$projectInfo$project
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
    args1$files$yml <- path
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
      files = list(),
      inputvars = list(), 
      cmdToCwl = list(),
      status = data.frame()
    )
    return(as(SYS.empty, "SYSargs2"))
  } else if (class == "SYSargsList") {
    SYS.empty <- list(
      stepsWF = list(),
      statusWF = list(),
      targetsWF = list(), 
      outfiles = list(),
      SEobj = list(), 
      dependency = list(),
      targets_connection=list(),
      projectInfo = list(),
      runInfo = list()
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
  return(output)
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
#  if (!file.exists(x)) warning("No such file or directory. Check the file PATH.")
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

# ##############################
# ## .sysconfigCheck function ##
# ##############################
# ## Function to check with the paths provided on the sysconfig file exists.
# .sprconfigCheck <- function(sysconfig) {
#   if (!file.exists(sysconfig) == TRUE) stop("Provide valid 'sysconfig' file. Check the file PATH.")
#   sysconfig <- yaml::read_yaml(sysconfig, eval.expr = TRUE)
#   project <- list(
#     project = sysconfig$project$path, data = sysconfig$data$path, param = sysconfig$param$path,
#     results = sysconfig$results$path
#   )
#   for (i in seq_along(project)) {
#     if (is.null(project[[i]])) {
#       warning(paste(names(project[i]), "directory is missing..."))
#     } else {
#       tryPath(path = project[[i]])
#     }
#   }
#   if (!is.null(sysconfig$targets$path)) tryPath(path = sysconfig$targets$path)
#   if (is.null(sysconfig$script)) warning("Workflow is missing...")
#   if (!is.null(sysconfig$script$path)) tryPath(path = sysconfig$script$path)
# }
# 
# ## Usage:
# # .sysconfigCheck(sysconfig="SYSconfig.yml")

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


