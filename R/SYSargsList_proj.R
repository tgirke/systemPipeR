##########################
## SPRproject function ##
##########################
## Detection and creation of the logs directory of the project
## This function detects an existing project or creates a project structure on the path provide

SPRproject <- function(projPath = getwd(), data = "data", param = "param", results = "results",
                      logs.dir= ".SPRproject",
                       restart = FALSE, sys.file=".SPRproject/SYSargsList.yml",
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
            return(sal)
        }
    }
    ## Return SYSargsList obj - empty
    yaml::write_yaml("", file= sys.file)
    dirProject <- c(dirProject, logs=logs.dir, sysargslist=normalizePath(sys.file))
    project_path <- list(class = "Directory", path = normalizePath(dirProject$project))
    data_path <- list(class = "Directory", path = normalizePath(dirProject$data))
    param_path <- list(class = "Directory", path = normalizePath(dirProject$param))
    results_path <- list(class = "Directory", path = normalizePath(dirProject$results))
    logs_path <- list(class = "Directory", path = normalizePath(dirProject$logs))
    sysargslist_path <- list(class = "File", path = normalizePath(sys.file))
    init <- as(SYScreate("SYSargsList"), "list")
    initProj <- list(project = project_path, data = data_path,
                     param = param_path, results = results_path, logs_path = logs_path, 
                     sysargslist=sysargslist_path)
    init$projectWF <- dirProject
    init$sprconfig <- initProj
    init <- as(init, "SYSargsList")
    write_SYSargsList(init, sys.file, silent=silent)
    return(init)
}



# projPath = getwd()
# data = "data"
# param = "param"
# results = "results"
# logs.dir= ".SPRproject"
# restart = FALSE
# sys.file=".SPRproject/SYSargsList.yml"
# silent = FALSE

##########################
## write_SYSargsList function ##
##########################
write_SYSargsList <- function(WF, sys.file=".SPRproject/SYSargsList.yml", silent=FALSE){
    if(!inherits(WF, "SYSargsList")) stop("args needs to be object of class 'SYSargsList'.") 
    WF2 <- sysargslist(WF)
    WF_comp <- sapply(WF2, function(x) list(NULL))
    steps <- names(stepsWF(WF))
    ## Simple yaml slots
    yaml_slots <- c("sprconfig", "projectWF", "SEobj")
    for(i in yaml_slots){
        WF_comp[[i]] <- yaml::as.yaml(WF2[[i]])
    }
    ## Yaml Slots + steps
    yaml_slots_S <- c("statusWF", "dependency","targets_connection")
    for(i in yaml_slots_S){
        steps_comp <- sapply(steps, function(x) list(NULL))
        for(j in steps){
            steps_comp[j] <- yaml::as.yaml(WF2[[i]][j])
        }
        WF_comp[[i]] <- steps_comp
    }
    ## DataFrame Slots
    df_slots <- c("targetsWF", "outfiles")
    for(i in df_slots){
      #  WF_comp[[i]] <- yaml::as.yaml(as.data.frame(WF2[[i]]$Mapping))
        steps_comp <- sapply(steps, function(x) list(NULL))
        for(j in steps){
            steps_comp[j] <- yaml::as.yaml(data.frame(WF2[[i]][[j]], check.names=FALSE))
        }
        WF_comp[[i]] <- steps_comp
    }
    ## SYSargs2 and LineWise
    steps_comp <- sapply(steps, function(x) list(NULL))
    for(j in steps){
        if(inherits(WF2[["stepsWF"]][[j]], "SYSargs2")){
            step_obj <- sysargs2(WF2[["stepsWF"]][[j]])
            steps_comp[[j]] <- yaml::as.yaml(step_obj)
        } else if(inherits(WF2[["stepsWF"]][[j]], "LineWise")){
            step_obj <- linewise(WF2[["stepsWF"]][[j]])
            step_obj$codeLine <- as.character(step_obj$codeLine)
            steps_comp[[j]] <- yaml::as.yaml(step_obj)
        }
    }
    WF_comp[["stepsWF"]] <- steps_comp
    ## Save file
    yaml::write_yaml(WF_comp, sys.file)
    if (silent != TRUE) cat("Creating file '", normalizePath(sys.file), "'", sep = "", "\n")
    return(normalizePath(sys.file))
}
# 
# write_SYSargsList(sal, sys.file, silent=FALSE)

read_SYSargsList <- function(sys.file){
    WF_comp_yml <- yaml::read_yaml(sys.file)
    WF_comp <- sapply(WF_comp_yml, function(x) list(NULL))
    steps <- names(WF_comp_yml$stepsWF)
    ## Simple yaml slots
    yaml_slots <- c("sprconfig", "projectWF", "SEobj")
    for(i in yaml_slots){
        WF_comp[[i]] <- yaml::yaml.load(WF_comp_yml[i])
    }
    ## Yaml Slots + steps
    yaml_slots_S <- c("statusWF", "dependency","targets_connection")
    for(i in yaml_slots_S){
        steps_comp <- sapply(steps, function(x) list(NULL))
        for(j in steps){
            steps_comp[j] <- yaml::yaml.load(WF_comp_yml[[i]][j])
        }
        WF_comp[[i]] <- steps_comp
    }
    ## DataFrame Slots
    df_slots <- c("targetsWF", "outfiles")
    for(i in df_slots){
        steps_comp <- sapply(steps, function(x) list(NULL))
        for(j in steps){
            steps_comp[[j]] <- DataFrame(yaml::yaml.load(WF_comp_yml[[i]][[j]]))
        }
        WF_comp[[i]] <- steps_comp
    }
    ## SYSargs2 and LineWise
    if(length(WF_comp_yml$stepsWF)>=1){
      steps_comp <- sapply(steps, function(x) list(NULL))
      for(j in steps){
        if("codeLine" %in% names(yaml::yaml.load(WF_comp_yml[["stepsWF"]][[j]]))){
          args <- yaml::yaml.load(WF_comp_yml[["stepsWF"]][[j]])
          if(length(args$codeLine)>=1) { args$codeLine <- parse(text=args$codeLine)}
          if(length(args$codeChunkStart)==0) args$codeChunkStart <- integer()
          if(length(args$rmdPath)==0) args$rmdPath <- character()
          if(length(args$dependency)==0) args$dependency <- character()
          args <- as(args, "LineWise")
          steps_comp[[j]] <- args
        } else {
          args <- as(yaml::yaml.load(WF_comp_yml[["stepsWF"]][[j]]), "SYSargs2")
          steps_comp[[j]] <- args
        }
        WF_comp[["stepsWF"]] <- steps_comp
      }
    } else if (length(WF_comp_yml$stepsWF)>=0){
      WF_comp[["stepsWF"]] <- list()
    }
    return(as(WF_comp,"SYSargsList"))
}

# sys.file=".SPRproject/SYSargsList.yml"
# sal3 <- read_SYSargsList(sys.file)


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

