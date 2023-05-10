####################################################################
## Functions to construct SYSargsList objects and other utilities ##
####################################################################

##########################
## SPRproject function ##
##########################
## Detection and creation of the logs directory of the project
## This function detects an existing project or creates a project structure on the path provide
SPRproject <- function(projPath = getwd(), data = "data", param = "param", results = "results",
                       logs.dir = ".SPRproject", sys.file = "SYSargsList.yml",
                       envir = new.env(),
                       restart = FALSE, resume = FALSE,
                       load.envir = FALSE,
                       overwrite = FALSE, silent = FALSE) {
    if (!file.exists(projPath)) stop("Provide valid 'projPath' PATH.")
    projPath <- normalizePath(projPath)
    on.exit(options(projPath = projPath), add = TRUE)
    ## Main folder structure
    dirProject <- .dirProject(projPath = projPath, data = data, param = param, results = results, silent = silent)
    ## sys.file full path
    sys.file <- file.path(logs.dir, sys.file)
    ## log Folder
    logs.dir <- file.path(logs.dir)
    if (file.exists(file.path(projPath, logs.dir)) == FALSE) {
        dir.create(file.path(projPath, logs.dir), recursive = TRUE)
        if (silent != TRUE) cat("Creating directory '", file.path(projPath, logs.dir), "'", sep = "", "\n")
        ## Return SYSargsList obj - empty
        yaml::write_yaml("", file = file.path(projPath, sys.file))
        dirProject <- c(dirProject, logsDir = logs.dir, sysargslist = sys.file)
        init <- as(SYScreate("SYSargsList"), "list")
        init$projectInfo <- dirProject
        init$runInfo <- list(env = envir)
        init <- as(init, "SYSargsList")
        ## used in `importWF`
        options(linewise_importing = FALSE)
    } else if (file.exists(file.path(projPath, logs.dir)) == TRUE) {
        ## reload object from sysarglist
        if (is.null(sys.file)) stop("Provide valid 'sys.file' PATH for restart/resume the project.")
        if (!file.exists(file.path(projPath, sys.file))) stop("Provide valid 'sys.file' PATH for restart/resume the project.")
        restart.sys.file <- normalizePath(file.path(projPath, sys.file))
        ## restart OR resume OR overwrite options...
        if (all(resume == TRUE && restart == TRUE && overwrite == TRUE)) {
            stop("Please select only one action:
        'resume' OR 'restart' OR 'overwrite'")
        } else if (all(resume == FALSE && restart == FALSE && overwrite == FALSE)) {
            stop(
                file.path(projPath, logs.dir),
                " already exists.", "\n",
                "The Workflow can be resume where it was stopped, using the argument 'resume=TRUE'",
                "\n", "OR ", "restart using 'restart=TRUE'. The latter will recover the workflow object, however, it will delete all the existing logs."
            )
        } else if (all(overwrite == TRUE && resume == FALSE && restart == FALSE)) {
            ## Return SYSargsList obj - empty
            unlink(logs.dir, recursive = TRUE)
            dir.create(file.path(projPath, logs.dir), recursive = TRUE)
            if (silent != TRUE) cat("Recreating directory '", file.path(projPath, logs.dir), "'", sep = "", "\n")
            yaml::write_yaml("", file = file.path(projPath, sys.file))
            dirProject <- c(dirProject, logsDir = logs.dir, sysargslist = sys.file)
            init <- as(SYScreate("SYSargsList"), "list")
            init$projectInfo <- dirProject
            init$runInfo <- list(env = envir)
            init <- as(init, "SYSargsList")
        } else if (all(overwrite == FALSE && any(resume, restart))) {
            ## read SYSargsList object
            init <- read_SYSargsList(restart.sys.file)
            init <- as(init, "list")
            if (load.envir == TRUE) {
                if (is.null(init$projectInfo$envir)) {
                    message("We could not find any environment saved...")
                } else {
                    envir <- readRDS(init$projectInfo$envir)
                    init$runInfo <- append(init$runInfo, envir, after = 0)
                    names(init$runInfo)[[1]] <- c("env")
                }
            } else {
                init$runInfo <- append(init$runInfo, envir, after = 0)
                names(init$runInfo)[[1]] <- "env"
                init$projectInfo$envir <- NULL
            }
            init$projectInfo$project <- projPath
            init$projectInfo$logsDir <- logs.dir
            init$projectInfo$sysargslist <- sys.file
            init <- as(init, "SYSargsList")
            if (all(resume == FALSE && restart == TRUE)) {
                ## remove old logs and create a empty folder
                unlink(file.path(projPath, logs.dir), recursive = TRUE)
                dir.create(file.path(projPath, logs.dir), recursive = TRUE)
                if (silent != TRUE) cat("Creating directory '", file.path(projPath, logs.dir), "'", sep = "", "\n")
            }
        } else {
            stop("Please select only one action:
        'resume' OR 'restart' OR 'overwrite'")
        }
    }
    write_SYSargsList(init, file.path(projPath, sys.file), silent = silent)
    ## Message about the paths
    if (normalizePath(getwd()) != projPath) {
        message(
            "Your current working directory is different from the directory chosen for the Project Workflow.",
            "\n",
            "For accurate location of the files and running the Workflow, please set the working directory to ",
            "\n",
            "'setwd('", projPath, "')'"
        )
    }
    return(init)
}

## Usage:
# sal <- SPRproject()
# sal <- SPRproject(projPath = tempdir())
# sal <- importWF(sal,
#                 file_path = system.file("extdata", "spr_simple_wf.Rmd", package = "systemPipeR"),
#                 verbose = FALSE)
# sal <- SPRproject(resume = TRUE)
# sal <- SPRproject(restart = TRUE)
# sal <- SPRproject(overwrite = TRUE)

##########################
## SYSargsList function ##
##########################
SYSargsList <- function(sysargs = NULL, step_name = "default",
                        targets = NULL, wf_file = NULL, input_file = NULL, dir_path = ".",
                        id = "SampleName",
                        inputvars = NULL,
                        rm_targets_col = NULL,
                        dir = TRUE,
                        dependency = NA,
                        run_step = "mandatory",
                        run_session = "management",
                        run_remote_resources = NULL,
                        silent = FALSE, projPath = getOption("projPath", getwd())) {
    ## step_name and dependency from importWF
    on.exit({
        options(replace_step = FALSE)
        options(spr_importing = FALSE)
    })
    ## check options
    run_step <- match.arg(run_step, c("mandatory", "optional"))
    run_session <- match.arg(run_session, c("management", "compute"))
    if (!is.null(run_remote_resources)) {
        if (!inherits(run_remote_resources, "list")) {
            stop("Argument 'run_remote_resources' needs to be assigned an object of class 'list'")
        }
    }
    ## sal object
    sal <- as(SYScreate("SYSargsList"), "list")
    ## Empty container
    if (all(is.null(sysargs) && is.null(wf_file) && is.null(input_file))) {
        sal <- sal ## This will not create a SPRproject.
        message("Please consider to use 'SPRproject()' function first")
        ## sal container based on a SYSargs2 container ##
    } else if (!is.null(sysargs)) {
        if (inherits(sysargs, "SYSargs2")) {
            if (length(cmdlist(sysargs)) == 0) stop("Argument 'sysargs' needs to be assigned an object of class 'SYSargs2' after 'renderWF()'.")
            if (step_name == "default") {
                step_name <- "Step_x"
            } else {
                .checkSpecialChar(step_name)
                step_name <- step_name
            }
            sal$stepsWF <- list(sysargs)
            ## Targets
            if (length(targets(sysargs)) == 0) {
                sal$targetsWF <- list(NULL)
                sal$SE <- list()
            } else {
                sal$targetsWF <- list(as(sysargs, "DataFrame"))
                row.names(sal$targetsWF[[1]]) <- sal$targetsWF[[1]][, sysargs$files$id]
                sal$SE <- list(SummarizedExperiment::SummarizedExperiment(
                    colData = sal$targetsWF,
                    metadata = sysargs$targetsheader
                ))
            }
            ## Status
            if (length(status(sysargs)) == 0) {
                sal$statusWF <- list(.statusPending(sysargs))
            } else {
                sal$statusWF <- list(status(sysargs))
            }
            sal$dependency <- list(dependency)
            sal$outfiles <- list(.outList2DF(sysargs))
            sal$targets_connection <- list(NULL)
            sal$runInfo <- list(runOption = list(list(
                directory = dir,
                run_step = run_step,
                run_session = run_session
            )))
            if (!is.null(run_remote_resources)) {
                if (run_session == "management") {
                    message("Please note that the '", step_name, "' run_session option '", run_session, "' was replaced with 'compute' because run_remote_resources was available.")
                    run_session <- "compute"
                }
            }
            names(sal$stepsWF) <- names(sal$targetsWF) <- names(sal$statusWF) <- names(sal$dependency) <- names(sal$outfiles) <- names(sal$targets_connection) <- names(sal$runInfo$runOption) <- names(sal$SE) <- step_name
        } else {
            stop("Argument 'sysargs' needs to be assigned an object of class 'SYSargs2'.")
        }
        ## Build the instance from files ##
    } else if (all(!is.null(wf_file) && !is.null(input_file))) {
        ## targets
        if (is.null(targets)) {
            targets <- targets
        } else if (inherits(targets, "character")) {
            if (all(is.fullPath(targets))) {
                targets <- targets
            } else if (all(all(file.exists(file.path(projPath, targets))) && length(targets) == 1)) {
                targets <- file.path(projPath, targets)
            } else { ## connection with previous steps
                targets_step <- targets
                targets <- NULL
            }
        } else if (inherits(targets, "SummarizedExperiment")) {
            se <- targets
            if (sum(dim(colData(targets))) == 0) {
                targets <- NULL
            } else {
                targets <- targets
            }
        } else {
            stop("Argument 'targets' needs to be assigned an object of class 'SummarizedExperiment', 'NULL' or 'character' PATH to tabular file.")
        }
        ## check param path
        if (is.fullPath(dir_path)) {
            dir_path <- dir_path
        } else {
            dir_path <- file.path(projPath, dir_path)
        }
        WF <- loadWF(
            targets = targets, wf_file = wf_file,
            input_file = input_file,
            dir_path = dir_path, id = id
        )
        ## Relative Path for targets and dir_path when these files are in the project folder,
        ## otherwise, keep the full path
        ## targets_path to projPath
        if (!is.null(WF$files$targets)) {
            if (!grepl(projPath, WF$files$targets)) {
                if (!is.na(WF$files$targets)) WF@files$targets <- gsub(projPath, "", WF$files$targets)
                if (all(!is.fullPath(WF$files$targets) && grepl("^/", WF$files$targets))) WF@files$targets <- sub("^(/|[A-Za-z]:|\\\\|~)", "", WF$files$targets)
            }
        }
        if (!grepl(projPath, WF$files$dir_path)) {
            ## dir_path
            if (!is.na(WF$files$dir_path)) WF@files$dir_path <- gsub(projPath, "", WF$files$dir_path)
            if (all(!is.fullPath(dir_path) && grepl("^/", WF$files$dir_path))) WF@files$dir_path <- sub("^(/|[A-Za-z]:|\\\\|~)", "", WF$files$dir_path)
        }
        WF <- renderWF(WF, inputvars = inputvars)
        if (step_name == "default") {
            step_name <- "Step_x"
        } else {
            .checkSpecialChar(step_name)
            step_name <- step_name
        }
        sal$stepsWF <- list(WF)
        names(sal$stepsWF) <- step_name
        ## Connection to previous outfiles
        if (exists("targets_step")) {
            targets_step_list <- list(targets_step = targets_step)
            new_targets_col <- names(inputvars)
            if (is.null(new_targets_col)) {
                stop("inputvars argument need to be assigned to the output column names from the previous step specified on the targets argument")
            }
            new_col <- list(new_targets_col = new_targets_col)
            if (!is.null(rm_targets_col)) {
                rm_col <- list(rm_targets_col = rm_targets_col)
            } else {
                rm_col <- list(rm_targets_col = NULL)
            }
            sal$targets_connection <- list(list(targets_step = targets_step_list, new_targets_col = new_col, rm_targets_col = rm_col))
            names(sal$targets_connection) <- step_name
        }
        if (length(sal$targets_connection) == 0) {
            sal$targets_connection <- list(NULL)
            names(sal$targets_connection) <- step_name
        }
        sal$dependency <- list(dependency)
        sal$statusWF <- list(.statusPending(WF))
        sal$runInfo <- list(runOption = list(list(
            directory = dir,
            run_step = run_step,
            run_session = run_session
        )))
        if (!is.null(run_remote_resources)) {
            if (run_session == "management") {
                message("Please note that the '", step_name, "' run_session option '", run_session, "' was replaced with 'compute' because run_remote_resources was available.")
                run_session <- "compute"
            }
            sal$runInfo <- list(runOption = list(list(
                directory = dir,
                run_step = run_step,
                run_session = run_session,
                run_remote_resources = run_remote_resources
            )))
        }
        ## names
        names(sal$statusWF) <- names(sal$dependency) <- names(sal$runInfo$runOption) <- step_name
        ## outfiles
        if (length(sal$stepsWF) > 0) {
            sal$outfiles <- .outList2DF(sal)
            ## targets
            if (length(targets(sal$stepsWF[[1]])) > 0) {
                sal$targetsWF <- list(as(sal$stepsWF[[1]], "DataFrame"))
                rownames(sal$targetsWF[[1]]) <- sal$targetsWF[[1]][, sal$stepsWF[[1]]$files$id]
                rownames(sal$outfiles[[1]]) <- sal$targetsWF[[1]][, sal$stepsWF[[1]]$files$id]
                if (exists("se", inherits = FALSE)) {
                    colData(se) <- sal$targetsWF[[1]]
                    metadata(se) <- list(metadata = metadata(se), SPRversion = utils::packageVersion("systemPipeR"), targetsheader = sal$stepsWF[[1]]$targetsheader)
                    sal$SE <- list(se)
                } else {
                    sal$SE <- list(SummarizedExperiment::SummarizedExperiment(
                        colData = sal$targetsWF[[1]],
                        metadata = list(SPRversion = utils::packageVersion("systemPipeR"), targetsheader = sal$stepsWF[[1]]$targetsheader)
                    ))
                }
            } else {
                sal$targetsWF <- list(S4Vectors::DataFrame())
                if (exists("se", inherits = FALSE)) {
                    sal$SE <- list(se)
                } else {
                    sal$SE <- list(SummarizedExperiment::SummarizedExperiment())
                }
            }
        }
        names(sal$targetsWF) <- names(sal$SE) <- step_name
    }
    sal <- as(sal, "SYSargsList")
    on.exit(return(sal))
}

####################
## runWF function ##
#####################
.renderMsg <- function(){
    warn_flag <- getOption("spr_render_msg")
    if(isTRUE(warn_flag)) return()
    cat(
        crayon::green$bold("Done with workflow running, now consider rendering logs & reports\n"),
        crayon::blue("To render logs, run:    "), "renderLogs(sal)\n",
        crayon::blue("From command-line:      "), 'Rscript -e "sal=systemPipeR::SPRproject(resume=TRUE);systemPipeR::renderLogs(sal)"\n',
        crayon::blue("To render reports, run: "), 'renderReport(sal)\n',
        crayon::blue("From command-line:      "), 'Rscript -e "sal=systemPipeR::SPRproject(resume=TRUE);systemPipeR::renderReport(sal)"\n',
        crayon::make_style("white")$bold("This message is displayed once per R session\n"),
        sep = ""
    )
    options(spr_render_msg = TRUE)
}

runWF <- function(sysargs, steps = NULL, targets = NULL,
                  force = FALSE, saveEnv = TRUE,
                  run_step = "ALL", ignore.dep = FALSE,
                  warning.stop = FALSE, error.stop = TRUE,
                  silent = FALSE, ...) {
    # Validations
    if (!inherits(sysargs, "SYSargsList")) stop("Argument 'sysargs' needs to be assigned an object of class 'SYSargsList'")
    if (length(sysargs) == 0) message("Workflow has no steps. Please add a step before trying to execute the workflow.")
    if (is.null(projectInfo(sysargs)$project)) {
        stop("'SYSargsList' instance was not initialized with the 'SPRproject' function
           and therefore there is missing the 'projectInfo' slot information.
           Please check 'SPRproject' help file.")
    }
    if (projectInfo(sysargs)$project != normalizePath(getwd())) {
        stop(
            "Your current working directory is different from the directory chosen for the Project Workflow.",
            "\n",
            "For accurate location of the files and running the Workflow, please set the working directory to ",
            "\n",
            "'setwd('", projectInfo(sysargs)$project, "')'"
        )
    }
    if (!dir.exists(projectInfo(sysargs)$logsDir)) stop("Project logsDir doesn't exist. Something went wrong...
        It is possible to restart the workflow saving the SYSargsList object with 'write_SYSargsList()' and restarting the project with 'SPRproject()'")
    sysproj <- projectInfo(sysargs)$logsDir
    ## targets selection
    targets_internal <- targets
    ## Check steps
    if (inherits(steps, "numeric")) {
        if (!all(steps %in% seq_along(sysargs))) {
            stop(
                "Please select the 'steps' number accordingly, options are: ", "\n",
                "        ", paste0(seq_along(sysargs), collapse = ", ")
            )
        }
    } else if (inherits(steps, "character")) {
        if (!all(steps %in% stepName(sysargs))) {
            stop(
                "Please select the 'steps' name accordingly, options are: ", "\n",
                "        ", paste0(stepName(sysargs), collapse = ", ")
            )
        }
        steps <- which(stepName(sysargs) %in% steps)
    }
    ## check dependency
    for (i in seq_along(dependency(sysargs))) {
        if (all(!is.na(dependency(sysargs)[[i]]))) {
            dep_names <- unlist(dependency(sysargs)[[i]])
            if (any(!dep_names %in% names(stepsWF(sysargs)))) {
                if (length(dep_names) == 1) {
                    stop(
                        "'sysargs' has a dependency on the following step:", "\n",
                        paste0("      ", paste0(dep_names)), "\n",
                        "However, this step is not available in the workflow."
                    )
                } else if (length(dep_names) > 1) {
                    stop(
                        "'sysargs' has dependencies on the following steps::", "\n",
                        paste0("      ", paste0(dep_names, collapse = " AND ")), "\n",
                        "However, these steps are not available in the workflow."
                    )
                }
            }
        }
    }
    ## Logs
    file_log <- file.path(sysproj, paste0("_logWF_", format(Sys.time(), "%b%d%Y_%H%M"), "_", paste(sample(0:9, 4), collapse = "")))
    sysargs[["projectInfo"]]$logsFile <- file_log
    ## Select steps for the loop based on the steps argument
    args2 <- sysargs
    if (is.null(steps)) steps <- 1:length(args2)
    ## Select steps for the loop based on the run_step
    if (run_step != "ALL") {
        run_step_wf <- sapply(sysargs$runInfo$runOption, function(x) x$run_step)[steps]
        run_step_names <- names(run_step_wf[run_step_wf %in% run_step])
        if (length(run_step_names) == 0) {
            stop("No step was selected. Please check 'steps', 'run_step' arguments.")
        } else {
            steps <- which(stepName(sysargs) %in% run_step_names)
        }
    }
    ## Selecting
    for (i in seq_along(stepsWF(args2))) {
        if (i %in% steps) {
            ## check single dependency
            if (!ignore.dep) {
                if (all(!is.na(dependency(args2)[[i]]))) {
                    dep_single <- sapply(dependency(args2)[[i]], function(x) args2$statusWF[[x]]$status.summary)
                    if ("Pending" %in% dep_single) {
                        message("Previous steps:", "\n", paste0(names(dep_single), collapse = " AND "), "\n", "have been not executed yet.")
                        break()
                    }
                }
            }
            ## Printing step name
            single.step <- names(stepsWF(args2)[i])
            cat(crayon::bgMagenta("Running Step: ", single.step), "\n")
            ## Printing single.step name at log files
            cat("# ", names(stepsWF(args2)[i]), "\n", file = file_log, fill = TRUE, append = TRUE)
            args.run <- stepsWF(args2)[[i]]
            run_location <- args2$runInfo$runOption[[i]]$run_session
            run_resorces <- args2$runInfo$runOption[[i]]$run_remote_resources
            envir <- args2$runInfo$env
            ## SYSargs2 STEP
            if (inherits(args.run, "SYSargs2")) {
                ## runC arguments
                dir <- args2$runInfo$runOption[[i]]$directory
                dir.name <- single.step
                ## check run_targets
                run_targets <- targets_internal
                ## run_targets input
                if (!is.null(run_targets)) {
                    if (inherits(run_targets, c("numeric", "integer"))) {
                        if (!all(run_targets %in% seq_along(cmdlist(args.run)))) {
                            stop(
                                "Please select the 'targets' number accordingly, options are: ", "\n",
                                "        ", paste0(seq_along(cmdlist(args.run)), collapse = ", ")
                            )
                        }
                    } else if (inherits(run_targets, "character")) {
                        if (!all(run_targets %in% SampleName(args.run))) {
                            stop(
                                "Please select the 'targets' name accordingly, options are: ", "\n",
                                "        ", paste0(SampleName(args.run), collapse = ", ")
                            )
                        }
                        run_targets <- which(SampleName(args.run) %in% run_targets)
                    }
                    run_targets <- run_targets
                } else {
                    run_targets <- seq_along(cmdlist(args.run))
                }
                ## assign to the envir
                assign(x = as.character(as.list(match.call())$sysargs), args2, envir = envir)
                if (run_location == "management") {
                    cat(crayon::bgMagenta("Running Session: Management"), "\n")
                    # RUN
                    args.run <- runCommandline(args.run,
                        dir = dir, dir.name = dir.name,
                        runid = paste0("_", dir.name),
                        force = force, input_targets = run_targets, ...
                    )
                } else if (run_location == "compute") {
                    cat(crayon::bgMagenta("Running Session: Compute Session"), "\n")
                    if (is.null(run_resorces)) stop("Resources are not available for this step.")
                    reg <- clusterRun(args.run[c(run_targets)],
                        FUN = runCommandline,
                        more.args = list(
                            args = args.run[c(run_targets)], dir = dir, dir.name = dir.name,
                            runid = paste0("_", dir.name),
                            force = force, ...
                        ),
                        conffile = run_resorces$conffile,
                        template = run_resorces$template,
                        Njobs = run_resorces$Njobs, runid = paste0("_", dir.name),
                        resourceList = run_resorces
                    )
                    batchtools::waitForJobs()
                    args.run <- .clusterRunResults(args.run, reg, nTargets = length(run_targets))
                    ## time
                    args.run@status$status.time$time_start <- as.POSIXct(args.run@status$status.time$time_start, origin = "1970-01-01")
                    args.run@status$status.time$time_end <- as.POSIXct(args.run@status$status.time$time_end, origin = "1970-01-01")
                    # cat(readLines(args.run$files$log),
                    #     file = file_log, sep = "\n",
                    #     append = TRUE
                    # )
                    ## moving file #reg$file.dir to new directory
                    output_path_args <- .getPath(args.run$output[[run_targets[1]]][[1]][1])
                    output_path_reg <- .getPath(reg$file.dir)
                    if (!identical(output_path_args, output_path_reg)) {
                        file.copy(reg$file.dir, file.path(output_path_args, "_logs"), recursive = TRUE)
                        message("\n", "Moving registry to '", file.path(output_path_args, "_logs"), "'.", "\n")
                        unlink(reg$file.dir, recursive = TRUE)
                        message("Jobs completed successfully!")
                        args.run@files[["log_cluster"]] <- file.path(output_path_args, "_logs", basename(reg$file.dir))
                        args.run@files[["log"]] <- file.path(output_path_args, "_logs", basename(reg$file.dir), basename(args.run@files$log))
                    }
                    ## Printing
                    cat(crayon::blue("---- Summary ----"), "\n")
                    print(args.run$status$status.completed)
                }
                run_targets <- NULL
                ## update object
                step.status.summary <- status(args.run)$status.summary
                statusWF(args2, single.step) <- args.run$status
                stepsWF(args2, single.step) <- args.run
                args2[["outfiles"]][[single.step]] <- .outList2DF(args.run)
                args2 <- .updateAfterRunC(args2, single.step)

                ## assign to the envir
                assign(x = as.character(as.list(match.call())$sysargs), args2, envir = envir)
                ## Stop workflow
                if (is.element("Warning", unlist(step.status.summary))) {
                    if (warning.stop == TRUE) {
                        on.exit(return(args2))
                        stop("Caught an warning, stop workflow!")
                    }
                } else if (is.element("Error", unlist(step.status.summary))) {
                    if (error.stop == TRUE) {
                        on.exit(return(args2))
                        stop("Caught an error, stop workflow!")
                    }
                }
                ## logs
                cat(readLines(args2[["stepsWF"]][[single.step]]$files$log),
                		file = file_log, sep = "\n",
                		append = TRUE
                )
                cat(status_color(step.status.summary)("Step Status: ", step.status.summary, "\n"))
                ## LineWise ##
            } else if (inherits(args.run, "LineWise")) {
                if (!is.null(targets_internal)) message("'targets' argument has been ignored since this step is 'LineWise' instance.")
                ## assign to the envir
                assign(x = as.character(as.list(match.call())$sysargs), args2, envir = envir)
                if (!dir.exists(file.path(sysproj, "Rsteps"))) {
                    dir.create(file.path(sysproj, "Rsteps"))
                }
                file_log_Rcode <- file.path(sysproj, "Rsteps", paste0("_logRstep_", single.step, "_", format(Sys.time(), "%b%d%Y_%H%M%S")))
                if (run_location == "management") {
                    cat(crayon::bgMagenta("Running Session: Management"), "\n")
                    # RUN
                    args.run <- runRcode(args.run,
                        step = single.step, file_log_Rcode = file_log_Rcode,
                        envir = envir, force = force
                    )
                } else if (run_location == "compute") {
                    cat(crayon::bgMagenta("Running Session: Compute Session"), "\n")
                    loaded_pkgs <- .packages()
                    tempImage <- file.path(sysproj, "Rsteps", "step_envir.RData")
                    save(list = viewEnvir(args2, silent = TRUE), file = tempImage, envir = envir)
                    list_return <- clusterRCode(args.run,
                        step = single.step, sysproj = sysproj,
                        file_log_Rcode = file_log_Rcode, force = force,
                        tempImage = tempImage, loaded_pkgs = loaded_pkgs,
                        resources = run_resorces
                    )
                    load(file.path(list_return$tempImage))
                    lapply(list_return$loaded_pkgs, require, character.only = TRUE)
                    args.run <- list_return$args
                    time_start <- args.run$status$status.time$time_start[!is.na(args.run$status$status.time$time_start)][1]
                    time_end_length <- length(args.run$status$status.time$time_end[!is.na(args.run$status$status.time$time_end)])
                    time_end <- args.run$status$status.time$time_end[!is.na(args.run$status$status.time$time_end)][time_end_length]
                    args.run@status$total.time <- list(
                        time_start = as.POSIXct(time_start, origin = "1970-01-01"),
                        time_end = as.POSIXct(time_end, origin = "1970-01-01")
                    )
                    unlink(tempImage)
                }
                ## time - double-check
                if (is.numeric(args.run@status$total.time$time_start)) {
                    args.run@status$total.time$time_start <- as.POSIXct(args.run@status$total.time$time_start, origin = "1970-01-01")
                }
                if (is.numeric(args.run@status$total.time$time_end)) {
                    args.run@status$total.time$time_end <- as.POSIXct(args.run@status$total.time$time_start, origin = "1970-01-01")
                }
                ## assign all the new object to the envir
                assign(
                    "args2", get(as.character(as.list(match.call())$sysargs), envir),
                    environment()
                )
                ## update object
                statusWF(args2, single.step) <- args.run$status
                stepsWF(args2, single.step) <- args.run
                cat(readLines(file_log_Rcode),
                    file = file_log, sep = "\n",
                    append = TRUE
                )
                assign(
                    x = as.character(as.list(match.call())$sysargs), args2,
                    envir = envir
                )
                ## Stop workflow
                if (is.element("Warning", unlist(args.run$status$status.summary))) {
                    if (warning.stop == TRUE) {
                        on.exit(return(args2))
                        stop("Caught an warning, stop workflow!")
                    }
                } else if (is.element("Error", unlist(args.run$status$status.summary))) {
                    if (error.stop == TRUE) {
                        on.exit(return(args2))
                        stop("Caught an error, stop workflow!")
                    }
                }
                cat(status_color(args.run$status$status.summary)("Step Status: ", args.run$status$status.summary, "\n"))
                run_targets <- NULL
            }
        } else {
            ## Printing step name
            single.step <- names(stepsWF(args2)[i])
            cat(status_color("Skipping")("Skipping Step: ", single.step), "\n")
            run_targets <- NULL
        }
    }
    ## Each time the workflow is executed, a new logs report should be generated (HTML format), via sal <- renderLogs(sal)
    if(!is.null(args2$projectInfo$Report)){
    	args2@projectInfo$Report_Logs <- NA
    }
    assign(
    	x = as.character(as.list(match.call())$sysargs), args2,
    	envir = envir
    )
    if (saveEnv == TRUE) {
        envPath <- file.path(sysproj, "sysargsEnv.rds")
        if (any(as.character(as.list(match.call())$sysargs) %in% ls(sysargs@runInfo$env, all.names = TRUE))) {
            rm(list = as.character(as.list(match.call())$sysargs), envir = sysargs@runInfo$env)
        }
        saveRDS(args2$runInfo$env, envPath)
        args2[["projectInfo"]][["envir"]] <- envPath
    }
    if (!silent) .renderMsg()
    args2 <- .check_write_SYSargsList(args2, TRUE)
    on.exit(return(args2))
}
## Usage:
## runWF(sal)

############################
## clusterRCode function ##
############################
clusterRCode <- function(args.run, step, sysproj, file_log_Rcode, force, tempImage, loaded_pkgs, resources) {
    checkPkg("batchtools", quietly = FALSE)
    if (!file.exists(tempImage)) stop("Something went wrong, temporary 'image' doesn't exist.", call. = FALSE)
    if (all(args.run$status$status.summary == "Success" && !force)) {
    	## Print at the log_file
    	cat(c(
    		paste0("Time: ", paste0(format(Sys.time(), "%b%d%Y_%H%Ms%S"))), "\n",
    		"## Code: ",
    		"```{r, eval=FALSE} ",
    		utils::capture.output(codeLine(args.run)),
    		"```", "\n",
    		"## Stdout: ",
    		"```{r, eval=FALSE}",
    		"The expected output file(s) already exist", "\n"
    	), file = file_log_Rcode, sep = "\n", fill = TRUE, append = TRUE)
        ## close R chunk
        cat("```", file = file_log_Rcode, sep = "\n", append = TRUE)
        args.run@files$log <- file_log_Rcode
        list_return <- list(args = args.run, loaded_pkgs = loaded_pkgs, tempImage = tempImage)
        return(list_return)
    }
    ## Function definition
    fct <- function(i, args.run, step, file_log_Rcode, force,
                    tempImage, loaded_pkgs, ...) {
        load(file.path(tempImage))
        ls_list1 <- ls()
        lapply(loaded_pkgs, require, character.only = TRUE)
        args <- runRcode(
            args = args.run, step = step, file_log_Rcode = file_log_Rcode,
            envir = environment(), force = force
        )
        loaded_pkgs <- .packages()
        save(list = ls(all.names = TRUE), file = tempImage, envir = environment())
        ls_list2 <- ls(envir = environment())
        return(list(args = args, loaded_pkgs = loaded_pkgs, tempImage = tempImage, ls1 = ls_list1, ls2 = ls_list2))
    }
    ##
    logdir1 <- paste0(sysproj, "/submitargs", 01, "_btdb_", paste(sample(0:9, 4), collapse = ""))
    conffile <- resources$conffile
    template <- resources$template
    reg <- batchtools::makeRegistry(file.dir = logdir1, conf.file = conffile, packages = "systemPipeR")
    ids <- batchtools::batchMap(fun = fct, 1, more.args = list(
        args.run = args.run, step = step, file_log_Rcode = file_log_Rcode,
        force = force, tempImage = tempImage,
        loaded_pkgs = loaded_pkgs
    ), reg = reg)
    chunk <- batchtools::chunk(ids$job.id, n.chunks = 1, shuffle = FALSE)
    ids$chunk <- chunk
    done <- batchtools::submitJobs(ids = ids, reg = reg, resources = resources)
    batchtools::waitForJobs()
    list_return <- batchtools::loadResult(reg = reg, id = ids)
    message("Jobs completed successfully!")
    return(list_return)
}

###########################
## status_color function ##
##########################
status_color <- function(x) {
    switch(x,
        "Pending" = crayon::bgBlue,
        "Warning" = crayon::bgYellow,
        "Error" = crayon::bgRed,
        "Success" = crayon::bgGreen,
        "Skipping" = crayon::bgCyan
    )
}

# cat(status_color("Pending")("test"))
# cat(status_color("Warning")("test"))
# cat(status_color("Error")("test"))
# cat(status_color("Success")("test"))
# cat(status_color("Skipping")("test"))

#######################
## runRcode function ##
#######################
runRcode <- function(args, step = stepName(args), file_log_Rcode = NULL, envir = globalenv(), force = FALSE) {
    ## Validation for 'args'
    if (!inherits(args, "LineWise")) stop("Argument 'args' needs to be assigned an object of class 'LineWise'")
    pb <- txtProgressBar(min = 0, max = length(args), style = 3)
    ## log_file
    if (is.null(file_log_Rcode)) {
        file_log_Rcode <- paste0("_logRcode_", format(Sys.time(), "%b%d%Y_%H%M"))
    }
    ## Print at the log_file
    cat(c(
        paste0("Time: ", paste0(format(Sys.time(), "%b%d%Y_%H%Ms%S"))), "\n",
        "## Code: ",
        "```{r, eval=FALSE} ",
        utils::capture.output(codeLine(args)),
        "```", "\n",
        "## Stdout: ",
        "```{r, eval=FALSE}"
    ), file = file_log_Rcode, sep = "\n", append = TRUE)
    ## Check status of step
    if (all(args$status$status.summary == "Success" && force == FALSE)) {
        args[["status"]]$status.time$time_start <- Sys.time()
        cat("The step status is 'Success' and it was skipped.", file = file_log_Rcode, fill = TRUE, append = TRUE, sep = "\n")
        args[["status"]]$status.time$time_end <- Sys.time()
    } else {
        ## Status and time register
        step_status <- list()
        time_status <- data.frame(Step = step, time_start = NA, time_end = NA)
        time_status$time_start <- Sys.time()
        ## Running the code
        stdout <- .tryRcode(args$codeLine, envir = envir)
        ## save stdout to file
        utils::capture.output(stdout$stdout, file = file_log_Rcode, append = TRUE)
        ## save error and warning messages
        if (!is.null(stdout$error)) {
            cat("## Error", file = file_log_Rcode, sep = "\n", append = TRUE)
            cat(stdout$error, file = file_log_Rcode, sep = "\n", append = TRUE)
            step_status[["status.summary"]] <- "Error"
        } else if (!is.null(stdout$warning)) {
            cat("## Warning", file = file_log_Rcode, sep = "\n", append = TRUE)
            cat(stdout$warning, file = file_log_Rcode, sep = "\n", append = TRUE)
            step_status[["status.summary"]] <- "Warning"
        } else if (all(is.null(stdout$error) && is.null(stdout$warning))) {
            step_status[["status.summary"]] <- "Success"
        }
        ## Saving the new status
        step_status[["status.completed"]] <- data.frame(Step = step, Status = step_status[[1]])
        time_status$time_end <- Sys.time()
        step_status[["status.time"]] <- time_status
        args[["status"]] <- step_status
        args[["status"]]$total.time <- list(time_start = args$status$status.time$time_start, time_end = args$status$status.time$time_end)
    }
    utils::setTxtProgressBar(pb, length(args))
    ## close R chunk
    cat("``` \n", file = file_log_Rcode, sep = "\n", append = TRUE)
    close(pb)
    args[["files"]] <- list(log = file_log_Rcode)
    return(args)
}

#############################
## output.as.df function ##
#############################
output.as.df <- function(x) {
    out_x <- output(x)
    out_x <- S4Vectors::DataFrame(matrix(unlist(out_x), nrow = length(out_x), byrow = TRUE))
    colnames(out_x) <- x$files$output_names
    return(out_x)
}

################################
## write_SYSargsList function ##
################################
write_SYSargsList <- function(sysargs, sys.file = ".SPRproject/SYSargsList.yml", silent = TRUE) {
    ## check logDir folder
    logDir <- .getPath(sys.file, warning = FALSE, normalizePath = FALSE)
    if (!file.exists(logDir)) stop("'logs.dir': No such file or directory. Check the file PATH.")
    if (!inherits(sysargs, "SYSargsList")) stop("sysargs needs to be object of class 'SYSargsList'.")
    args2 <- sysargslist(sysargs)
    args_comp <- sapply(args2, function(x) list(NULL))
    steps <- names(stepsWF(sysargs))
    ## special case for "runInfo" slot
    yaml_slots <- c("runInfo")
    for (i in yaml_slots) {
        args_comp[[i]] <- yaml::as.yaml(args2[[i]]$runOption)
    }
    ## Simple yaml slots
    yaml_slots <- c("projectInfo")
    for (i in yaml_slots) {
        args_comp[[i]] <- yaml::as.yaml(args2[[i]])
    }
    ## Yaml Slots + steps
    yaml_slots_S <- c("statusWF", "dependency", "targets_connection")
    for (i in yaml_slots_S) {
        steps_comp <- sapply(steps, function(x) list(NULL))
        for (j in steps) {
            steps_comp[j] <- yaml::as.yaml(args2[[i]][j])
        }
        args_comp[[i]] <- steps_comp
    }
    ## DataFrame Slots
    df_slots <- c("targetsWF", "outfiles")
    for (i in df_slots) {
        #  args_comp[[i]] <- yaml::as.yaml(as.data.frame(args2[[i]]$Mapping))
        steps_comp <- sapply(steps, function(x) list(NULL))
        for (j in steps) {
            steps_comp[j] <- yaml::as.yaml(data.frame(args2[[i]][[j]], check.names = FALSE))
        }
        args_comp[[i]] <- steps_comp
    }
    ## SYSargs2 and LineWise
    steps_comp <- sapply(steps, function(x) list(NULL))
    for (j in steps) {
        if (inherits(args2[["stepsWF"]][[j]], "SYSargs2")) {
            step_obj <- sysargs2(args2[["stepsWF"]][[j]])
            steps_comp[[j]] <- yaml::as.yaml(step_obj)
        } else if (inherits(args2[["stepsWF"]][[j]], "LineWise")) {
            step_obj <- linewise(args2[["stepsWF"]][[j]])
            step_obj$codeLine <- as.character(step_obj$codeLine)
            steps_comp[[j]] <- yaml::as.yaml(step_obj)
        }
    }
    args_comp[["stepsWF"]] <- steps_comp
    ## SE slot
    path <- file.path(.getPath(sys.file, full_path = TRUE, normalizePath = FALSE, warning = FALSE), "SE")
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    steps_comp <- sapply(steps, function(x) list(NULL))
    for (j in steps) {
        if (!is.null(args2[["SE"]][[j]])) {
            writeSE(args2[["SE"]][[j]], dir.path = path, dir.name = j, overwrite = TRUE, silent = TRUE)
            steps_comp[[j]] <- yaml::as.yaml(list(dir.path = path, dir.name = j))
        } else {
            steps_comp[j] <- yaml::as.yaml(args2[["SE"]][j])
        }
    }
    args_comp[["SE"]] <- steps_comp
    ## Save file
    yaml::write_yaml(args_comp, sys.file)
    if (silent != TRUE) cat("Creating file '", file.path(sys.file), "'", sep = "", "\n")
    return(sys.file)
}

# ## Usage:
# write_SYSargsList(sal, sys.file, silent=FALSE)

################################
## read_SYSargsList function ##
################################
read_SYSargsList <- function(sys.file) {
    args_comp_yml <- yaml::read_yaml(sys.file)
    args_comp <- sapply(args_comp_yml, function(x) list(NULL))
    steps <- names(args_comp_yml$stepsWF)
    ## Simple yaml slots
    yaml_slots <- c("projectInfo")
    for (i in yaml_slots) {
        args_comp[[i]] <- yaml::yaml.load(args_comp_yml[i])
        if (!is.null(args_comp[[i]]$rmd_lines)) {
            args_comp[[i]]$rmd_lines <- data.frame(args_comp[[i]]$rmd_lines, check.names = FALSE)
        }
    }
    ## runInfo yaml slots
    yaml_slots <- c("runInfo")
    for (i in yaml_slots) {
        args_comp[[i]] <- list(runOption = yaml::yaml.load(args_comp_yml[i]))
    }
    ## Yaml Slots + steps
    yaml_slots_S <- c("dependency", "targets_connection")
    for (i in yaml_slots_S) {
        steps_comp <- sapply(steps, function(x) list(NULL))
        for (j in steps) {
            steps_comp[j] <- yaml::yaml.load(args_comp_yml[[i]][j])
        }
        args_comp[[i]] <- steps_comp
    }
    ## SYSargs2 and LineWise
    if (length(args_comp_yml$stepsWF) >= 1) {
        steps_comp <- sapply(steps, function(x) list(NULL))
        for (j in steps) {
            if ("codeLine" %in% names(yaml::yaml.load(args_comp_yml[["stepsWF"]][[j]]))) {
                args <- yaml::yaml.load(args_comp_yml[["stepsWF"]][[j]])
                if (length(args$codeLine) >= 1) {
                    args$codeLine <- parse(text = args$codeLine)
                }
                if (length(args$codeChunkStart) == 0) args$codeChunkStart <- integer()
                if (length(args$files) == 0) args$files$rmdPath <- character()
                if (length(args$dependency) == 0) args$dependency <- character()
                args <- as(args, "LineWise")
                steps_comp[[j]] <- args
            } else {
                args <- yaml::yaml.load(args_comp_yml[["stepsWF"]][[j]])
                args[["status"]][[2]] <- data.frame(args[["status"]][[2]], check.names = FALSE)
                args[["status"]][[3]] <- data.frame(args[["status"]][[3]])
                if (length(names(args$targets)) != 0) {
                    df_rownames <- names(args$targets)
                } else {
                    df_rownames <- args[["status"]][[3]]$Targets
                }
                rownames(args[["status"]][[2]]) <- df_rownames
                rownames(args[["status"]][[3]]) <- df_rownames
                steps_comp[[j]] <- as(args, "SYSargs2")
            }
            args_comp[["stepsWF"]] <- steps_comp
        }
    } else if (length(args_comp_yml$stepsWF) >= 0) {
        args_comp[["stepsWF"]] <- list()
    }
    ## targetsWF Slots
    df_slots <- c("targetsWF", "outfiles")
    for (i in df_slots) {
        steps_comp <- sapply(steps, function(x) list(NULL))
        for (j in steps) {
            steps_comp[[j]] <- S4Vectors::DataFrame(yaml::yaml.load(args_comp_yml[[i]][[j]]), check.names = FALSE)
        }
        args_comp[[i]] <- steps_comp
    }
    ## status
    yaml_slots_Status <- c("statusWF")
    for (i in yaml_slots_Status) {
        steps_comp <- sapply(steps, function(x) list(NULL))
        for (j in steps) {
            steps_comp[j] <- yaml::yaml.load(args_comp_yml[[i]][j])
            steps_comp[[j]]$status.completed <- data.frame(steps_comp[[j]]$status.completed, check.names = FALSE)
            steps_comp[[j]]$status.time <- data.frame(steps_comp[[j]]$status.time)
            if (length(steps_comp[[j]]$status.time) != 0) {
                steps_comp[[j]]$status.time$time_start <- .POSIXct(steps_comp[[j]]$status.time$time_start)
                steps_comp[[j]]$status.time$time_end <- .POSIXct(steps_comp[[j]]$status.time$time_end)
            }
            if (length(steps_comp[[j]]$total.time) != 0) {
                steps_comp[[j]]$total.time$time_start <- .POSIXct(steps_comp[[j]]$total.time$time_start)
                steps_comp[[j]]$total.time$time_end <- .POSIXct(steps_comp[[j]]$total.time$time_end)
            }
        }
        args_comp[[i]] <- steps_comp
    }
    ## rownames
    for (j in steps) {
        if (inherits(args_comp$stepsWF[[j]], "SYSargs2")) {
            rownames(args_comp$targetsWF[[j]]) <- args_comp$targetsWF[[j]][[args_comp$stepsWF[[j]]$files$id]]
            rownames(args_comp$outfiles[[j]]) <- args_comp$targetsWF[[j]][[args_comp$stepsWF[[j]]$files$id]]
            rownames(args_comp$statusWF[[j]]$status.completed) <- args_comp$targetsWF[[j]][[args_comp$stepsWF[[j]]$files$id]]
            rownames(args_comp$statusWF[[j]]$status.time) <- args_comp$targetsWF[[j]][[args_comp$stepsWF[[j]]$files$id]]
        }
    }
    ## SE slot
    steps_comp <- sapply(steps, function(x) list(NULL))
    for (j in steps) {
        steps_comp[[j]] <- yaml::yaml.load(args_comp_yml[["SE"]][[j]])
        if (!is.null(steps_comp[[j]][[1]])) {
            dir.path <- steps_comp[[j]][[1]]
            dir.name <- steps_comp[[j]][[2]]
            SE <- readSE(dir.path = dir.path, dir.name = dir.name)
            steps_comp[[j]] <- list(SE)
        } else {
            steps_comp[[j]] <- steps_comp[[j]]
        }
    }
    args_comp[["SE"]] <- sapply(steps_comp, function(x) x[[1]])
    return(as(args_comp, "SYSargsList"))
}

# ## Usage:
# sys.file=".SPRproject/SYSargsList.yml"
# sal3 <- read_SYSargsList(sys.file)

################################
## writeSE function ##
################################
writeSE <- function(SE, dir.path, dir.name, overwrite = FALSE, silent = FALSE) {
    # Validations
    if (!inherits(SE, "SummarizedExperiment")) stop("Argument 'SE' needs to be assigned an object of class 'SummarizedExperiment'")
    if (!dir.exists(dir.path)) stop("'dir.path' doesn't exist.")
    if (all(dir.exists(file.path(dir.path, dir.name)) & overwrite == FALSE)) stop("'dir.name' directory already exist. Please delete existing directory: ", dir.name, " or set 'overwrite=TRUE'")
    if (!dir.exists(file.path(dir.path, dir.name))) {
        dir.create(file.path(dir.path, dir.name))
    }
    path <- file.path(dir.path, dir.name)
    ## Counts
    if (length(SummarizedExperiment::assays(SE)) > 0) {
        for (i in length(SummarizedExperiment::assays(SE))) {
            write.table(SummarizedExperiment::assays(SE)[[i]], file.path(path, paste0("counts_", i, ".csv")),
                quote = FALSE, row.names = TRUE,
                col.names = TRUE, sep = "\t"
            )
        }
    }
    ## Metadata
    yaml::write_yaml(S4Vectors::metadata(SE), file.path(path, paste0("metadata.yml")))
    ## colData
    write.table(SummarizedExperiment::colData(SE), file.path(path, paste0("colData.csv")),
        quote = FALSE, row.names = TRUE,
        col.names = NA, sep = "\t"
    )
    ## RowRanges
    if (!is.null(SummarizedExperiment::rowRanges(SE))) {
        write.table(as.data.frame(SummarizedExperiment::rowRanges(SE)),
            file = file.path(path, paste0("rowRanges.csv")),
            sep = "\t", quote = FALSE, row.names = FALSE
        )
    }
    ## Final message
    if (silent != TRUE) cat("\t", "Written content of 'SE' to directory:", path, "\n")
}

# writeSE(rse, dir.path = getwd(), dir.name = "seobj")
#
# dir.path <- getwd()
# dir.name <- "seobj"

################################
## readSE function ##
################################
readSE <- function(dir.path, dir.name) {
    path <- file.path(dir.path, dir.name)
    if (!dir.exists(path)) stop("'dir.path' doesn't exist.")
    ## Counts
    files_counts <- list.files(path, pattern = "counts")
    if (length(files_counts) > 0) {
        counts_ls <- S4Vectors::SimpleList()
        for (i in files_counts) {
            counts_ls1 <- S4Vectors::SimpleList(as.matrix(read.table(file.path(path, i), check.names = FALSE, header = TRUE)))
            counts_ls <- append(counts_ls, counts_ls1)
        }
    } else {
        counts_ls <- S4Vectors::SimpleList()
    }
    ## Metadata
    metadata <- yaml::read_yaml(file.path(path, paste0("metadata.yml")))
    ## colData
    colData <- tryCatch(read.table(file.path(path, paste0("colData.csv")), check.names = FALSE, sep = "\t", header = TRUE),
        error = function(e) NULL
    )
    if (is.null(colData)) colData <- data.frame()
    ## rowRanges
    files_counts <- list.files(path, pattern = "rowRanges")
    if (length(files_counts) > 0) {
        rowRanges_df <- read.table(file.path(path, paste0("rowRanges.csv")), check.names = FALSE, header = TRUE)
        rowRanges <- makeGRangesFromDataFrame(rowRanges_df, keep.extra.columns = TRUE)
    } else {
        rowRanges <- GRangesList()
    }
    if (length(rowRanges) == 0) {
        SE <- SummarizedExperiment::SummarizedExperiment(
            assays = counts_ls,
            # rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    } else {
        SE <- SummarizedExperiment::SummarizedExperiment(
            assays = counts_ls,
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    }
    return(SE)
}

# nrows <- 200; ncols <- 6
# counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#                      IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#                      strand=sample(c("+", "-"), 200, TRUE),
#                      feature_id=sprintf("ID%03d", 1:200))
# colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
#                      row.names=LETTERS[1:6])
# rse <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(counts=counts),
#                             rowRanges=rowRanges, colData=colData)
# rse
# writeSE(rse, dir.path = getwd(), dir.name = "seobj")
# SE <- readSE(dir.path = getwd(), dir.name = "seobj")
# library(diffobj)
# diffPrint(target=SE, current=rse)

# writeSE(sal$SE$gzip, dir.path = getwd(), dir.name = "seobj")
# SE <- readSE(dir.path = getwd(), dir.name = "seobj")
# library(diffobj)
# diffPrint(target=SE, current=rse)


###########################
## writeTargets function ##
###########################
writeTargets <- function(sysargs, step, file = "default", silent = FALSE, overwrite = FALSE) {
    if (all(!inherits(sysargs, "SYSargsList"))) stop("Argument 'sysargs' needs to be assigned an object of class 'SYSargsList")
    if (!step %in% stepName(sysargs)) stop("It was not possible to find the 'step' name in the 'sysargs' object")
    ## Workflow and Step Name
    if (file == "default") {
        file <- paste("targets_", step, ".txt", sep = "")
    } else {
        file <- file
    }
    ## SYSargsList class
    if (file.exists(file) & overwrite == FALSE) stop("I am not allowed to overwrite files; please delete existing file: ", file, " or set 'overwrite=TRUE'")
    targets <- targetsWF(sysargs)[[step]]
    targetslines <- c(paste(colnames(targets), collapse = "\t"), apply(targets, 1, paste, collapse = "\t"))
    headerlines <- sysargs$stepsWF[[step]][["targetsheader"]]
    writeLines(c(headerlines$targetsheader, targetslines), file)
    if (silent != TRUE) cat("\t", "Written content of 'targetsout(x)' to file:", file, "\n")
}

########################
## ConfigWF function ##
########################
configWF <- function(x, input_steps = "ALL", exclude_steps = NULL, silent = FALSE, ...) {
    ## Validations
    if (!inherits(x, "SYSargsList")) stop("Argument 'x' needs to be assigned an object of class 'SYSargsList'")
    utils::capture.output(steps_all <- subsetRmd(Rmd = x$sysconfig$script$path), file = ".SYSproject/.NULL")
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
    utils::capture.output(steps_all <- subsetRmd(
        Rmd = x$sysconfig$script$path, input_steps = input_steps,
        exclude_steps = exclude_steps, save_Rmd = save_rmd, Rmd_outfile = Rmd_outfile
    ), file = ".NULL")
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

###########################
## renderReport function ##
###########################
## type: c("pdf_document", "html_document")
renderReport <- function(sysargs,
                         fileName = "SPR_Report",
                         rmd_title = "SPR workflow Template - Report",
                         rmd_author = "Author",
                         rmd_date = "Last update: `r format(Sys.time(), '%d %B, %Y')`",
                         type = c("html_document"),
                         desc = "This is a workflow template.",
                         quiet = FALSE,
                         open_file = TRUE) {
    stopifnot(is.logical(quiet) && length(quiet) == 1)
    stopifnot(is.logical(open_file) && length(open_file) == 1)
    stopifnot(inherits(sysargs, "SYSargsList"))

    if (is.null(sysargs$projectInfo$rmd_file)) return(
        sal2rmd(sysargs,
                out_path = out_path_rmd,
                rmd_title = rmd_title,
                rmd_author = rmd_author,
                rmd_date = rmd_date,
                rmd_output = type,
                desc = desc,
                verbose = !quiet
        )
    )

    ext <- sub("_.*", "", type)
    out_path <- file.path(paste0(fileName, ".", ext))
    out_path_rmd <- file.path(paste0(fileName, ".Rmd"))

    file_path <- sysargs$projectInfo$rmd_file

    lines <- readLines(file_path) ## original file
    ## get info from sysargs
    step_names <- names(sysargs$stepsWF)
    deps <- sysargs$dependency
    t_connects <- sysargs$targets_connection
    opts <- sysargs$runInfo$runOption
    ## get the first line
    chunk_start <- lines %>% stringr::str_which("^```\\{.*\\}.*")
    df <- parseRmd(file_path, verbose = FALSE)[, 1:3]
    # check steps
    df_in_sal <- df$step_name %in% step_names
    if(!all(df_in_sal)) {
        warning("Some steps exist in template but not in SYSargsList. ",
                "Step code and text of these steps will be removed in report\n",
                "Consider using `sal <- importWF(sal, 'xx.Rmd', update = TRUE)` ",
                "to sync the template and SYSargsList before rendering report\n",
                paste(df$step_name[!df_in_sal], collapse = " "), immediate. = TRUE)
    }
    # find text in between
    text_start <- c(NA, df$end + 1)
    df$text_start <- text_start[-length(text_start)]
    df$text_end <- df$start - 1
    df$text_end[1] <- NA
    # get starter lines
    rmd_text <- lines[seq(1, df$start[1] -1)]
    for (i in seq_along(step_names)) {
        appStep <- switch(class(sysargs$stepsWF[[i]]),
            "LineWise" = .sal2rmd_rstep(
                sysargs,
                con = NULL, i = i, step_name = step_names[i],
                dep = deps[[i]], req = opts[[i]]$run_step,
                prepro = opts[[i]]$prepro_lines,
                session = opts[[i]]$run_session, return = "object"
            ),
            "SYSargs2" = .sal2rmd_sysstep(
                sysargs,
                con = NULL, i = i,
                step_name = step_names[i],
                dep = deps[[i]],
                dir = opts[[i]]$directory,
                req = opts[[i]]$run_step,
                session = opts[[i]]$run_session,
                prepro = opts[[i]]$prepro_lines,
                t_con = t_connects[[i]], return = "object"
            )
        )
        # find sal step matched steps from df lines
        df_step <- df[df$step_name == step_names[i],]
        if(nrow(df_step) == 0) {
            message(crayon::yellow$bold(
                "Step name ", step_names[i], " exist in SYSargsList but not in template.\n ",
                "It will be added as a new step in report, consider add it to your template first, ",
                "and then use `sal <- importWF(sal, 'xx.Rmd', update = TRUE)`",
                "to update the template first before rendering report."
            ))
        } else if (nrow(df_step) > 1){
            stop(crayon::red$bold("Step name ", step_names[i], " has more than one match in template, something wrong?"))
        } else {
            # add text
            if(!is.na(df_step$text_start[1])) rmd_text <- c(rmd_text, lines[seq(df_step$text_start[1], df_step$text_end[1])])
        }
        # add code
        rmd_text <- c(rmd_text, if(nrow(df_step) == 0) appStep else appStep[-1])
    }
    # add tail
    rmd_text <- c(rmd_text, lines[seq(df$end[nrow(df)] + 1, length(lines))])
    writeLines(rmd_text, out_path_rmd)

    rmarkdown::render(input = out_path_rmd, c(paste0("BiocStyle::", type)), quiet = quiet, envir = new.env())
    detach("package:BiocStyle", unload = TRUE)
    if (!quiet) message(crayon::green$bold("Success! Report created at", out_path))
    sysargs <- as(sysargs, "list")
    sysargs$projectInfo[["Report"]] <- normalizePath(file.path(out_path))
    if (open_file) try(utils::browseURL(file.path(out_path)), TRUE)
    sysargs <- as(sysargs, "SYSargsList")
    write_SYSargsList(sysargs, silent = TRUE)
    return(sysargs)
}

## Usage:
# sal <- SPRproject(overwrite = TRUE)
# file_path <- system.file("extdata", "spr_simple_wf.Rmd", package = "systemPipeR")
# sal <- importWF(sal, file_path = file_path, verbose = FALSE)
# targetspath <- system.file("extdata/cwl/example/targets_example.txt", package = "systemPipeR")
# appendStep(sal) <- SYSargsList(step_name = "echo",
#                               targets = targetspath, dir = TRUE,
#                               wf_file = "example/workflow_example.cwl", input_file = "example/example.yml",
#                               dir_path = system.file("extdata/cwl", package = "systemPipeR"),
#                               inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
# sal <- renderReport(sal)

###########################
## renderLogs function ##
###########################
renderLogs <- function(sysargs,
                       type = c("html_document", "pdf_document"),
                       fileName = "default",
                       quiet = FALSE,
                       open_file = TRUE) {
    if (!inherits(sysargs, "SYSargsList")) stop("`sysargs` must be a 'SYSargsList' object.")
    type <- match.arg(type, c("html_document", "pdf_document"))
    stopifnot(is.character(fileName) && length(fileName) == 1)
    stopifnot(is.logical(quiet) && length(quiet) == 1)
    stopifnot(is.logical(open_file) && length(open_file) == 1)
    wd <- getwd()
    if (wd != projectInfo(sysargs)$project) {
        setwd(projectInfo(sysargs)$project)
        on.exit({
            try(setwd(wd), TRUE)
        })
    }
    dir_log <- projectInfo(sysargs)$logsDir
    if (!file.exists(dir_log) == TRUE) stop("Provide valid 'SYSargsList' object. Check the initialization of the project.")
    if (is.null(projectInfo(sysargs)$logsFile)) {
        log <- ""
    } else {
        file <- projectInfo(sysargs)$logsFile
        log <- readLines(file)
    }
    if (fileName == "default") {
        fileName <- file.path(projectInfo(sysargs)$project, paste0("logs_", format(Sys.time(), "%b%d%Y_%H%M"), ".Rmd"))
    } else {
        fileName <- file.path(paste0(fileName, ".Rmd"))
    }
    if (type == "html_document") plot_path <- normalizePath(.prepareRmdPlot(sysargs, dir_log))
    writeLines(c(
        "---",
        "title: SPR Workflow Technical Report",
        paste0("date: 'Last update: ", format(Sys.time(), "%d %B, %Y"), "'"),
        "output:",
        paste0("  BiocStyle::", type, ":"),
        if (type == "html_document") paste0("    includes:\n      before_body: ", plot_path),
        "    toc: true",
        "    toc_float: true",
        "    code_folding: show",
        "package: systemPipeR",
        "fontsize: 14pt",
        "---",
        "",
        log
    ),
    con = fileName
    )
    # rmarkdown::render(input = fileName, c(paste0(type)), quiet = TRUE, envir = new.env())
    rmarkdown::render(input = fileName, c(paste0("BiocStyle::", type)), quiet = TRUE, envir = new.env())
    detach("package:BiocStyle", unload = TRUE)
    file_path <- .getPath(fileName)
    file_out <- .getFileName(fileName)
    ext <- if (type == "html_document") "html" else "pdf"
    sysargs <- as(sysargs, "list")
    sysargs$projectInfo[["Report_Logs"]] <- file.path(file_path, paste(file_out, ext, sep = "."))
    if (!quiet) cat("Written content of 'Report' to file:", "\n", paste(file_out, ext, sep = "."), "\n")
    if (open_file) try(utils::browseURL(file.path(file_path, paste(file_out, ext, sep = "."))), TRUE)
    return(as(sysargs, "SYSargsList"))
}

## Usage
# sal <- runWF(sal)
# sal <- renderLogs(sal)

########################
## subsetRmd function ##
########################
subsetRmd <- function(Rmd, input_steps = NULL, exclude_steps = NULL, Rmd_outfile = NULL, save_Rmd = TRUE) {
    . <- NULL
    # function start, check inputs
    if (!file.exists(Rmd) == TRUE) stop("Provide valid 'Rmd' file.")
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
        if (chunk_start[i + 1] <= chunk_end[i]) stop("A code chunk does not end: chunk line ", chunk_start[i + 1])
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
    # if (!assertthat::not_empty(input_steps)) {
    if (!is.null(input_steps)) {
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
    cat("File write to ", file.path(Rmd_outfile), "\n")
    return(rmd_df)
}

# # Usage:
# Rmd <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
# newRmd <- subsetRmd(Rmd=Rmd, input_steps="1:2.1, 3.2:4, 4:6", exclude_steps="3.1", Rmd_outfile="test_out.Rmd", save_Rmd=TRUE)

###########################
## config.param function ##
###########################
config.param <- function(input_file = NULL, param, file = "default", silent = FALSE) {
    ## In the case of 'input_file' == character (file)
    if (inherits(input_file, "character")) {
        if (!file.exists(input_file)) {
            stop("Provide valid 'input_file' file. Check the file PATH.")
        }
        input <- yaml::read_yaml(file.path(input_file), eval.expr = TRUE)
        input <- out_obj <- .replace(input = input, param = param)
        path_file <- normalizePath(input_file)
        out_msg <- c("input_file")
    } else if (inherits(input_file, "list")) {
        if (is.null(names(param))) {
            stop("for each element of the 'param' list need to assign a name.")
        }
        input <- out_obj <- .replace(input = input_file, param = param)
        path_file <- normalizePath(file)
        out_msg <- c("input_file")
    } else if (inherits(input_file, "SYSargs2")) {
        input <- .replace(input = yamlinput(input_file), param = param)
        dir_path <- .getPath(files(input_file)[["yml"]])
        if (is.na(files(input_file)[["targets"]])) {
            targets <- NULL
        } else {
            targets <- targets(input_file)
        }
        args1 <- loadWorkflow(
            targets = targets, wf_file = basename(files(input_file)[["cwl"]]),
            input_file = basename(files(input_file)[["yml"]]), dir_path = dir_path
        )
        args1 <- as(args1, "list")
        args1$yamlinput <- input
        if ("ModulesToLoad" %in% names(param)) {
            for (i in seq_along(param$ModulesToLoad)) {
                args1$modules[names(param$ModulesToLoad[i])] <- param$ModulesToLoad[[i]]
            }
        }
        args1 <- out_obj <- as(args1, "SYSargs2")
        out_msg <- c("yamlinput(args1)")
        path_file <- files(input_file)[["yml"]]
    } else if (inherits(input_file, "SYSargsList")) {
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
    ## Write YML file
    yaml::write_yaml(x = input, file = path)
    if (silent != TRUE) cat("\t", "All the new param + ", out_msg, "were written to:", "\n", path, "\n")
    if (inherits(input_file, "SYSargs2")) {
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
            status = data.frame(),
            internal_outfiles = list()
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
            targets_connection = list(),
            projectInfo = list(),
            runInfo = list()
        )
        return(as(SYS.empty, "SYSargsList"))
    } else if (!class %in% c("SYSargs2", "SYSargsList")) {
        stop("Argument 'args' needs to be assigned an character 'SYSargs2' OR 'SYSargsList'")
    }
}
## Usage:
# list <- SYScreate("SYSargsList")
# list <- SYScreate("SYSargs2")

#########################################################################################
## Function to check if the command line / Software is installed and set in your PATH ##
#########################################################################################
tryCMD <- function(command, silent = FALSE) {
    if (command == "fastqc") command <- "fastqc --version"
    if (command == "gunzip") command <- "gunzip -h"
    if (command == "gzip") command <- "gzip -h"
    tryCatch(
        {
            system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
            if (!silent) print("All set, proceed!")
            if (silent) (return("proceed"))
        },
        warning = function(w) {
            if (silent) invisible(return("error"))
            if (!silent) {
                cat(
                    "ERROR: ", "\n", command, ": COMMAND NOT FOUND. ", "\n",
                    "Please make sure to configure your PATH environment variable according to the software in use.",
                    "\n"
                )
            }
        }
    )
}

## Usage:
# tryCMD(command="R")
# tryCMD(command="blastp")
# tryCMD(command="fastqc")
# tryCMD(command="gzip")

tryCL <- tryCMD
########################################################
## Function to check if the Path (file or dir) exists ##
########################################################
tryPath <- function(path) {
    tryCatch(normalizePath(path),
        warning = function(w) message(path, ": ", "No such file or directory"),
        error = function(e) message(path, ": ", "Please provide a valid Path")
    )
}
## Usage:
# tryPath(path="./")


########################################################################
## Function to list/check if all cmd tools are in path for a workflow ##
########################################################################
listCmdTools <- function(
    sal,
    check_path = FALSE,
    check_module = FALSE
){
    if(!inherits(sal, "SYSargsList")) stop("ERROR: Please provide a `SYSargsList` object")
    # list of cmd steps
    all_steps <- stepsWF(sal) %>% lapply(class) %>% unlist()
    cmd_steps <- which(all_steps == "SYSargs2")
    if(length(cmd_steps) == 0) return(cat(crayon::make_style("orange")$bold("There is no commandline (SYSargs) step in this workflow, skip.\n")))
    base_cmd <- lapply(cmd_steps, function(step) {
        # step_module <- stepsWF(sal)[[step]]@modules
        # if(length(step_module) > 0)
        step_clt <- stepsWF(sal)[[step]]@clt
       lapply(step_clt, function(clt) {
            clt$baseCommand[1]
        })
    })

    base_cmd <- data.frame(
        step_name = rep(names(base_cmd), (lapply(base_cmd, length) %>% unlist(use.names = FALSE))),
        tool = base_cmd %>% unlist(use.names = FALSE),
        in_path = NA
    )
    if(!check_path) {
        cat(crayon::blue$bold("Following tools are used in steps in this workflow:\n"))
        print(base_cmd)
    }

    if(!check_path) return(invisible(base_cmd))
    cat(crayon::blue$bold("Check if they are in path:\n"))
    unique_tools <- unique(base_cmd$tool)
    path_res <- lapply(unique_tools, function(i) {
        cat(crayon::blue$bold(paste0("Checking path for ", i, "\n")))
        try_res <- tryCMD(i, silent = TRUE)
        if(try_res == "proceed") {
            cat(crayon::blue$green$bold("PASS\n"))
            return(TRUE)
        }
        cat(crayon::blue$red$bold("ERROR\n"))
        FALSE
    }) %>% unlist()

    base_cmd$in_path <- path_res[match(base_cmd$tool, unique_tools)]
    print(base_cmd)

    if (all(path_res) && check_module) return(cat(crayon::green$bold(
        "All required tools in PATH, skip module check. If you want to check modules use `listCmdModules`"
    )))

    if(!check_module) return({
        if(!all(path_res)) {
            cat(crayon::blue$bold(
                c("Not all tools required are in PATH.",
                  "If you have a modular system, turn `load_module = TRUE`",
                  "It will call `listCmdModules` to check the required tools are in your modular system.\n")))
        }
        invisible(base_cmd)
    })

    cat(crayon::blue$bold("Now check if modular system:\n"))
    listCmdModules(sal, check_module)
    invisible(base_cmd)
}

listCmdModules <- function(
    sal,
    check_module = FALSE
) {
    if(!inherits(sal, "SYSargsList")) stop("ERROR: Please provide a `SYSargsList` object")
    # list of cmd steps
    all_steps <- stepsWF(sal) %>% lapply(class) %>% unlist()
    cmd_steps <- which(all_steps == "SYSargs2")
    if(length(cmd_steps) == 0) return(cat(crayon::make_style("orange")$bold("There is no commandline (SYSargs) step in this workflow, skip.\n")))
    base_mod <- lapply(cmd_steps, function(step) {
       stepsWF(sal)[[step]]@modules
    }) %>% unlist()
    if(is.null(base_mod)) return(cat(crayon::make_style("orange")$bold("No module is listed, check your CWL yaml configuration files, skip.\n")))
    base_mod <- data.frame(
        step_name = stringr::str_remove(names(base_mod), "\\.module[0-9]{1,}$"),
        module = unname(base_mod),
        avail = NA
    )

    if(!check_module) {
        cat(crayon::make_style("orange")$bold("Here are modules used in this workflow:\n"))
        print(base_mod)
    }

    if(!check_module) return(invisible(base_mod))
    module_avail <- is.modules.avail()
    if(is.null(module_avail)) return({
        cat(crayon::make_style("orange")$bold("Modular system is not installed, skip\n"))
        invisible(base_mod)
    })

    all_mods <- moduleAvail()$available_modules %>% unlist(use.names = FALSE)
    base_mod$avail <-  base_mod$module %in% all_mods
    cat(crayon::blue$bold("Module availability:\n"))
    print(base_mod)
}


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

###################################
## Update CWL description files  ##
###################################
cwlFilesUpdate <- function(destdir, force = FALSE, verbose = TRUE) {
    tempDir <- tempdir()
    download.file(
        url = "https://raw.githubusercontent.com/systemPipeR/cwl_collection/master/cwl/repo_version.txt",
        file.path(tempDir, "repo_version_new.txt"), quiet = TRUE
    )
    ## option for old versions
    if (force) {
        download.file(
            url = "https://github.com/systemPipeR/cwl_collection/archive/refs/heads/master.zip",
            file.path(tempDir, "cwl_collection-master.zip")
        )
        unzip(file.path(tempDir, "cwl_collection-master.zip"), exdir = tempDir)
        file.copy(file.path(tempDir, "cwl_collection-master", "cwl"), to = file.path(destdir), overwrite = TRUE, recursive = TRUE)
        file.copy(file.path(tempDir, "cwl_collection-master", "docopt.R"), to = file.path(destdir), overwrite = TRUE, recursive = TRUE)
        if (verbose) {
            cat(crayon::magenta(file.path(destdir, "cwl"), "folder was updated successfully!"))
        }
    } else {
        ## Get version
        new <- readLines(file.path(tempDir, "repo_version_new.txt"))
        if (!file.exists(file.path(destdir, "cwl", "repo_version.txt"))) {
            if (verbose) {
                cat(crayon::magenta(
                    "We expect a file called:", file.path(destdir, "cwl", "repo_version.txt"),
                    "\n", "OR please use the argument `force = TRUE`."
                ))
            }
        } else {
            current <- readLines(file.path(destdir, "cwl", "repo_version.txt"))
            ## Download CWL repo
            if (as.numeric(sub(".*\\.", "", current)) < as.numeric(sub(".*\\.", "", new))) {
                download.file(
                    url = "https://github.com/systemPipeR/cwl_collection/archive/refs/heads/master.zip",
                    file.path(tempDir, "cwl_collection-master.zip"), quiet = TRUE
                )
                unzip(file.path(tempDir, "cwl_collection-master.zip"), exdir = tempDir)
                file.copy(file.path(tempDir, "cwl_collection-master", "cwl"),
                    to = file.path(destdir), overwrite = TRUE, recursive = TRUE
                )
                file.copy(file.path(tempDir, "cwl_collection-master", "docopt.R"),
                    to = file.path(destdir), overwrite = TRUE, recursive = TRUE
                )
                if (verbose) {
                    cat(crayon::magenta(
                        file.path(destdir, "cwl"), "folder was updated successfully!", "\n",
                        new
                    ))
                }
            } else {
                if (verbose) {
                    cat(crayon::magenta("Already up to date.", "\n", new))
                }
            }
        }
    }
}
## Usage
# destdir <- "param/"
# options(timeout = max(3000, getOption("timeout")))
# cwlFilesUpdate(destdir)

#################################
## Unexported helper functions ##
#################################

##################################
## Return the path of the file  ##
##################################
## [x] A character vector or an object containing file PATH.
.getPath <- function(x, normalizePath = TRUE, full_path = TRUE, warning = TRUE) {
    if (warning) {
        if (!any(file.exists(x))) warning("No such file or directory. Check the file PATH.")
    }
    if (normalizePath) {
        x <- normalizePath(x)
    }
    if (full_path) {
        for (i in seq_along(x)) {
            path_un <- unlist(strsplit(x[i], "/|\\\\"))
            path <- path_un[path_un != basename(x[i])]
            x[i] <- paste0(path, collapse = "/")
        }
    } else {
        for (i in seq_along(x)) {
            path_un <- unlist(strsplit(x[i], "/|\\\\"))
            path <- path_un[path_un != basename(x[i])]
            x[i] <- paste0(path[length(path)], collapse = "/")
        }
    }
    return(x)
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
    ext <- Biostrings::strsplit(basename(x), split = "\\.")[[1]]
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
    filename <- Biostrings::strsplit(basename(x), split = ".([^.]*)$")[[1]]
    # filename <- filename[[-2]]
    return(filename)
}

## Usage:
# targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
# .getFileName(targetspath)

#################################################
## Return the logical, if the path is absolute ##
#################################################
is.fullPath <- function(x) {
    grepl("^(/|[A-Za-z]:|\\\\|~)", x)
}

## Usage:
# targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
# is.fullPath(targetspath)
# is.fullPath("./results")

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

##########################
## .parse_step function ##
##########################
## Internal parse function used in the subsetRmd function
.parse_step <- function(t_lvl, input_steps) {
    . <- NULL
    t_lvl_name <- names(t_lvl)
    input_steps <- unlist(input_steps %>% stringr::str_remove_all(" ") %>% stringr::str_split(",") %>% list())
    # single steps
    nocolon_steps <- input_steps[stringr::str_which(input_steps, "^[^:]+$")]
    lapply(nocolon_steps, function(x) if (!any(t_lvl_name %in% x)) stop("Step ", x, " is not found"))
    # dash linked steps
    dash_list <- NULL
    for (i in stringr::str_which(input_steps, ":")) {
        dash_step <- unlist(stringr::str_split(input_steps[i], ":"))
        dash_parse <- unlist(lapply(dash_step, function(x) {
            which(t_lvl_name %in% x) %>% ifelse(length(.) > 0, ., stop("Step ", x, " is not found"))
        })) %>%
            {
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

#############################
## .tryCatch function ##
#############################
.tryCatch <- function(x, file = NULL) {
    if (is.null(file)) file <- tempfile()
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

################################
## .checkSpecialChar function ##
################################
.checkSpecialChar <- function(x) {
    chunk_names_bad <- stringr::str_detect(x, "\\W")
    if (any(chunk_names_bad)) {
        stop(
            "Only letters, numbers, and '_' allowed for step_name. Invalid name:\n",
            paste0(x[chunk_names_bad], collapse = ", "),
            call. = FALSE
        )
    }
}

## Usage:
# .checkSpecialChar("name@")
# .checkSpecialChar("name test")

############################
## .renderMsg function ##
############################
.renderMsg <- function() {
    warn_flag <- getOption("spr_render_msg")
    if (isTRUE(warn_flag)) {
        return()
    }
    cat(
        crayon::green$bold("Done with workflow running, now consider rendering logs & reports\n"),
        crayon::blue("To render logs, run:    "), "sal <- renderLogs(sal)\n",
        crayon::blue("From command-line:      "), 'Rscript -e "sal = systemPipeR::SPRproject(resume = TRUE); sal = systemPipeR::renderLogs(sal)"\n',
        crayon::blue("To render reports, run: "), "sal <- renderReport(sal)\n",
        crayon::blue("From command-line:      "), 'Rscript -e "sal= s ystemPipeR::SPRproject(resume = TRUE); sal = systemPipeR::renderReport(sal)"\n',
        crayon::make_style("white")$bold("This message is displayed once per R session\n"),
        sep = ""
    )
    options(spr_render_msg = TRUE)
}

#################################
## .clusterRunResults function ##
################################
.clusterRunResults <- function(args.run, reg, nTargets) {
    logdir <- reg$file.dir
    file_log <- file.path(logdir, paste0(
        "sysargs2_log_",
        format(Sys.time(), "%b%d%Y_%H%Ms%S"),
        paste(sample(0:9, 4), collapse = "")
    ))
    id_save <- character()
    for (i in seq_along(1:nTargets)) {
        ## ind object created by the batchtools
        newsysargs <- batchtools::loadResult(reg = reg, id = i)
        id_sysargs <- names(newsysargs$cmdlist)
        id_save <- c(id_save, id_sysargs)
        ## update output slot
        args.run@output[id_sysargs] <- newsysargs$output[id_sysargs]
        ## Update status slot
        args.run@status$status.completed[id_sysargs, ] <- newsysargs$status$status.completed[id_sysargs, ]
        args.run@status$status.time[id_sysargs, ] <- newsysargs$status$status.time[id_sysargs, ]
        ## time
        args.run@status$status.time[id_sysargs, ]$time_start <- as.POSIXct(args.run@status$status.time[id_sysargs, ]$time_start, origin = "1970-01-01")
        args.run@status$status.time[id_sysargs, ]$time_end <- as.POSIXct(args.run@status$status.time[id_sysargs, ]$time_end, origin = "1970-01-01")
        ## Update files logs --> combining
        logs <- readLines(newsysargs$files$log)
        write(logs, file_log, append = TRUE, sep = "\n")
    }
    time_start <- sort(args.run$status$status.time[id_save, ]$time_start)[1]
    time_end <- tail(sort(args.run$status$status.time[id_save, ]$time_end), 1)
    args.run@status$total.time <- list(
        time_start = as.POSIXct(time_start, origin = "1970-01-01"),
        time_end = as.POSIXct(time_end, origin = "1970-01-01")
    )
    args.run@status$status.summary <- .statusSummary(args.run)
    args.run@files$log <- file_log
    return(args.run)
}


########################
## .tryRcode function ##
########################
.tryRcode <- function(command, envir) {
    warning <- error <- NULL
    value <- withCallingHandlers(
        tryCatch(
            eval(command, envir = envir),
            error = function(e) {
                error <<- conditionMessage(e)
                NULL
            }
        ),
        warning = function(w) {
            warning <<- append(warning, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    list(stdout = value, warning = warning, error = error)
}

###############################
## .updateAfterRunC function ##
###############################
.updateAfterRunC <- function(args, step) {
    args2 <- args
    conList <- args2$targets_connection[lengths(args2$targets_connection) != 0]
    conList_step <- sapply(conList, "[[", 1)
    for (l in seq_along(conList_step)) {
        if (step %in% conList_step[[l]]) {
            requiredUP <- names(conList)[[l]]
            for (s in requiredUP) {
                WF <- args2[s]
                WFstep <- names(stepsWF(WF))
                new_targets <- WF$targetsWF[[1]]
                col_out <- lapply(outfiles(args2), function(x) colnames(x))
                col_out_l <- sapply(names(col_out), function(x) list(NULL))
                for (i in names(col_out)) {
                    col_out_l[[i]] <- col_out[[i]][col_out[[i]] %in% WF$targets_connection[[WFstep]]$new_targets_col[[1]]]
                }
                col_out_l <- col_out_l[lapply(col_out_l, length) > 0]
                if (all(sapply(col_out_l, function(x) length(x) == 1))) {
                    col_out_df <- lapply(names(col_out_l), function(x) getColumn(args2, step = x, position = "outfiles", column = col_out_l[[x]]))
                    names(col_out_df) <- col_out_l
                    new_targets[as.character(col_out_l)] <- col_out_df
                } else {
                    col_out_df <- data.frame(args2[step]$outfiles[[step]][, col_out_l[[1]]])
                    names(col_out_df) <- col_out_l[[1]]
                    new_targets[as.character(col_out_l[[1]])] <- col_out_df
                }
                WF2 <- stepsWF(WF)[[1]]
                WF2 <- updateWF(WF2, new_targets = targets.as.list(data.frame(new_targets)), inputvars = WF2$inputvars, write.yaml = FALSE)
                ## Preserve outfiles
                WF2[["output"]] <- WF$stepsWF[[s]]$output
                args2 <- sysargslist(args2)
                args2$stepsWF[[WFstep]] <- WF2
                args2$targetsWF[[WFstep]] <- as(WF2, "DataFrame")
                rownames(args2$targetsWF[[WFstep]]) <- rownames(args$targetsWF[[WFstep]])
                args2$outfiles[[WFstep]] <- output.as.df(WF2)
                rownames(args2$outfiles[[WFstep]]) <- rownames(args$outfiles[[WFstep]])
                args2$statusWF[[WFstep]] <- WF2$status
                rownames(args2$statusWF[[WFstep]]$status.completed) <- rownames(args$outfiles[[WFstep]])
                rownames(args2$statusWF[[WFstep]]$status.time) <- rownames(args$outfiles[[WFstep]])
                args2 <- as(args2, "SYSargsList")
            }
        } else {
            do <- "donothing"
        }
    }
    return(args2)
}

#############################
## .statusSummary function ##
#############################
.statusSummary <- function(args) {
    if (inherits(args, "SYSargs2")) {
        step.status.summary <- args$status$status.completed
    } else if (inherits(args, "data.frame")) {
        step.status.summary <- args[5:ncol(args)]
    } else {
        stop("Argument 'args' needs to be assigned an object of class 'SYSargs2' or 'data.frame'.")
    }
    if ("Error" %in% unlist(unique(step.status.summary))) {
        step.status <- "Error"
    } else if ("Warning" %in% unlist(unique(step.status.summary))) {
        step.status <- "Warning"
    } else if ("Success" %in% unlist(unique(step.status.summary))) {
        step.status <- "Success"
    } else if ("Pending" %in% unlist(unique(step.status.summary))) {
        step.status <- "Pending"
    } else if (is.null(step.status.summary)) {
        step.status <- "Pending"
    }
    return(step.status)
}

########################
## .prepareRmdPlot ##
########################
.prepareRmdPlot <- function(sysargs, dir_log) {
    out_path <- file.path(dir_log, "log_plot.html")
    plotWF(sysargs,
           out_format = "html", out_path = out_path, rmarkdown = TRUE,
           in_log = TRUE, rstudio = TRUE, plot_ctr = FALSE
    )
    # modify HTML content
    if (!file.exists(out_path)) stop("Cannot create the workflow plot for logs at\n", out_path)
    plot_content <- readLines(out_path)
    out_path
}

################################
## .dirProject function ##
################################
.dirProject <- function(projPath, data, param, results, silent) {
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
                "It is required to have the project structure in place. We can create the directories now, and you can find more information by reading the `?SPRProject` help file. ", "\n",
                "\n", "There is no directory called", "\n", paste(names(create), collapse = " OR ", sep = "\n"), "\n", "\n",
                "Would you like to create this directory now? Type a number: \n 1. Yes \n 2. No \n"
            ))
        } else {
            ## For an non-interactive session
            dir_create <- "1"
        }
        for (i in seq_along(create)) {
            if (dir_create == "1") {
                dir.create(create[[i]], recursive = TRUE)
                if (silent != TRUE) cat("Creating directory: ", create[[i]], "\n")
            } else if (dir_create == 2) {
                stop("Aborting project creation. Find more information by reading the `?SPRProject` help file.", call. = FALSE)
            }
        }
    }
    project <- list(
        project = projPath,
        data = file.path(data),
        param = file.path(param),
        results = file.path(results)
    )
    return(project)
}

#############################
## .statusPending function ##
#############################
.statusPending <- function(args) {
    status.pending <- check.output(args)
    ## function
    .statusSYSargs2 <- function(args, status.pending) {
        pending <- sapply(args$files$steps, function(x) list(x = "Pending"))
        pending <- data.frame(matrix(unlist(pending), ncol = length(pending), byrow = TRUE), stringsAsFactors = TRUE)
        colnames(pending) <- args$files$steps
        status.pending <- cbind(status.pending, pending)
        rownames(status.pending) <- status.pending$Targets
        status.pending[c(2:4)] <- sapply(status.pending[c(2:4)], as.numeric)
        status.pending[c(5:ncol(status.pending))] <- sapply(status.pending[c(5:ncol(status.pending))], as.character)
        # status.time <- data.frame(matrix(, nrow=nrow(status.pending),  ncol = length(pending)))
        status.time <- status.pending[, 1, drop = FALSE]
        # status.time <- cbind(status.time)
        status.time <- cbind(status.time, time_start = NA, time_end = NA)
        status.time$time_start <- .POSIXct(status.time$time_start)
        status.time$time_end <- .POSIXct(status.time$time_end)
        rownames(status.time) <- status.pending$Targets
        pendingList <- list(
            status.summary = .statusSummary(status.pending),
            status.completed = status.pending, status.time = status.time
        )
    }
    if (inherits(args, "SYSargsList")) {
        for (i in seq_along(status.pending)) {
            # print(i)
            if (inherits(args$stepsWF[[i]], "SYSargs2")) {
                status.pending[[i]] <- .statusSYSargs2(args$stepsWF[[i]], status.pending[[i]])
                # print(status.pending[[i]])
            }
        }
        pendingList <- status.pending
    } else if (inherits(args, "SYSargs2")) {
        pendingList <- .statusSYSargs2(args, status.pending)
    }
    return(pendingList)
}

##########################
## .outList2DF function ##
##########################
.outList2DF <- function(args) {
    if (inherits(args, "list")) {
        args <- as(args, "SYSargsList")
        out <- sapply(names(stepsWF(args)), function(x) list(NULL))
        for (i in seq_along(stepsWF(args))) {
            l_out <- output(stepsWF(args)[[i]])
            out[[i]] <- S4Vectors::DataFrame(matrix(unlist(l_out), nrow = length(l_out), byrow = TRUE))
            # out <- S4Vectors::DataFrame(as.data.frame(do.call(rbind, l_out)))
            colnames(out[[i]]) <- stepsWF(args)[[i]]$files$output_names
        }
    } else if (inherits(args, "SYSargs2")) {
        l_out <- output(args)
        # sapply(l_out, function(x) length(x))
        out <- S4Vectors::DataFrame(matrix(unlist(l_out), nrow = length(l_out), byrow = TRUE))
        # out <- S4Vectors::DataFrame(as.data.frame(do.call(rbind, l_out)))
        colnames(out) <- args$files$output_names
    }
    return(out)
}

#############################
## .outputTargets function ##
#############################
.outputTargets <- function(args, fromStep, index = 1, toStep, replace = c("FileName")) {
    if (!inherits(args, "SYSargsList")) stop("Argument 'args' needs to be assigned an object of class 'SYSargsList'")
    outputfiles <- outfiles(args[fromStep])[[1]]
    if (length(targetsWF(args)[[fromStep]]) > 0) {
        df <- targetsWF(args)[[fromStep]]
        df[replace] <- outputfiles
    }
    return(df)
}
