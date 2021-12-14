#######################
## LineWise function ##
#######################
LineWise <- function(code, step_name = "default", codeChunkStart = integer(),
                     rmdPath = character(), dependency = "",
                     run_step = "mandatory",
                     run_session = "management") {
    ## used in `importWF`
    on.exit({options(linewise_importing = FALSE)})
    ## check options
    run_step <- match.arg(run_step, c("mandatory", "optional"))
    run_session <- match.arg(run_session, c("management", "compute"))
    ## Step name
    if (step_name == "default") {
        step_name <- "Step_x"
    } else {
        .checkSpecialChar(step_name)
        step_name <- step_name
    }
    ## dependency
    dependency <- list(dependency)
    names(dependency) <- step_name
    ## RunInfo
    runInfo <- list(runOption = list(list(directory = FALSE, run_step = run_step,
                                          run_session = run_session)))
    names(runInfo$runOption) <- step_name
    ## status
    step_status <- list(
        status.summary = "Pending",
        status.completed = data.frame(Step = step_name, status.summary = "Pending"),
        status.time = data.frame()
    )
    ## codeLine
    if (!getOption("linewise_importing", TRUE)) {
        codeLine <- parse(text = gsub("^\\{|\\}$", "", deparse(substitute(code))))
    } else if(getOption("linewise_importing", TRUE)) {
        codeLine <- parse(text = gsub("^\\{|\\}$", "", deparse(substitute(code))))
    }
    line <- list(
        codeLine = codeLine,
        codeChunkStart = codeChunkStart,
       # rmdPath = rmdPath,
        stepName = step_name,
        dependency = dependency,
        status = step_status,
        files = list(rmdPath = rmdPath),
        runInfo = runInfo
    )
    return(as(line, "LineWise"))
}
## Usage:
# xx <- LineWise({
#   x <- 1+1;  y <- 100
# })
# plot <- LineWise({
#   plot(iris)
# })


########################
## importRmd function ##
########################
importWF <- function(sysargs, file_path, ignore_eval = TRUE, verbose = TRUE) {
    on.exit({
        options(linewise_importing = FALSE)
        options(spr_importing = FALSE)
        options(importwf_options = NULL)
    })
    stopifnot(is.character(file_path) && length(file_path) == 1)
    if (!stringr::str_detect(file_path, "\\.[Rr]md$")) stop("File must be .Rmd, or .rmd ending.")
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    stopifnot(is.logical(ignore_eval) && length(ignore_eval) == 1)
    ## start
    if (verbose) cat(crayon::blue$bold("Reading Rmd file"))
    df <- parseRmd(file_path, ignore_eval = ignore_eval, verbose = verbose)
    df$dep <- lapply(df$dep, function(x) ifelse(x == "", NA, x))
    names(df$dep) <- df$step_name
    ## create a new env for sysarges to eval
    sysargs_env <- new.env()
    sal_imp <- sysargs
    sal_imp <- as(sal_imp, "list")
    ## adding steps
    for (i in seq_along(df$spr)) {
        if (verbose) cat(crayon::blue$bold("Now importing step '", df$step_name[i], "' \n", sep = ""))
        # if (df$spr[i] == "r") {
        #     options(linewise_importing = TRUE)
        #     line_obj <- LineWise(code = df$code[i], step_name = df$step_name[i],
        #                          codeChunkStart = df$start[i], rmdPath = file_path,
        #                          dependency = df$dep[[i]])
        #     sal_imp$stepsWF[[df$step_name[i]]] <- line_obj
        #     sal_imp$statusWF[[df$step_name[i]]] <- list(
        #         status.summary = "Pending",
        #         status.completed = data.frame(Step = df$step_name[i], status.summary = "Pending"),
        #         status.time = data.frame()
        #     )
        #     sal_imp$targetsWF[[df$step_name[i]]] <- S4Vectors::DataFrame()
        #     sal_imp$outfiles[[df$step_name[i]]] <- S4Vectors::DataFrame()
        #     sal_imp$SE[df$step_name[i]] <- list(NULL)
        #     sal_imp$dependency[[df$step_name[i]]] <- df$dep[[i]]
        #     sal_imp$targets_connection[df$step_name[i]] <- list(NULL)
        #     sal_imp$runInfo[["runOption"]][df$step_name[i]] <- list(
        #         list(directory = FALSE, run_step = df$req[i], run_session = df$session[i], 
        #              rmd_line = paste(df[i,2:3], collapse = ":")))
        # } else if (df$spr[i] == "sysargs") {
            options(spr_importing = TRUE)
            options(importwf_options = c(df$step_name[i], df$dep[i]))
            sal_imp <- as(sal_imp, "SYSargsList")
            salname <- sub("[\\)].*", "", sub(".*(appendStep\\()", "", df$code[i]))
            assign(salname, sal_imp, sysargs_env)
            args <- eval(parse(text = df$code[i]), envir = sysargs_env)
          #  if (!inherits(args, "SYSargsList")) stop("Cannot import this step. It is not returning a `SYSargsList` object.")
            appendStep(sal_imp) <- args
            sal_imp <- as(sal_imp, "list")
            sal_imp$runInfo[["runOption"]][[df$step_name[i]]][['rmd_line']] <- paste(df[i,2:3], collapse = ":")
        #}
    }
    sal_imp[["projectInfo"]]$rmd_file <- file_path
    # sal_imp[["projectInfo"]]$rmd_lines <- df[,1:3]
    sal_imp <- as(sal_imp, "SYSargsList")
    sysargslist <- file.path(sal_imp$projectInfo$project, sal_imp$projectInfo$sysargslist)
    write_SYSargsList(sal_imp, sysargslist, silent = TRUE)
    return(sal_imp)
}

## Usage:
# file_path <- system.file("extdata", "spr_simple_wf.Rmd", package="systemPipeR")
# sal <- SPRproject(overwrite = TRUE)
# sal <- importWF(sal, file_path)

########################
## parseRmd function ##
########################
########################
parseRmd <- function(file_path, ignore_eval = TRUE, verbose = FALSE) {
    # assertions
    stopifnot(is.character(file_path) && length(file_path) == 1)
    stopifnot(is.logical(ignore_eval) && length(ignore_eval) == 1)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    # read file
    lines <- readLines(file_path)
    # get chunk numbers --------
    chunk_start <- lines %>% stringr::str_which("^```\\{.*\\}.*")
    chunk_end <- lines %>% stringr::str_which("^```[[:blank:]]{0,}$")
    if (length(chunk_start) != length(chunk_end)) {
        stop("unmatched number of code chunk starts and ends")
    }
    for (i in seq_along(chunk_start)[-length(chunk_end)]) {
        if (chunk_start[i + 1] <= chunk_end[i]) {
            stop(paste(
                "A code chunk does not end, chunk line:\n",
                chunk_start[i + 1]
            ))
        }
    }
    # create df
    df <- data.frame(
        # cols export
        step_name = "",
        start = chunk_start,
        end = chunk_end,
        dep = NA,
        code = NA,
        spr = NA,
        req = NA,
        session = NA,
        has_run = FALSE,
        success = FALSE,
        sample_pass = 0,
        sample_warn = 0,
        sample_error = 0,
        sample_total = 0,
        log_path = "",
        time_start = Sys.time(),
        time_end = Sys.time() + 1,
        # cols will be removed?
        opt_text = ""
    )
    # get r chunks -----
    r_chunk <- lines[chunk_start] %>% stringr::str_detect("^```\\{r")
    # get chunk start and end
    if (verbose) message(crayon::blue("\n", "---- Actions ----"))
    if (verbose && length(chunk_start[!r_chunk])) message("Ignore none-R chunks at line: ", paste0(chunk_start[!r_chunk], collapse = ", "))
    df <- df[r_chunk, ]
    if (nrow(df) == 0) stop("No valid R code chunk left")
    # get chunk text
    chunk_text <- lines[df$start] %>%
        stringr::str_remove("^```\\{r[ ]?") %>%
        stringr::str_remove("\\}[[:blank:]]{0,}$")
    # get chunk names
    chunk_names <- stringr::str_match(chunk_text, "^([^,])+")[, 1]
    # chunks with no name, no option will be ""
    chunk_names[is.na(chunk_names)] <- ""
    # chunks with no name but only options
    opt_only_chunk <- stringr::str_which(chunk_names, "=")
    opt_only_chunk_text <- chunk_names[opt_only_chunk]
    # safe to clean chunk names of chunks with no name but only options
    chunk_names[opt_only_chunk] <- ""
    # get opt string
    opt_text <- stringr::str_match(chunk_text, ",(.*)")[, 2] %>% stringr::str_remove("^[ ]{0,}")
    # preappend  no name but only options chunks options to other options
    opt_text[opt_only_chunk] <- paste0(opt_only_chunk_text, ", ", opt_text[opt_only_chunk])
    # clean ", " end options
    opt_text <- stringr::str_remove(opt_text, ",[[:blank:]]{0,}$")
    # clean empty options
    opt_text[is.na(opt_text)] <- ""
    # update chunk name
    df$step_name <- chunk_names
    df$opt_text <- opt_text
    # get spr option
    if (verbose) message("Checking chunk SPR option")
    opt_spr <- stringr::str_match(df$opt_text, "spr[ ]{0,}=[ ]{0,}(['\"]\\w+['\"])")[, 2] %>% stringr::str_remove_all('\'|"')
    ## update df
    if (verbose && any(is.na(opt_spr))) message("Ignore non-SPR chunks: ", paste0(df$start[is.na(opt_spr)], collapse = ", "))
    df <- df[!is.na(opt_spr), ]
    df$spr <- opt_spr[!is.na(opt_spr)]
    ## enforce sysarg or r option only for spr
    bad_spr <- !df$spr %in% c("sysargs", "r")
    if (any(bad_spr)) {
        stop(
            "Bad SPR option for chunk at line:\n",
            paste(df$start[bad_spr], df$spr[bad_spr], collapse = ", "),
            "\nOnly 'sysargs' or 'r' allowed"
        )
    }
    # get reqiurement options
    if (verbose) message("Checking chunk spr.req option")
    spr_req <- stringr::str_match(df$opt_text, "spr\\.req[ ]{0,}=[ ]{0,}['\"]((\\w+[;]{0,}[ ]{0,})+)['\"]")[, 2] %>%
        stringr::str_remove_all('\'|"') %>%
        stringr::str_remove_all("^spr\\.req[ ]{0,}=[ ]{0,}")
    ## update df
    spr_req[is.na(spr_req)] <- "mandatory"
    spr_req[spr_req == "m"] <- "mandatory"
    spr_req[spr_req == "o"] <- "optional"

    df$req <- spr_req
    ## enforce sysarg or r option only for spr
    bad_req <- !df$req %in% c("mandatory", "optional")
    if (any(bad_req)) {
        stop(
            "Bad spr.req option for chunk at line:\n",
            paste(df$start[bad_req], collapse = ", "),
            "\nOnly 'm' or 'o' allowed"
        )
    }
    # get session options
    if (verbose) message("Checking chunk spr.ses option")
    spr_ses <- stringr::str_match(df$opt_text, "spr\\.ses[ ]{0,}=[ ]{0,}['\"]((\\w+[;]{0,}[ ]{0,})+)['\"]")[, 2] %>%
        stringr::str_remove_all('\'|"') %>%
        stringr::str_remove_all("^spr\\.ses[ ]{0,}=[ ]{0,}")
    ## update df
    spr_ses[is.na(spr_ses)] <- "management"
    spr_ses[spr_ses == "r"] <- "management"
    spr_ses[spr_ses == "c"] <- "compute"

    df$session <- spr_ses
    ## enforce sysarg or r option only for spr
    bad_ses <- !df$session %in% c("management", "compute")
    if (any(bad_ses)) {
        stop(
            "Bad spr.req option for chunk at line:\n",
            paste(df$start[bad_ses], collapse = ", "),
            "\nOnly 'c' or 'r' allowed"
        )
    }
    # get eval
    if (verbose) message("Checking chunk eval values")
    opt_eval <- stringr::str_match(df$opt_text, "eval[ ]{0,}=[ ]{0,}(\\w+)")[, 2] %>% stringr::str_detect("^TRUE|T")
    # eval unspecified will be TRUE
    opt_eval[is.na(opt_eval)] <- TRUE
    # overwrite if ignore eval param is TRUE
    if (ignore_eval) opt_eval <- rep(TRUE, length(opt_eval))
    ## update df
    if (verbose && any(!opt_eval)) message("Ignore chunks with 'eval' are off: ", paste0(df$start[!opt_eval], collapse = ", "))
    df <- df[opt_eval, ]
    if (nrow(df) == 0) stop("No valid SPR code chunk left")
    # check step names
    if (verbose) message("Resolve step names")
    chun_name_empty <- df$step_name %in% c("", NA)
    if (verbose && any(chun_name_empty)) message("Give default names to chunks with empty names at line: ", paste0(df$start[chun_name_empty], collapse = ", "))
    df$step_name[chun_name_empty] <- paste0("unnamed_step", seq(sum(chun_name_empty)))
    if (verbose) message("Check duplicated step names")
    if (any(duplicated(df$step_name))) {
        stop("Duplicated step names not allowed: ", paste0(df$step_name[duplicated(df$step_name)], collapse = ", "))
    }
    chunk_names_bad <- stringr::str_detect(df$step_name, "\\W")
    if (any(chunk_names_bad)) {
        stop(
            "Only letters, numbers, and '_' allowed for step names. Invalid names:\n",
            paste0(df$step_name[chunk_names_bad], collapse = ", ")
        )
    }
    # get dependencies
    if (verbose) message("Checking chunk dependencies")
    ## first use previous step as dep for all and use "" for first step
    df$dep[seq(2, nrow(df))] <- df$step_name[seq(1, nrow(df) - 1)]
    df$dep[1] <- ""
    ## parse user provided dep
    spr_dep <- stringr::str_match(df$opt_text, "spr\\.dep[ ]{0,}=[ ]{0,}['\"]((\\w+[;]{0,}[ ]{0,})+)['\"]")[, 2] %>%
        stringr::str_remove_all('\'|"') %>%
        stringr::str_remove_all("^spr\\.dep[ ]{0,}=[ ]{0,}")
    if (verbose && any(is.na(spr_dep))) message("Use the previous step as dependency for steps without 'spr.dep' options: ", paste0(df$start[is.na(spr_dep)], collapse = ", "))
    spr_dep <- stringr::str_split(spr_dep, ";")
    spr_dep <- lapply(spr_dep, function(x) stringr::str_remove_all(x, "^[ ]{0,}") %>% stringr::str_remove_all("[ ]{0,}$"))

    lapply(seq_along(spr_dep), function(i) {
        lapply(spr_dep[[i]], function(x) {
            if (!x %in% c(df$step_name, NA)) stop("Step '", df$step_name[i], "'s dependency '", x, "' not found.")
        })
    })
    df$dep <- lapply(seq_along(spr_dep), function(i) {
        if (is.na(spr_dep[[i]][1])) {
            df$dep[i]
        } else {
            spr_dep[[i]]
        }
    })
    if (df$dep[[1]] != "") {
        warning("First step has the dependency of '", paste0(df$dep[[1]], collapse = ", "), "'. Usually the first step has no depenency", immediate. = TRUE)
        message("Clear dependency of step 1")
        df$dep[[1]] <- ""
    }
    # get code
    if (verbose) message("Parse chunk code")
    code <- lapply(seq_along(df$start), function(x) {
        paste0(lines[(df$start[x] + 1):(df$end[x] - 1)], collapse = "\n")
    }) %>% unlist()
    df$code <- code
    # create log path
    df$log_path <- paste0("#", tolower(stringr::str_replace_all(df$step_name, "[ ]", "-")))
    if (verbose) message(crayon::blue("---- Succes! Create output ----"))
    df
}
## Usage:
# rmdpath <- system.file("extdata", "spr_simple_wf.Rmd", package="systemPipeR")
# df <- parseRmd(rmdpath)
