#######################
## LineWise function ##
#######################
LineWise <- function(code, step_name = "default", codeChunkStart = integer(),
                     rmdPath = character(), dependency = NA,
                     run_step = "mandatory",
                     run_session = "management",
                     run_remote_resources = NULL) {
    ## used in `importWF`
    on.exit({
        options(linewise_importing = FALSE)
    })
    ## check options
    run_step <- match.arg(run_step, c("mandatory", "optional"))
    run_session <- match.arg(run_session, c("management", "compute"))
    if (!is.null(run_remote_resources)) {
        if (!inherits(run_remote_resources, "list")) {
            stop("Argument 'run_remote_resources' needs to be assigned an object of class 'list'")
        }
    }
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
    if (!is.null(run_remote_resources)) {
        if (run_session == "management") {
            message("Please note that the '", step_name, "' run_session option '", run_session, "' was replaced with 'compute' because run_remote_resources was available.")
            run_session <- "compute"
        }
    }
    runInfo <- list(runOption = list(list(
        directory = FALSE, run_step = run_step,
        run_session = run_session,
        run_remote_resources = run_remote_resources
    )))
    names(runInfo$runOption) <- step_name
    ## status
    step_status <- list(
        status.summary = "Pending",
        status.completed = data.frame(Step = step_name, status.summary = "Pending"),
        status.time = data.frame()
    )
    ## codeLine
    if (length(deparse(substitute(code))) == 1) {
        codeLine <- parse(text = sub("^\\{|\\}$", "", deparse(substitute(code))))
    } else {
        codeLine <- deparse(substitute(code))
        codeLine[1] <- sub("^\\{", "", codeLine[1])
        codeLine[length(codeLine)] <- sub("\\}$", "", codeLine[length(codeLine)])
        codeLine <- parse(text = codeLine)
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
    ## create a new env for sysargs to eval
    sysargs_env <- new.env()
    sal_imp <- sysargs
    sal_imp <- as(sal_imp, "list")
    ## adding steps
    for (i in seq_along(df$spr)) {
        if (verbose) cat(crayon::blue$bold("Now importing step '", df$step_name[i], "' \n", sep = ""))
        options(spr_importing = TRUE)
        sal_imp <- as(sal_imp, "SYSargsList")
        salname <- sub("[\\)].*", "", sub(".*(appendStep\\()", "", df$code[i]))
        assign(salname, sal_imp, sysargs_env)
        args <- eval(parse(text = df$code[i]), envir = sysargs_env)
        appendStep(sal_imp) <- args
        sal_imp <- as(sal_imp, "list")
        sal_imp$runInfo[["runOption"]][[df$step_name[i]]][["rmd_line"]] <- paste(df[i, 2:3], collapse = ":")
    }
    sal_imp[["projectInfo"]]$rmd_file <- file_path
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
            stop(
                "A code chunk does not end, chunk line:\n",
                chunk_start[i + 1]
            )
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
    # check eval
    if (verbose) message("Checking chunk eval values")
    opt_eval <- stringr::str_match(df$opt_text, "eval[ ]{0,}=[ ]{0,}(['\"]{0,1}\\w+['\"]{0,1})")[, 2]
    # eval unspecified will be TRUE
    opt_eval[is.na(opt_eval)] <- TRUE
    if (any(stringr::str_detect(opt_eval, "['\"]"))) {
        stop(
            "Line ",
            paste(df$start[stringr::str_detect(opt_eval, "['\"]")], collapse = ", "),
            " has a character value for `eval` option, only TRUE, T, FALSE, F are valid"
        )
    }
    opt_eval <- stringr::str_detect(opt_eval, "^TRUE|T")
    # overwrite if ignore eval param is TRUE
    if (ignore_eval) opt_eval <- rep(TRUE, length(opt_eval))
    ## update df
    if (verbose && any(!opt_eval)) message("Ignore chunks with 'eval' are FALSE or invalid values: ", paste0(df$start[!opt_eval], collapse = ", "))
    df <- df[opt_eval, ]
    if (nrow(df) == 0) stop("No evaluated code chunk detected")
    # get spr option
    if (verbose) message("Checking chunk SPR option")
    opt_spr <- stringr::str_match(df$opt_text, "spr[ ]{0,}=[ ]{0,}(\\w+)")[, 2] %>% stringr::str_detect("^TRUE|T")
    if (all(is.na(opt_spr))) stop("No valid SPR code chunk detected, remember to add `spr=TRUE` in chunk options")
    ## update df
    if (verbose && any(is.na(opt_spr))) message("Ignore non-SPR chunks: ", paste0(df$start[is.na(opt_spr)], collapse = ", "))
    # eval unspecified will be TRUE
    opt_spr[is.na(opt_spr)] <- FALSE
    df <- df[opt_spr, ]
    if (nrow(df) == 0) stop("No SPR code chunk left")
    if (verbose) message("Parse chunk code")
    code <- lapply(seq_along(df$start), function(x) {
        paste0(lines[(df$start[x] + 1):(df$end[x] - 1)], collapse = "\n")
    }) %>% unlist()
    df$code <- code
    df
}

## Usage:
# rmdpath <- system.file("extdata", "spr_simple_wf.Rmd", package="systemPipeR")
# df <- parseRmd(rmdpath)
