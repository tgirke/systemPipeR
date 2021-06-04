# library(magrittr)

#' parseRmd
#'
#' @param file_path string, file path of the workflow file
#' @param ignore_eval bool, treat all R chunks' `eval` option as `TRUE` in workflow
#' Rmd file even if some chunks have `eval=FALSE`?
#' @param verbose bool, print out verbose message while function running?
#'
#' @return a dataframe
#' step_name: string, step ID, unique names
#' start: numeric, chunk start line number
#' end: numeric, chunk end line number
#' dep: list of character vector, other step(s) this step is depend on
#' code: character string of length 1, code inside current step, all code is collapse in one string by "\n" for multiple lines
#' spr: one of 'sysargs' or 'r'
#' run:
#' succes:
#' @export
#'
#' @examples
#' @importFrom stringr str_which str_detect str_remove str_match str_remove_all str_split
#' @importFrom magrittr %>%
parseRmd <- function(file_path, ignore_eval = TRUE, verbose = TRUE){
    # assertions
    stopifnot(is.character(file_path) && length(file_path) == 1)
    stopifnot(is.logical(ignore_eval) && length(ignore_eval) == 1)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    # read file
    lines <- readLines(file_path)
    
    # get chunk numbers --------
    chunk_start <- lines %>% stringr::str_which("^```\\{.*\\}.*")
    chunk_end <- lines %>% stringr::str_which("^```[[:blank:]]{0,}$")
    if (length(chunk_start) != length(chunk_end))
        stop("unmatched number of code chunk starts and ends")
    
    for (i in seq_along(chunk_start)[-length(chunk_end)]){
        if (chunk_start[i+1] <= chunk_end[i])
            stop(paste("A code chunk does not end: chunk line",
                       chunk_start[i+1]))
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
        run = FALSE,
        success = FALSE,
        # cols will be removed
        opt_text = ""
    )
    
    # get r chunks -----
    r_chunk <- lines[chunk_start] %>% stringr::str_detect("^```\\{r")
    # get chunk start and end
    if (verbose) message(crayon::blue("---- Actions ----"))
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
    
    # get spr chunks
    if(verbose) message("Checking chunk SPR option")
    opt_spr <- stringr::str_match(df$opt_text, "spr[ ]{0,}=[ ]{0,}(['\"]\\w+['\"])")[, 2] %>% stringr::str_remove_all('\'|"')
    ## update df
    if (verbose && any(is.na(opt_spr))) message("Ignore non-SPR chunks: ", paste0(df$start[is.na(opt_spr)], collapse = ", "))
    df <- df[!is.na(opt_spr), ]
    df$spr <- opt_spr[!is.na(opt_spr)]
    ## enforce sysargs or r option only for spr
    bad_spr <- !df$spr %in% c("sysargs", "r")
    if(any(bad_spr)) stop(
        "Bad SPR option for chunk at line:\n",
        paste(df$chunk_start[bad_spr], df$spr[bad_spr], collapse = ", "),
        "\nOnly 'sysargs' or 'r' allowed"
    )
    
    # get eval
    if(verbose) message("Checking chunk eval values")
    opt_eval <- stringr::str_match(df$opt_text, "eval[ ]{0,}=[ ]{0,}(\\w+)")[, 2] %>% stringr::str_detect("^TRUE|T")
    # eval unspecified will be TRUE
    opt_eval[is.na(opt_eval)] <- TRUE
    # overwrite if ignore eval param is TRUE
    if(ignore_eval) opt_eval <- rep(TRUE, length(opt_eval))
    ## update df
    if (verbose && any(!opt_eval)) message("Ignore chunks with 'eval' are off: ", paste0(df$start[!opt_eval], collapse = ", "))
    df <- df[opt_eval, ]
    if (nrow(df) == 0) stop("No valid SPR code chunk left")
    
    # check step names
    if(verbose) message("Resolve step names")
    chun_name_empty <- df$step_name %in% c("", NA)
    if(verbose && any(chun_name_empty)) message("Give default names to chunks with empty names at line: ", paste0(df$start[chun_name_empty], collapse = ", "))
    df$step_name[chun_name_empty] <- paste0("unnamed_step", seq(sum(chun_name_empty)))
    if(verbose) message("Check duplicated step names")
    if (any(duplicated(df$step_name))) {
        stop("Duplicated step names not allowed: ", paste0(df$step_name[duplicated(df$step_name)], collapse = ", "))
    }
    chunk_names_bad <- stringr::str_detect(df$step_name, "\\W")
    if (any(chunk_names_bad)) {
        stop("Only letters, numbers, and '_' allowed for step names. Invalid names:\n",
             paste0(df$step_name[chunk_names_bad], collapse = ", "))
    }
    
    # get dependencies
    if(verbose) message("Checking chunk dependencies")
    ## first use previous step as dep for all and use "" for first step
    df$dep[seq(2, nrow(df))] <- df$step_name[seq(1, nrow(df) -1)]
    df$dep[1] <- ""
    ## parse user provided dep
    spr_dep <- stringr::str_match(df$opt_text, "spr\\.dep[ ]{0,}=[ ]{0,}['\"]((\\w+[;]{0,}[ ]{0,})+)['\"]")[, 2] %>%
        stringr::str_remove_all('\'|"') %>%
        stringr::str_remove_all('^spr\\.dep[ ]{0,}=[ ]{0,}')
    if(verbose && any(is.na(spr_dep))) message("Use the previous step as dependency for steps without 'spr.dep' options: ", paste0(df$start[is.na(spr_dep)], collapse = ", "))
    spr_dep <- stringr::str_split(spr_dep, ";")
    spr_dep <- lapply(spr_dep, function(x) stringr::str_remove_all(x, "^[ ]{0,}") %>% stringr::str_remove_all("[ ]{0,}$"))
    
    lapply(seq_along(spr_dep), function(i){
        lapply(spr_dep[[i]], function(x){
            #print(x)
            if(!x %in% c(df$step_name, NA)) stop("Step '", df$step_name[i], "'s dependency '", x, "' not found.")
        })
    })
    
    df$dep <- lapply(seq_along(spr_dep), function(i){
        if(is.na(spr_dep[[i]][1])) df$dep[i]
        else spr_dep[[i]]
    })
    
    if(df$dep[[1]] != "") {
        warning("First step has the dependency of '", paste0(df$dep[[1]], collapse = ", "), "'. Usually the first step has no depenency", immediate. = TRUE)
        message("Clear dependency of step 1")
        df$dep[[1]] = ""
    }
    # get code
    if(verbose) message("Parse chunk code")
    code <- lapply(seq_along(df$start), function(x) {
        paste0(lines[(df$start[x] + 1): (df$end[x] -1)], collapse = "\n")
    }) %>% unlist()
    
    if(verbose) message(crayon::blue("---- Succes! Create output ----"))
    
    df$code <- code
    df
}


#

## Usage:
# library(systemPipeRdata)
# rmdPath <- system.file("extdata/workflows/rnaseq", "systemPipeRNAseq.Rmd", package="systemPipeRdata")
# targets <- system.file("extdata", "targets.txt", package="systemPipeR")
# rmdPath <- "systemPipeRNAseq.Rmd"
# parseRmd(rmdPath, ignore_eval = FALSE)
# 
# parseRmd("../../rnaseq/rna1.Rmd")
# # a <- parseRmd("rna1.Rmd")
# file_path <- "rna1.Rmd"

LineWise <- function(code, stepName="default", codeChunkStart=integer(), rmdPath=character(), dependency=character()){
    if(stepName=="default"){
        stepName <- "Step_x"
    } else {
        stepName <- stepName
    }
    line <- list(
        codeLine = parse(text=code),
        codeChunkStart = codeChunkStart,
        rmdPath = rmdPath,
        stepName=stepName,
        dependency=dependency
    )
    return(as(line, "LineWise"))
}


importRmd <- function(file_path, ignore_eval = TRUE, verbose = TRUE){
    df <- parseRmd(file_path, ignore_eval = ignore_eval, verbose = verbose)
    sal <- SYSargsList(silent = TRUE)
    sal <- sysargslist(sal)
    for(i in seq_along(df$spr)){
        if(df$spr[i]=="r"){
            line_obj <- LineWise(df$code[i], codeChunkStart=df$start[i], file_path)
            sal$stepsWF[[df$step_name[i]]] <- line_obj
            sal$dependency[[df$step_name[i]]] <- df$dep[[i]]
            sal$statusWF[[df$step_name[i]]] <- list(status="Pending")
            } else if(df$spr[i]=="sysargs"){
            args <- eval(parse(text =df$code[i]), envir = globalenv())
            # suppressMessages(
            # args <- eval(parse(text=df$code[i])))
            sal$stepsWF[[df$step_name[i]]] <- args$stepsWF[[1]]
            sal$dependency[[df$step_name[i]]] <- df$dep[[i]]
            sal$statusWF[[df$step_name[i]]] <- list(status="Pending")
            sal$outfiles[[df$step_name[i]]] <- args$outfiles[[1]]
            sal$targetsWF[[df$step_name[i]]] <- args$targetsWF[[1]]
     }
    }
    return(as(sal, "SYSargsList"))
}
# file_path <- "../systemPipeR/inst/extdata/systemPipeTEST.Rmd"
# sal <- importRmd(file_path)
