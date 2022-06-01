############################
## ## Version 2 Functions  ##
############################

############# Exports ###########

# createParam2 # the same as createParam

printParam2 <- function(
        sysargs,
        base = FALSE, args = FALSE, inputs = FALSE,
        outputs = FALSE, stdout = FALSE, raw_cmd = FALSE, all = TRUE
    ) {
    obj <- if(inherits(sysargs, "SYSargs2")) cmdToCwl(sysargs) else if(inherits(sysargs, "cwlParse")) sysargs else {
        stop("sysargs must be `SYSargs2` or `cwlParse` object")
    }
    if(inherits(obj, "v1")) stop("input sysargs built on syntax version 1, please use `printParam")
    .catCmdAll2(
        obj, base = base, args = args, inputs = inputs,
        outputs = outputs, stdout = stdout, raw_cmd = raw_cmd, all = all
    )
}

appendParam2 <- function(
        sysargs, x,
        position = c("inputs", "args", "outputs"),
        after = NULL,
        verbose = FALSE
    ) {
    obj <- if(inherits(sysargs, "SYSargs2")) cmdToCwl(sysargs) else if(inherits(sysargs, "cwlParse")) sysargs else {
        stop("sysargs must be `SYSargs2` or `cwlParse` object")
    }
    if(inherits(obj, "v1")) stop("input sysargs built on syntax version 1, please use `appendParam`")
    .cmdAppend2(
        obj, x = x,
        position = position,
        after = after,
        verbose = verbose
    )
}

replaceParam2 <- function(
        sysargs, x, index=NULL,
        position = c("inputs", "baseCommand", "args", "outputs", "stdout"),
        verbose = FALSE
    ) {
    obj <- if(inherits(sysargs, "SYSargs2")) cmdToCwl(sysargs) else if(inherits(sysargs, "cwlParse")) sysargs else {
        stop("sysargs must be `SYSargs2` or `cwlParse` object")
    }
    if(inherits(obj, "v1")) stop("input sysargs built on syntax version 1, please use `replaceParam`")
    .cmdReplace2(
        obj, x = x,
        index=index,
        position = position,
        verbose = verbose
    )
}

removeParam2 <- function(
        sysargs, index=NULL,
        position = c("inputs", "args", "outputs", "stdout"),
        verbose = FALSE
    ) {
    obj <- if(inherits(sysargs, "SYSargs2")) cmdToCwl(sysargs) else if(inherits(sysargs, "cwlParse")) sysargs else {
        stop("sysargs must be `SYSargs2` or `cwlParse` object")
    }
    if(inherits(obj, "v1")) stop("input sysargs built on syntax version 1, please use `subsetParam`")

    .cmdRemove2(
        obj, index=index,
        position = position,
        verbose = verbose
    )
}


renameParam2 <- function(
        sysargs, index=NULL,
        new_names,
        position = c("inputs", "args", "outputs", "stdout"),
        verbose = FALSE
        ) {
    obj <- if(inherits(sysargs, "SYSargs2")) cmdToCwl(sysargs) else if(inherits(sysargs, "cwlParse")) sysargs else {
        stop("sysargs must be `SYSargs2` or `cwlParse` object")
    }
    if(inherits(obj, "v1")) stop("input sysargs built on syntax version 1, please use `renameParam`")
    .cmdRename2(
        obj, index=index,
        new_names = new_names,
        position = position,
        verbose = verbose
    )
}

############# Internal ###########
## actual cmd
# mycmd \
#     -s sample1.txt \
#     -s sample2.txt \
#     --c \
#     -o myout.txt \
#     a.fasta \
#     --n mouse \
#     > abc.txt

## format ##
# mycmd \                                   -> base command
# p: prefixx;     ;                       \ -> args
# p: prefixx; type; default_value         \ -> inpputs
# name      ; type; default_value         \ -> inpputs no prefix
# p: prefixx; type; out: default_value    \ -> inpputs, outputs
# name      ; type; default_value         \ -> inpputs
# name      ; type; stdout: default_value   -> stdout

# 1. Each line specifies one argument and its default value.
# 2. text before first `;` will be will used as prefix/names. If it starts with
#    keyword "p:", anything after "p:" and before the first `;` will be used as
#    prefix, and the name of this position will be the prefix but with leading
#    dash(s) '-', '--' removed. If there is any duplication, a number index will be
#    added to the end.
# 3. If there is no keyword "p:" before first `;`, all text before first `;` will
#    the name.
# 4. If there is  keyword "p:" before first `;` but nothing before and after the second `;`,
#    this position will be treated as CWL argument instead of input.
# 5. Text between first and second `;` is type. Must be one of  File, Directory, string,
#    int, double, float, long, boolean.
# 6. Text after second `;` and before `\` or end of the line is the default value. If it starts with
#    keywrod "out" or "stdout", this position will also be added to outputs or
#    standard output.
# 7. There is only 1 position with "stdout" allowed and usually it is the last position.
# 8. Ending with "\" is recommended but not required.
##### Test code
# cmd <- '
# mycmd \
#     p: -s; File; sample1.txt \
#     p: -s; File; sample2.txt \
#     p: --c; ; \
#     p: -o; File; out: myout.txt \
#     ref_genome; File; a.fasta \
#     p: --n; string; mouse \
#     mystdout; File; stdout: abc.txt
# '
# library(magrittr)
# cmd <- '
# mycmd \
#     p: -s; File; sample1.txt \
#     p: -s; File; sample2.txt \
#     p: --c; ; \
#     p: -o; File; out: myout.txt \
#     ref_genome; File; a.fasta \
#     p: --n; string; mouse \
#     mystdout; File; stdout: abc.txt
# '
# parsed <- .cmd2cwl2(cmd)
# .catCmdAll2(parsed)
#
# cmd <- ' p: -s; File; sample1.txt'

#' @param cmd_str length 1 commandline string.
#' @param verbose bool
#' @return `cmdParse` object
.cmd2cwl2 <- function(cmd_str, verbose = TRUE) {
    stopifnot(is.character(cmd_str) && length(cmd_str) == 1)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    cmd_lines <- stringr::str_split(cmd_str, "\n", simplify = TRUE) %>% # split args by line
        stringr::str_remove_all("^[ ]+") %>% # remove all leading, ending spaces
        stringr::str_remove_all("[ ]+$")
    # remove empty lines
    cmd_lines <- cmd_lines[cmd_lines != ""]
    cmd_base <- cmd_lines[1]
    cmd_lines <- cmd_lines[-1]
    # split by ; and remove all leading/ending spaces
    cmd_split <- stringr::str_split(cmd_lines, ";") %>%
        lapply(stringr::str_remove_all, "^[ ]{0,}") %>%
        lapply(stringr::str_remove_all, "[ ]{0,}$")
    valid_types <- c("File", "Directory", "string", "int", "double", "float", "long", "boolean")
    cmd_types <- lapply(seq_along(cmd_split), function(x){
        if(length(cmd_split[[x]]) != 3) stop("Each param line needs to have exact 2 ';' to separate prefix/name with type and value. Error in:\n", cmd_lines[[x]])
        has_pre_name <- cmd_split[[x]][1] != ""
        has_type <-  cmd_split[[x]][2] != ""
        has_value <-  cmd_split[[x]][3] != ""
        is_out <- stringr::str_detect(cmd_split[[x]][3], "^out:")
        is_stdout <- stringr::str_detect(cmd_split[[x]][3], "^stdout:")
        if(has_pre_name && !has_type && !has_value) return("arg")
        if(!has_type && has_value) stop(cmd_lines[[x]], ": If default value is provided, type must be provided as well.")
        if(!cmd_split[[x]][2] %in% valid_types) stop("Input, output require a valid type, should be one of:\n", paste(valid_types, collapse = ", "))
        if(is_stdout) {
            if(stringr::str_detect( cmd_split[[x]][1], "^p:")) stop("standard out cannot have prefix")
            return("stdout")
        }
        if(is_out) return("out")
        return("input")
    }) %>% unlist()
    cmd_split <- mapply(x = cmd_split, i = seq_along(cmd_split), FUN = function(x, i) c(x, i), SIMPLIFY = FALSE)
    cmd_input <- cmd_split[cmd_types %in% c("input", "out")]
    cmd_arg <- cmd_split[cmd_types == "arg"]
    cmd_out <- cmd_split[cmd_types == "out"]
    cmd_stdout <- cmd_split[cmd_types == "stdout"]
    if(length(cmd_stdout) > 1) stop("You can only have one standard out")

    # arguments
    args <- lapply(cmd_arg, function(x){
        preF <- if(stringr::str_detect(x[1], "^p:")) stringr::str_remove(x[1], "^p:[ ]{0,}") %>% .removeSpaces2()
        else stop("command with no type and value will be arguments and must have prefix.")
        list(preF = preF, index = x[4])
    })
    if(length(args)) names(args) <- paste0("argument", seq(length(args)))
    # inputs
    inputs <<- lapply(cmd_input, function(x){
        preF <- if(stringr::str_detect(x[1], "^[ ]{0,}p:")) stringr::str_remove(x[1], "^[ ]{0,}p:") %>% .removeSpaces2() else ""
        list(
            preF = preF,
            arg_name =  if(preF != "") stringr::str_remove(preF, "^[-]{0,}") else x[1] %>% .removeSpaces2(),
            type = x[2] %>% .removeSpaces2(),
            value = stringr::str_remove(x[3], "^out:[ ]{0,}") %>% .removeSpaces2(),
            index = x[4]
        )
    })
    if(length(inputs)){
        input_names <- lapply(inputs, `[[`, 'arg_name') %>% unlist()
        inputs <- lapply(inputs, `[`, c('preF', 'type', 'value', 'index'))
        # fix dup names
        input_names <- ave(
            input_names, input_names,
            FUN=function(i) if (length(i) > 1) paste0(i[1], seq_along(i)) else i[1])
        names(inputs) <- input_names
    }
    # outputs
    outputs <- lapply(cmd_out, function(x) list(
        type = x[2] %>% .removeSpaces2(),
        value = stringr::str_remove(x[3], "^out:[ ]{0,}") %>% .removeSpaces2()
    ))
    if(length(outputs)) names(outputs) <- paste0("output", seq(length(outputs)))
    # stdout
    stdout <- if(length(cmd_stdout)) list(list(
        type = cmd_stdout[[1]][2] %>% .removeSpaces2(),
        value = stringr::str_remove(cmd_stdout[[1]][3], "^stdout:[ ]{0,}") %>% .removeSpaces2()
    )) else list()
    if(length(stdout)) names(stdout) <- cmd_stdout[[1]][1]
    commandline <- structure(list(baseCommand = cmd_base, args = args, inputs = inputs, outputs = outputs, stdout = stdout), class = c("list", "cwlParse", "v2"))
    if (verbose) .catCmdAll2(commandline)
    invisible(commandline)
}

## appending methods

##


#' @param cmd `cmdParse` object
#' @param x either a list of values of a param or a length 1 character string
#' which can be translated to the param structure.
#' @param position one of "inputs", "args", "outputs"
#' @param after numeric, where to append the new param, default is after all current params in
#' the `position` selected, `0` means before all current params.
#' @param verbose bool, verbose mode
#' @return `cmdParse` object
#' @details
#' different position has different requirements for the `x`
#' #### List format
#' inputs position require the list contains following items:
#' "name", "preF", "type", "value", "index".
#' args position require the list contains following items:
#' "name", "preF", "index".
#' outputs position require the list contains following items:
#' "name", "type", "value".
#'
#' #### String format
#' Please follow this format, append and replace methods
#' needs 3 semi-colons to separate 4 columns.
#'
#' p: prefixx;     ;                      ; index  -> args
#' p: prefixx; type; default_value        ; index  -> inputs
#' name      ; type; default_value        ; index  -> inputs
#' p: prefixx; type; out: default_value   ; index  -> outputs
#' name      ; type; out: default_value   ; index  -> outputs
.cmdAppend2 <- function(
        cmd, x,
        position = c("inputs", "args", "outputs"),
        after = NULL,
        verbose = FALSE
){
    position <- match.arg(position, c("inputs", "args", "outputs"))
    if (!inherits(cmd, "cwlParse"))  stop("Expecting a 'cwlParse' object")
    if(!is.null(after)) stopifnot(is.numeric(after) && length(after) == 1)
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    # find append index
    after <- if(is.null(after)) length(cmd[[position]]) else {
        if(length(after) > 1) stop("`after` can only be one number")
        as.integer(after)
    }
    # parse string or check list type, names
    if(is.character(x)) {
        x <- .parseInstr2(x, verbose) %>% .checkInstr2(position, cmd, check_dup = TRUE, verbose)
    } else if (is.list(x)) {
        x <- .checkInlist2(x, position, cmd, check_dup = TRUE, verbose)
    } else stop("x can only be a length 1 character string or a list")
    # build required list
    new_list <- .buildInlist2(x, position, verbose)
    # append the list
    if(verbose) message("Appnding to the ", after + 1, "th position in the list")
    cmd[[position]] <- append(cmd[[position]], new_list, after = after)
    if(verbose) .catCmdAll2(cmd)
    if(position == "outputs") message("Note: New output appended. However, outputs are often come in pairs with inputs. ",
                                      "Make sure the corresponding input is there, or append it.")
    invisible(cmd)
}

#' @param cmd `cmdParse` object
#' @param x either a list of values of a param or a length 1 character string
#' which can be translated to the param structure.
#' @param position one of "baseCommand","inputs", "args", "outputs", "stdout"
#' @param index numeric or character, where or the name to replace the new param
#' @param verbose bool, verbose mode
#' @return `cmdParse` object
#' @details
#' different position has different requirements for the `x`
#' #### List format
#' inputs position require the list contains following items:
#' "name", "preF", "type", "value", "index".
#' args position require the list contains following items:
#' "name", "preF", "index".
#' outputs position require the list contains following items:
#' "name", "type", "value".#
#' stdout position require the list contains following items:
#' "name", "type", "value".
#' baseCommand is **not supported** to use list
#'
#' #### String format
#' Please follow this format, append and replace methods
#' needs 3 semi-colons to separate 4 columns.
#' ""
#' p: prefixx;     ;                      ; index  -> args
#' p: prefixx; type; default_value        ; index  -> inputs
#' name      ; type; default_value        ; index  -> inputs
#' p: prefixx; type; out: default_value   ; index  -> outputs
#' name      ; type; out: default_value   ; index  -> outputs
#' name      ; type; stdout: default_value; index  -> stdout
.cmdReplace2 <- function(
        cmd, x, index=NULL,
        position = c("inputs", "baseCommand", "args", "outputs", "stdout"),
        verbose = FALSE
){
    # assertions
    position <- match.arg(position,  c("inputs", "baseCommand", "args", "outputs", "stdout"))
    if(position == "baseCommand") {
        stopifnot(is.character(x) && length(x) == 1)
        cmd$baseCommand <- x
        if(verbose) .catCmdAll2(cmd, base = TRUE, raw_cmd = TRUE, all = FALSE)
        return(cmd)
    }
    if (!inherits(cmd, "cwlParse"))  stop("Expecting a 'cwlParse' object")
    if(!(is.numeric(index) || is.character(index)) || length(index) != 1)
        stop("Index must be a single integer or a character name to indicate the name of param to replace.")
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    # find the position
    if(cmd[[position]][index] %>% names() %>% is.na()) {
        stop("Cannot find index '", index, "' at position ", position)
    }
    # parse string or check list type, names
    if(is.character(x)) {
        x <- .parseInstr2(x, verbose) %>% .checkInstr2(position, cmd, check_dup = TRUE, verbose)
    } else if (is.list(x)) {
        x <- .checkInlist2(x, position, cmd, check_dup = TRUE, verbose)
    } else stop("x can only be a length 1 character string or a list")
    # build required list
    new_list <- .buildInlist2(x, position, verbose)

    if(verbose) message("Replacing to the ", index, " at position ", position)
    pos_names <- names(cmd[[position]])
    after <- which(names(cmd[[position]][index]) %in% pos_names)
    cmd[[position]] <- append(cmd[[position]], new_list, after = after)
    cmd[[position]][index] <-  NULL
    if(verbose) .catCmdAll2(cmd)
    if(position == "outputs") message("Note: New output replaced. However, outputs are often come in pairs with inputs. ",
                                      "Make sure the corresponding input is also updated.")
    invisible(cmd)
}



#' @param cmd `cmdParse` object
#' @param index numeric or character vector, where or the name(s) of param(s) to remove
#' @param position one of "inputs", "args", "outputs", "stdout"
#' @param verbose bool, verbose mode
#' @return `cmdParse` object
.cmdRemove2 <- function(
        cmd, index=NULL,
        position = c("inputs", "args", "outputs", "stdout"),
        verbose = FALSE
){
    # assertions
    position <- match.arg(position,  c("inputs", "args", "outputs", "stdout"))
    if (!inherits(cmd, "cwlParse"))  stop("Expecting a 'cwlParse' object")
    if(!(is.numeric(index) || is.character(index)))
        stop("Index must be integer(s) or a character name(s) to indicate the name of param to remove.")
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    # find the position
    if(cmd[[position]][index] %>% names() %>% is.na() %>% any()) {
        stop("Cannot find at least one of the indices '", index, "' at position ", position)
    }
    # parse string or check list type, names
    if(verbose) message("Remove the ", index, " at position ", position)
    cmd[[position]][index] <-  NULL
    if(verbose) .catCmdAll2(cmd)
    invisible(cmd)
}

#' @param cmd `cmdParse` object
#' @param index numeric or character vector, where or the name(s) of param(s) to remove
#' @param new_names character vector, new names to rename, must be equal length with `index`
#' @param position one of "inputs", "args", "outputs", "stdout"
#' @param verbose bool, verbose mode
#' @return `cmdParse` object
.cmdRename2 <- function(
        cmd, index=NULL,
        new_names,
        position = c("inputs", "args", "outputs", "stdout"),
        verbose = FALSE
){
    # assertions
    position <- match.arg(position,  c("inputs", "args", "outputs", "stdout"))
    if (!inherits(cmd, "cwlParse"))  stop("Expecting a 'cwlParse' object")
    if(length(index) != length(new_names)) stop("new names must be equal length with `index`")
    if(!(is.numeric(index) || is.character(index)))
        stop("Index must be integer(s) or a character name(s) to indicate the name of param to remove.")
    stopifnot(is.logical(verbose) && length(verbose) == 1)
    if(length(unique(new_names)) != length(new_names)) stop("All values in `new_names` must be unique")
    if(any(new_names %in% .getAllCmdNames2(cmd))) stop("At least one of the new names exists in current param names, try different names.")
    # find the position
    if(cmd[[position]][index] %>% names() %>% is.na() %>% any()) {
        stop("Cannot find at least one of the indices '", index, "' at position ", position)
    }
    # parse string or check list type, names
    if(verbose) message("Rename the ", index, " at position ", position)
    pos_names <- names(cmd[[position]])
    index_names <- names(cmd[[position]][index])
    for(i in seq_along(index_names)) {
        pos_names[pos_names == index_names[i]] <- new_names[i]
    }
    names(cmd[[position]]) <- pos_names
    if(verbose) .catCmdAll2(cmd)
    invisible(cmd)
}

## New param list format parser and checker
.buildInlist2 <- function(x, position, verbose = FALSE) {
    return(switch(position,
                  "args" = {
                      cmd_list <- list()
                      cmd_list[[x$name]] <- list(preF = x$preF, index = x$index)
                      cmd_list
                  },
                  "inputs" = {
                      cmd_list <- list()
                      cmd_list[[x$name]] <- list(preF = x$preF, type = x$type, value = x$value, index = x$index)
                      cmd_list
                  },
                  "outputs" = {
                      cmd_list <- list()
                      cmd_list[[x$name]] <- list(type = x$type, value = x$value)
                      cmd_list
                  },
                  "stdout" = {
                      cmd_list <- list()
                      cmd_list[[x$name]] <- list(type = x$type, value = x$value)
                      cmd_list
                  }
    ))
}
.checkInlist2 <- function(x, position, cmd, check_dup = FALSE, verbose = FALSE) {
    if (verbose) message("Detected the input to be a list, checking items in the list...")
    lapply(x[[1]], function(i) {
        if(length(i) != 1 || !is.character(i)) stop("Each item in new param list must be a length 1 character string.")
    })
    cmd_names <- names(x)
    switch(position,
           "args" = {
               if(!all(c("name", "preF", "index") %in% cmd_names))
                   stop('"name", "preF", "index" are required for argument position')
           },
           "inputs" = {
               if(!all(c("name", "preF", "type", "value", "index") %in% cmd_names))
                   stop('"name", "preF", "type", "value", "index" are required for input position')
           },
           "outputs" = {
               if(!all(c("name", "type", "value") %in% cmd_names))
                   stop('"name", "type", "value" are required for output position')
           },
           "stdout" = {
               if(!all(c("name", "type", "value") %in% cmd_names))
                   stop('"name", "type", "value" are required for standard out position')
           }
    )
    if(!check_dup) return()
    pos_names <- names(cmd[[position]])
    if(x$name %in% .getAllCmdNames2(cmd)) {
        new_name <- paste0(c(x$name, "_", sample(c(seq(9), letters), 3)), collapse = "")
        message("Param name ", x$name, " exists, it will be renamed to: ", new_name)
        x$name <- new_name
    }
    x
}

## New param string format parser and checkers
.parseInstr2 <- function(x, verbose=FALSE) {
    cmd_lines <- stringr::str_split(x, "\n", simplify = TRUE) %>% # split args by line
        stringr::str_remove_all("^[ ]+") %>% # remove all leading, ending spaces
        stringr::str_remove_all("[ ]+$")
    if(length(cmd_lines) > 1) stop("Only one new command can be modified or added a time")
    # split by ; and remove all leading/ending spaces
    cmd_split <- stringr::str_split(cmd_lines, ";") %>%
        lapply(stringr::str_remove_all, "^[ ]{0,}") %>%
        lapply(stringr::str_remove_all, "[ ]{0,}$")
    cmd_split <- cmd_split[[1]]

    valid_types <- c("File", "Directory", "string", "int", "double", "float", "long", "boolean")
    cmd_type <- (function(){
        if(length(cmd_split) != 4) stop("New param line needs to have exact 3 ';' to separate prefix/name, type, value, and index. Error in:\n", cmd_lines)
        has_pre_name <- cmd_split[1] != ""
        has_type <-  cmd_split[2] != ""
        has_value <-  cmd_split[3] != ""
        is_out <- stringr::str_detect(cmd_split[3], "^out:")
        is_stdout <- stringr::str_detect(cmd_split[3], "^stdout:")
        if(has_pre_name && !has_type && !has_value) return("arg")
        if(!has_type && has_value) stop(cmd_lines, ": If default value is provided, type must be provided as well.")
        if(!cmd_split[2] %in% valid_types) stop("Input, output require a valid type, should be one of:\n", paste(valid_types, collapse = ", "))
        if(is_stdout) {
            if(stringr::str_detect( cmd_split[1], "^p:")) stop("standard out cannot have prefix")
            return("stdout")
        }
        if(is_out) return("out")
        return("input")
    })()
    if(verbose) message("Detected type is: ", cmd_type)

    switch(cmd_type,
           'arg' = {
               preF <- if(stringr::str_detect(cmd_split[1], "^p:")) stringr::str_remove(cmd_split[1], "^p:[ ]{0,}") %>% .removeSpaces2()
               else stop("New string is detected to be commandline argument type, but it has noe prefix, not allowed")
               arg_name <- if(preF != "") stringr::str_remove(preF, "^[-]{0,}") %>% .removeSpaces2()
               else cmd_split[1] %>% .removeSpaces2()
               list(
                   pos = "arg", name = arg_name, preF = preF,
                   index = if(!is.na(cmd_split[4])) cmd_split[4] %>% .removeSpaces2()
                   else stop("Index required for argument param but not found in the string.")
               )
           },
           'input' = {
               preF <- if(stringr::str_detect(cmd_split[1], "^p:")) stringr::str_remove(cmd_split[1], "^p:[ ]{0,}") %>% .removeSpaces2() else ""
               arg_name <- if(preF != "") stringr::str_remove(preF, "^[-]{0,}") %>% .removeSpaces2()
               else cmd_split[1] %>% .removeSpaces2()
               list(
                   pos = "input", name = arg_name,
                   preF = preF, type = cmd_split[2] %>% .removeSpaces2(),
                   index = if(!is.na(cmd_split[4])) cmd_split[4] %>% .removeSpaces2()
                   else stop("Index required for input param but not found in the string."),
                   value = stringr::str_remove(cmd_split[3], "^out:[ ]{0,}") %>% .removeSpaces2()
               )
           },
           'out' = list(
               pos = "out", name = if(!is.na(cmd_split[1])) cmd_split[1] %>% .removeSpaces2() else "",
               type = cmd_split[2] %>% .removeSpaces2(),
               value = stringr::str_remove(cmd_split[3], "^out:[ ]{0,}") %>% .removeSpaces2()
           ),
           'stdout' = list(
               pos = "stdout", name = if(!is.na(cmd_split[1])) cmd_split[1] %>% .removeSpaces2() else "",
               type = cmd_split[2] %>% .removeSpaces2(),
               value = stringr::str_remove(cmd_split[3], "^stdout:[ ]{0,}") %>% .removeSpaces2()
           )
    )
}
.checkInstr2 <- function(x, position, cmd, check_dup = FALSE, verbose = FALSE) {
    cmd_pos <- x$pos
    position <- switch(position,
                       "inputs" = "input",
                       "args" = "arg",
                       "outputs" = "out",
                       "stdout" = "stdout",
                       position
    )
    if(cmd_pos != position) stop("The parsed param string is good for ",
                                 cmd_pos,
                                 " position, but the position you have provided is ",
                                 position)
    if(!check_dup) return()
    pos_names <- names(cmd[[position]])
    if(x$name %in% .getAllCmdNames2(cmd)) {
        new_name <- paste0(c(x$name, "_", sample(c(seq(9), letters), 3)), collapse = "")
        message("Param name ", x$name, " exists, it will be renamed to: ", new_name)
        x$name <- new_name
    }
    x
}

## printing methods----
#' @param cmd `cmdParse` object
#' @param base bool, base command
#' @param args bool, arguments
#' @param inputs bool, inputs
#' @param outputs bool, outputs
#' @param stdout bool, standard out
#' @param raw_cmd bool, raw parsed command
#' @param all bool, print all positions?
#'
#' @return invisible of `cmdParse` object
.catCmdAll2 <- function(cmd, base = FALSE, args = FALSE, inputs = FALSE,
                       outputs = FALSE, stdout = FALSE, raw_cmd = FALSE, all = TRUE) {
    if(all) base <- args <- inputs <- outputs <- stdout <- raw_cmd <- TRUE
    if(base) .catBase2(cmd$baseCommand)
    if(args) .catArgs2(cmd$args)
    if(inputs) .catInputs2(cmd$inputs)
    if(outputs) .catOutputs2(cmd$outputs)
    if(stdout) .catStdout2(cmd$stdout)
    if(raw_cmd) .catCmdRaw2(cmd)
    return(invisible(cmd))
}

.catBase2 <- function(baseCommand) {
    cat(crayon::bgBlue("*****BaseCommand*****\n"))
    cat(baseCommand, "\n")
}

.catArgs2 <- function(args) {
    cat(crayon::bgBlue("*****Arguments*****\n"))
    if(!length(args)) return(invisible())
    for (i in seq_along(args)) {
        cat(
            names(args)[i], ":\n",
            "    prefix: ", args[[i]][["preF"]], "\n",
            "    position: ", args[[i]][["index"]], "\n",
            sep = ""
        )
    }
}

.catInputs2 <- function(inputs) {
    cat(crayon::bgBlue("*****Inputs*****\n"))
    if(!length(inputs)) return(invisible())
    input_names <- names(inputs)
    for (i in seq_along(inputs)) {
        cat(
            input_names[i], ":\n",
            "    type: ", inputs[[i]][["type"]], "\n",
            "    prefix: ", inputs[[i]][["preF"]], "\n",
            "    default value: ", inputs[[i]][["value"]], "\n",
            "    position: ", inputs[[i]][["index"]], "\n",
            sep = ""
        )
    }
}

.catOutputs2 <- function(outputs) {
    cat(crayon::bgBlue("*****Outputs*****\n"))
    if(!length(outputs)) return(invisible())
    out_names <- names(outputs)
    for (i in seq_along(outputs)) {
        cat(
            out_names[i], ":\n",
            "    type: ", outputs[[i]][["type"]], "\n",
            "    default value: ", outputs[[i]][["value"]], "\n",
            sep = ""
        )
    }
}

.catStdout2 <- function(stdout) {
    cat(crayon::bgBlue("*****Standard Outputs*****\n"))
    if(!length(stdout)) return(invisible())
    cat(
        names(stdout)[1], ":\n",
        "    type: ", stdout[[1]][["type"]], "\n",
        "    default value: ", stdout[[1]][["value"]], "\n",
        sep = ""
    )
}

.catCmdRaw2 <- function(commandline) {
    arg_in_list <- append(commandline$inputs,  commandline$args)
    arg_in_index <- lapply(arg_in_list, `[[`, "index") %>% unlist() %>% sort()
    extract_order <- names(arg_in_index)
    arg_in_ordered <- arg_in_list[extract_order]
    cmd_string <- lapply(arg_in_ordered, `[`, c("preF", "value")) %>%
        unlist() %>%
        paste0(collapse = " ")
    stdout <- if(length(commandline$stdout)) paste(">", commandline$stdout[[1]]$value) else ""
    cat(crayon::bgBlue("*****Parsed raw command line*****\n"))
    cat(commandline$baseCommand, cmd_string, stdout, "\n")
}

.removeSpaces2 <- function(x){
    x %>% stringr::str_remove("^[ ]{0,}") %>% stringr::str_remove("[ ]{0,}$")
}

.getAllCmdNames2 <- function(cmd) {
    lapply(cmd, names) %>% unlist()
}

