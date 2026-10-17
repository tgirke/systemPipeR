## Tests that run external command-line tools (hisat2, bowtie2, trimmomatic,
## macs2, ...) are off by default, so R CMD check passes on builders that
## lack these tools (e.g. r-universe). Note: r-universe unsets CI, so
## skip_on_ci() does not skip there. To run them locally:
##   Sys.setenv(SPR_TEST_CMDLINE = "true")    # within R, or
##   SPR_TEST_CMDLINE=true R CMD check ...    # from the shell
skip_if_no_cmdline <- function(tools) {
    skip_if_not(isTRUE(as.logical(Sys.getenv("SPR_TEST_CMDLINE", "false"))),
                "command-line tool tests off (set SPR_TEST_CMDLINE=true)")
    missing <- tools[!nzchar(Sys.which(tools))]
    skip_if(length(missing) > 0,
            paste("not in PATH:", paste(missing, collapse = ", ")))
}
