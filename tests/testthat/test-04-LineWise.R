## It will check the LineWise class and methods
library(systemPipeR)
library(systemPipeRdata)

genWorkenvir("new", mydirname = file.path(tempdir(), "newtest"))
setwd(file.path(tempdir(), "newtest"))

test_that("check_LineWise_test", {
    ## build instance 
    rmdPath <- system.file("extdata", "systemPipeTEST.Rmd", package="systemPipeR")
    sal <- SPRproject(overwrite = TRUE)
    sal <- importWF(sal, rmdPath, overwrite=TRUE)
    lw <- stepsWF(sal[3])[[1]]
    show(lw)
    ## Class
    expect_s4_class(sal, "SYSargsList")
    expect_length(sal, 14)
    expect_s4_class(lw, "LineWise")
    expect_length(lw, 3)
    expect_length(lw$codeLine, 3)
    ## Methods
     expect_type(names(lw), "character")
     expect_type(length(lw), "integer")
     expect_type(lw$codeLine, "expression")
     expect_length(codeLine(lw), 0)
     expect_type(codeChunkStart(lw), "integer")
     #expect_type(rmdPath(lw), "character")
     ## Coerce
     lw2 <- linewise(lw)
     expect_type(lw2, "list")
     ## Constructor
     expect_s4_class(LineWise({1+1}), "LineWise")
     ## Subset 
     expect_length(lw[2], 1)
     expect_length(lw[-2], 2)
    ## Replacements LineWise
    replaceCodeLine(lw, line=2) <- "5+5"
    expect_equal(as.character(lw[2]$codeLine), "5 + 5")
    appendCodeLine(lw) <- "6+7"
    expect_equal(as.character(lw[4]$codeLine), "6 + 7")
   ## Error
    expect_error(codeLine(sal[4]))
    ## Replacements SYSargsList
    replaceCodeLine(sal, step=1, line=1) <- "5+5"
    codeLine(sal[1])
    expect_equal(as.character(stepsWF(sal[1])[[1]]$codeLine[1]), "5 + 5")
    appendCodeLine(sal, step=1) <- "66+55"
    codeLine(sal[1])
    expect_equal(as.character(stepsWF(sal[1])[[1]]$codeLine[2]), "66 + 55")
   #  ##runWF()
   # # sal <- runWF(sal)
   # # expect_setequal(statusWF(sal)[[1]][[1]], "DONE")
   # # check <- check.output(sal)
   # # expect_setequal(check$Step_1$Existing_Files, "1")
   # # 
   # ## replacement methods
   # expect_warning(appendStep(sal) <- sal)
   # renameStep(sal, 2) <- "newStep"
   # ## `+` method
   # sal <- sal[1] + sal[2] 
   # expect_length(sal, 2)
})
