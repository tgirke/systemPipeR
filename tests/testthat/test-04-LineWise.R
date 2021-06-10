## It will check the SYSargs2 class and methods
library(systemPipeR)
# library(systemPipeRdata)
setwd(file.path(tempdir(), "rnaseq"))

test_that("check_LineWise_test", {
    ## build instance 
    rmdPath <- system.file("extdata/", "systemPipeTEST.Rmd", package="systemPipeR")
    sal <- importWF(rmdPath, overwrite=TRUE)
    lw <- stepsWF(sal[3])[[1]]
    show(lw)
    length(lw)
    names(lw)
    codeLine(lw)
    codeChunkStart(lw)
    rmdPath(lw)
    l <- lw[2]
    l_sub <- lw[-2]
    codeLine(l)
    codeLine(l_sub)
    ## Methods
     expect_type(names(lw), "character")
     expect_type(length(lw), "integer")
     expect_type(lw$codeLine, "expression")

     
    replaceCodeLine(lw, line=2) <- "5+5"
    expect_length(lw, 3)
    codeLine(lw)
    appendCodeLine(lw) <- "6+7"
    codeLine(lw)
    expect_length(lw, 4)
   #  
    expect_error(codeLine(sal[4]))
    
    replaceCodeLine(sal, step=1, line=2) <- "5+5"
    codeLine(sal[1])

    appendCodeLine(sal, step=1) <- "66+55"
    codeLine(sal[1])

    appendCodeLine(sal, step=10, after=1) <- "66+55"
    codeLine(sal[10])

   #  ## Class
    expect_s4_class(sal, "SYSargsList")
    expect_s4_class(lw, "LineWise")
    expect_length(sal, 14)
    expect_length(lw, 4)
    expect_length(stepsWF(sal[1])[[1]], 3)
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
