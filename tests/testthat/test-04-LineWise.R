## It will check the LineWise class and methods
library(systemPipeR)
library(systemPipeRdata)

genWorkenvir("new", mydirname = file.path(tempdir(), "newtest"))
setwd(file.path(tempdir(), "newtest"))

test_that("check_LineWise_test", {
    ## build simple example
     xx <- LineWise({
      x <- 1+1; y <- 100
    })
    ## build instance from Rmd
    rmdPath <- system.file("extdata", "spr_simple_wf.Rmd", package="systemPipeR")
    sal <- SPRproject(overwrite = TRUE)
    sal <- importWF(sal, rmdPath)
    lw <- stepsWF(sal[2])[[1]]
    show(lw)
    ## Class
    expect_s4_class(sal, "SYSargsList")
    expect_length(sal, 5)
    expect_s4_class(lw, "LineWise")
    expect_length(lw, 1)
    expect_length(lw$codeLine, 1)
    ## Methods
     expect_type(names(lw), "character")
     expect_type(length(lw), "integer")
     expect_type(lw$codeLine, "expression")
     expect_length(codeLine(lw), 0)
     expect_type(codeChunkStart(lw), "integer")
     expect_type(rmdPath(lw), "character")
     expect_type(dependency(lw), "list")
     ## Coerce
     lw2 <- linewise(lw)
     expect_type(lw2, "list")
     ## Constructor
     expect_s4_class(LineWise({1+1}), "LineWise")
     ## Subset 
     expect_length(lw[2], 1)
     lw <- stepsWF(sal[5])[[1]]
     expect_length(lw[-2], 4)
    ## Replacements LineWise
    replaceCodeLine(lw, line = 2) <- "5+5"
    expect_equal(as.character(lw[2]$codeLine), "5 + 5")
    appendCodeLine(lw) <- "6+7"
    expect_equal(as.character(lw[6]$codeLine), "6 + 7")
   ## Error
    expect_error(codeLine(sal, step=6))
    ## Replacements SYSargsList
    replaceCodeLine(sal, step=1, line=1) <- LineWise(code={
                                                        5+5
                                                      })
    codeLine(sal[1])
    expect_equal(as.character(stepsWF(sal[1])[[1]]$codeLine[1]), "5 + 5")
    appendCodeLine(sal, step=1) <- "66+55"
    codeLine(sal[1])
    expect_equal(as.character(stepsWF(sal[1])[[1]]$codeLine[2]), "66 + 55")
   #runWF()
   sal <- runWF(sal, steps = 1:2)
   expect_setequal(sal$statusWF[[1]]$status.summary, "Success")
   ## Methods
   expect_setequal(status(stepsWF(sal)[[1]])$status.summary, "Success")
   expect_type(files(stepsWF(sal)[[1]])[[1]], "character")
   ## replacement methods
   expect_error(appendStep(sal) <- sal[1])
   renameStep(sal, 2) <- "newStep"
   ## `+` method
   sal <- sal[1] + sal[2]
   expect_length(sal, 2)
})
